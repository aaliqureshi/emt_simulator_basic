module PowerFlow

export PowerFlowCache, powerflow!, solve_power_flow!, solve_power_flow_continuation!

using ..Models

using ForwardDiff: pickchunksize
using LinearAlgebra: mul!
using NonlinearSolve
using PreallocationTools: DiffCache, get_tmp
using SparseArrays

"""
    PowerFlowCache(sys, u0)

Workspace for `powerflow!`, initialized once for a system and its initial guess.
The unknowns are PQ voltage magnitudes followed by non-slack voltage angles.
DiffCache buffers support both ordinary residuals and ForwardDiff Jacobians.
Each concurrent solve needs its own cache.
"""
struct PowerFlowCache{M, YT, BI, VI, DC}
    models::M
    Y::YT
    non_slack_buses::BI
    v_update_buses::VI
    balance_reactive_power::UnitRange{Int}
    balance_real_power::UnitRange{Int}
    V::DC
    theta::DC
    P::DC
    Q::DC
    v_storage::DC
    current_storage::DC
end

function PowerFlowCache(sys, u0)
    n_bus = length(sys.models.bus.idx)
    n_v = length(sys.v_update_buses)
    n_theta = length(sys.non_slack_buses)
    length(u0) == n_v + n_theta || throw(DimensionMismatch("invalid power-flow state length"))
    T = eltype(u0)
    chunk_size = pickchunksize(max(1, length(u0)))
    buffer(n) = DiffCache(zeros(T, n), chunk_size)
    Y = sparse(sys.Y)


    return PowerFlowCache(
        sys.models, Y, sys.non_slack_buses, sys.v_update_buses,
        1:n_v, (n_v + 1):(n_v + n_theta),
        buffer(n_bus), buffer(n_bus), buffer(n_bus), buffer(n_bus),
        buffer(2n_bus), buffer(2n_bus),
    )
end

function powerflow!(du, u, cache::PowerFlowCache)
    bus = cache.models.bus
    generator = cache.models.generator
    load = cache.models.load

    V = get_tmp(cache.V, u)
    theta = get_tmp(cache.theta, u)
    P = get_tmp(cache.P, u)
    Q = get_tmp(cache.Q, u)
    # Real backing storage lets DiffCache select Dual buffers before we view
    # each real/imaginary pair as a complex number, without copying.
    v = reinterpret(Complex{eltype(u)}, get_tmp(cache.v_storage, u))
    current = reinterpret(Complex{eltype(u)}, get_tmp(cache.current_storage, u))

    copyto!(V, bus.v)
    copyto!(theta, bus.theta)
    for (i, b) in zip(cache.balance_reactive_power, cache.v_update_buses)
        V[b] = u[i]
    end
    for (i, b) in zip(cache.balance_real_power, cache.non_slack_buses)
        theta[b] = u[i]
    end

    @. v = V * cis(theta)
    mul!(current, cache.Y, v)
    for i in eachindex(v)
        S = v[i] * conj(current[i])
        P[i] = real(S)
        Q[i] = imag(S)
    end

    # Power balance: network injection + load - generation = 0.
    for i in eachindex(generator.bus)
        P[generator.bus[i]] -= generator.p_m[i]
    end
    for i in eachindex(load.bus)
        P[load.bus[i]] += load.p[i]
        Q[load.bus[i]] += load.q[i]
    end

    for (i, b) in zip(cache.balance_reactive_power, cache.v_update_buses)
        du[i] = Q[b]
    end
    for (i, b) in zip(cache.balance_real_power, cache.non_slack_buses)
        du[i] = P[b]
    end
    return nothing
end

function solve_power_flow!(sys)
    # Assumes bus numbering is 1:n_bus (sys.Y and bus vectors indexed by bus number).
    models = sys.models
    non_slack_buses = sys.non_slack_buses
    v_update_buses = sys.v_update_buses

    # initial guess from bus data
    u0 = vcat(models.bus.v[v_update_buses], models.bus.theta[non_slack_buses])
    # u0 = vcat(ones(length(v_update_buses)), zeros(length(non_slack_buses)))
    cache = PowerFlowCache(sys, u0)

    # FullSpecialize preserves ForwardDiff's normal derivative batching. The
    # default AutoSpecialize wrapper restricts ForwardDiff to chunk size 1,
    # making large Jacobians slower despite the allocation-free residual.
    f = NonlinearFunction{true, NonlinearSolve.SciMLBase.FullSpecialize}(powerflow!)
    prob = NonlinearProblem{true}(f, u0, cache)
    stats = @timed solve(prob, NewtonRaphson(); abstol = 1e-9, reltol = 1e-9, maxiters = 100)

    sol = stats.value
    if Int(sol.retcode) == 1
        println("power flow solved in $(round(stats.time/1e-3, digits = 4)) ms")
    else
        println("[!] power flow failed with retcode: $(sol.retcode)")
        # return nothing
        return sol
    end
    # populate models with the solution
    @views models.bus.v[v_update_buses] .= sol.u[cache.balance_reactive_power]
    @views models.bus.theta[non_slack_buses] .= sol.u[cache.balance_real_power]

    # convert to d-q components
    phasor2DP!(models.bus)

    # return nothing
    return sol
end

include("power_flow_continuation.jl")

end # module
