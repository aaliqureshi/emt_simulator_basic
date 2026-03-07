module PowerFlow

export solve_power_flow!

include("Models.jl"); using .Models

using NonlinearSolve

function power_flow_traditional!(du, u, p)
    T = eltype(u)

    address, models, Y, non_slack_buses, v_update_buses = p

    bus = models.bus
    generator = models.generator
    load = models.load

    V = Vector{T}(bus.v)
    theta = Vector{T}(bus.theta)

    V[v_update_buses] = @. u[address["balance_reactive_power"]]
    theta[non_slack_buses] = @. u[address["balance_real_power"]]

    P = zeros(T, length(bus.idx))
    Q = zeros(T, length(bus.idx))

    v = V .* cis.(theta)
    S = v .* conj.(Y * v)

    P = real.(S)
    Q = imag.(S)

    ## power balance : (power injection into network) + (load) - (generation) = 0
    P[generator.bus] -= @. generator.p_m
    P[load.bus] += @. load.p
    Q[load.bus] += @. load.q

    du[address["balance_reactive_power"]] .= @. Q[v_update_buses]
    du[address["balance_real_power"]] .= @. P[non_slack_buses]
end

function powerflow(u, p)
    T = eltype(u)
    du = zeros(T, length(u))
    power_flow_traditional!(du, u, p)
    return du
end

function solve_power_flow!(sys)
    # Assumes bus numbering is 1:n_bus (sys.Y and bus vectors indexed by bus number).
    models = sys.models
    Y = sys.Y
    non_slack_buses = sys.non_slack_buses
    v_update_buses = sys.v_update_buses

    # build address dict for power flow variables
    idx_balance_reactive_power = 1 : length(v_update_buses)
    idx_balance_real_power = idx_balance_reactive_power[end]+1 : idx_balance_reactive_power[end]+(length(non_slack_buses))

    address = Dict(
        "balance_reactive_power" => idx_balance_reactive_power,
        "balance_real_power" => idx_balance_real_power
    )

    # initial guess from bus data
    u0 = vcat(models.bus.v[v_update_buses], models.bus.theta[non_slack_buses])

    p = (address, models, Y, non_slack_buses, v_update_buses)

    prob = NonlinearProblem(powerflow, u0, p)
    sol = solve(prob, NewtonRaphson(); abstol = 1e-9, reltol = 1e-9, maxiters = 100)

    # populate models with the solution
    models.bus.v[v_update_buses] .= sol.u[address["balance_reactive_power"]]
    models.bus.theta[non_slack_buses] .= sol.u[address["balance_real_power"]]

    # convert to d-q components
    phasor2DP!(models.bus)

    return sol
end

end # module
