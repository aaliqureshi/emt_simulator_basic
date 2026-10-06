using Test
using LinearAlgebra
using SparseArrays
using ForwardDiff
using NonlinearSolve

# Load the package normally so dependency declarations are tested too.
using Barq
using Barq.PowerFlow: PowerFlowCache, powerflow!

function test_system(; sparse_y = false, pq = true)
    bus = Bus{Float64}(3)
    bus.idx .= 1:3
    bus.v .= [1.04, 1.02, 1.0]
    bus.theta .= [0.03, 0.0, 0.0]
    Y = ComplexF64[2-8im -1+4im -1+4im; -1+4im 2-8im -1+4im; -1+4im -1+4im 2-8im]
    models = (bus = bus, generator = (bus = [2], p_m = [0.4]),
              load = (bus = [3], p = [0.6], q = [0.2]))
    return (; models, Y = sparse_y ? sparse(Y) : Y,
            non_slack_buses = [2, 3], v_update_buses = pq ? [3] : Int[])
end

initial_state(sys) = vcat(sys.models.bus.v[sys.v_update_buses],
                          sys.models.bus.theta[sys.non_slack_buses])

# Independent, allocating formulation of the power-balance equations.
function reference_powerflow(u, sys)
    V = eltype(u).(sys.models.bus.v)
    theta = eltype(u).(sys.models.bus.theta)
    n_v = length(sys.v_update_buses)
    V[sys.v_update_buses] = u[1:n_v]
    theta[sys.non_slack_buses] = u[n_v+1:end]
    v = V .* cis.(theta)
    S = v .* conj.(sys.Y * v)
    P, Q = real.(S), imag.(S)
    for (b, p) in zip(sys.models.generator.bus, sys.models.generator.p_m)
        P[b] -= p
    end
    for (b, p, q) in zip(sys.models.load.bus, sys.models.load.p, sys.models.load.q)
        P[b] += p
        Q[b] += q
    end
    return vcat(Q[sys.v_update_buses], P[sys.non_slack_buses])
end

function residual_allocations(du, u, cache)
    powerflow!(du, u, cache)
    return @allocated powerflow!(du, u, cache)
end

@testset "Power-flow cache" begin
    for sparse_y in (false, true), pq in (false, true)
        @testset "sparse=$sparse_y, PQ=$pq" begin
            sys = test_system(; sparse_y, pq)
            u = initial_state(sys)
            cache = PowerFlowCache(sys, u)
            du = similar(u)
            original_v, original_theta = copy(sys.models.bus.v), copy(sys.models.bus.theta)

            @test powerflow!(du, u, cache) === nothing
            @test du ≈ reference_powerflow(u, sys)
            @test residual_allocations(du, u, cache) == 0
            u_changed = u .+ 0.025
            powerflow!(du, u_changed, cache)
            @test du ≈ reference_powerflow(u_changed, sys)
            powerflow!(du, u, cache)
            @test du ≈ reference_powerflow(u, sys)
            @test sys.models.bus.v == original_v
            @test sys.models.bus.theta == original_theta

            f! = (out, x) -> powerflow!(out, x, cache)
            config = ForwardDiff.JacobianConfig(f!, du, u)
            J = zeros(length(u), length(u))
            ForwardDiff.jacobian!(J, f!, du, u, config)
            @test J ≈ ForwardDiff.jacobian(x -> reference_powerflow(x, sys), u)
            # Exercise the same cached dual buffers used during Jacobian evaluation.
            @test residual_allocations(similar(config.duals[2]), config.duals[2], cache) == 0
            powerflow!(du, u, cache)
            @test du ≈ reference_powerflow(u, sys)

            reference = solve(NonlinearProblem((x, p) -> reference_powerflow(x, p), u, sys),
                              NewtonRaphson(); abstol=1e-9, reltol=1e-9, maxiters=100)
            sol = solve_power_flow!(sys)
            @test sol.u ≈ reference.u atol=1e-9
            @test norm(reference_powerflow(sol.u, sys), Inf) < 1e-9
            @test sys.models.bus.v[sys.v_update_buses] == sol.u[cache.balance_reactive_power]
            @test sys.models.bus.theta[sys.non_slack_buses] == sol.u[cache.balance_real_power]
            @test sys.models.bus.v[1] == original_v[1]
            @test sys.models.bus.theta[1] == original_theta[1]
            @test sys.models.bus.vd ≈ sys.models.bus.v .* cos.(sys.models.bus.theta)
            @test sys.models.bus.vq ≈ sys.models.bus.v .* sin.(sys.models.bus.theta)
        end
    end

    @testset "Multiple devices at one bus" begin
        sys = test_system()
        models = merge(sys.models, (generator = (bus = [2, 2], p_m = [0.1, 0.3]),
                                   load = (bus = [3, 3], p = [0.2, 0.4], q = [0.1, 0.1])))
        repeated = merge(sys, (; models))
        u = initial_state(sys)
        du = similar(u)
        powerflow!(du, u, PowerFlowCache(repeated, u))
        @test du ≈ reference_powerflow(u, sys)
    end
end

@testset "Repository power-flow cases" begin
    for filename in ("ieee14_fault_barq.xlsx", "ieee39_fault.xlsx")
        @testset "$filename" begin
            models = load_data(joinpath(@__DIR__, "..", "cases", "Fault_Cases", filename))
            sys = build_system(models)
            u = initial_state(sys)
            cache = PowerFlowCache(sys, u)
            du = similar(u)
            powerflow!(du, u, cache)
            @test du ≈ reference_powerflow(u, sys)
            @test residual_allocations(du, u, cache) == 0

            f! = (out, x) -> powerflow!(out, x, cache)
            config = ForwardDiff.JacobianConfig(f!, du, u)
            J = zeros(length(u), length(u))
            ForwardDiff.jacobian!(J, f!, du, u, config)
            @test J ≈ ForwardDiff.jacobian(x -> reference_powerflow(x, sys), u)
            @test residual_allocations(similar(config.duals[2]), config.duals[2], cache) == 0

            reference = solve(NonlinearProblem((x, p) -> reference_powerflow(x, p), u, sys),
                              NewtonRaphson(); abstol=1e-9, reltol=1e-9, maxiters=100)
            sol = solve_power_flow!(sys)
            @test sol.u ≈ reference.u atol=1e-9
            @test norm(reference_powerflow(sol.u, sys), Inf) < 1e-9

            # Check the solver's actual differentiation configuration, since
            # testing powerflow! directly misses AutoSpecialize's chunk-size-1 fallback.
            solver_cache = init(sol.prob, NewtonRaphson())
            solver_config = solver_cache.jac_cache.di_extras.config
            @test ForwardDiff.chunksize(solver_config) == ForwardDiff.pickchunksize(length(u))
        end
    end
end
