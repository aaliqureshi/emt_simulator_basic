# Also runs the existing residual/cache and direct-solver regression checks.
include("test_power_flow.jl")
using Barq.PowerFlow: _powerflow_jacobian

@testset "Sparse power-flow Jacobian" begin
    for pq in (false, true)
        sys = test_system(; sparse_y=true, pq)
        # Exercise shunts and an asymmetric admittance, as well as non-flat angles.
        sys.Y[1, 2] *= cis(0.1)
        sys.Y[2, 2] += 0.03im
        u = initial_state(sys) .+ 0.02
        J = _powerflow_jacobian(u, PowerFlowCache(sys, u))
        @test issparse(J)
        @test Matrix(J) ≈ ForwardDiff.jacobian(x -> reference_powerflow(x, sys), u)
    end
end

@testset "Flat-start continuation" begin
    for pq in (false, true)
        sys = test_system(; sparse_y=true, pq)
        original_v, original_theta = copy(sys.models.bus.v), copy(sys.models.bus.theta)
        original_p, original_q = copy(sys.models.load.p), copy(sys.models.load.q)
        reference = solve(NonlinearProblem((x, p) -> reference_powerflow(x, p),
                                           initial_state(sys), sys), NewtonRaphson();
                          abstol=1e-10, reltol=1e-10)
        sol = solve_power_flow_continuation!(sys; verbose=false)
        @test sol.converged
        @test sol.retcode == :Success
        @test sol.lambda == 1.0
        @test sol.residual_norm <= 1e-9
        @test norm(reference_powerflow(sol.u, sys), Inf) <= 1e-9
        @test sol.u ≈ reference.u atol=1e-9
        @test sys.models.bus.v[setdiff(1:3, sys.v_update_buses)] == original_v[setdiff(1:3, sys.v_update_buses)]
        @test sys.models.bus.theta[1] == original_theta[1]
        @test sys.models.load.p == original_p
        @test sys.models.load.q == original_q
        @test sys.models.bus.vd ≈ sys.models.bus.v .* cos.(sys.models.bus.theta)
        @test sys.models.bus.vq ≈ sys.models.bus.v .* sin.(sys.models.bus.theta)

        # Arbitrary unknown guesses in the bus table must not affect the path.
        other = test_system(; sparse_y=true, pq)
        other.models.bus.v[other.v_update_buses] .= 0.4
        other.models.bus.theta[other.non_slack_buses] .= 2.0
        other_sol = solve_power_flow_continuation!(other; verbose=false)
        @test other_sol.u ≈ sol.u atol=1e-12
        @test other_sol.history == sol.history
    end

    @testset "Failure leaves system unchanged" begin
        for kwargs in ((; maxsteps=1), (; maxiters=0, initial_step=0.2, min_step=0.1))
            sys = test_system()
            phasor2DP!(sys.models.bus)
            before = deepcopy(sys.models)
            sol = solve_power_flow_continuation!(sys; kwargs..., verbose=false)
            @test !sol.converged
            @test sol.lambda < 1
            @test sol.retcode in (:MaxSteps, :StepSizeTooSmall)
            @test sol.residual_norm ≈ norm(reference_powerflow(sol.u, sys), Inf)
            @test sys.models.bus.v == before.bus.v
            @test sys.models.bus.theta == before.bus.theta
            @test sys.models.bus.vd == before.bus.vd
            @test sys.models.bus.vq == before.bus.vq
            if haskey(kwargs, :maxiters)
                @test !any(h.accepted for h in sol.history)
                @test sol.history[2].step == sol.history[1].step / 2
            end
        end
    end

    @testset "Invalid options" begin
        sys = test_system()
        @test_throws ArgumentError solve_power_flow_continuation!(sys; initial_step=0)
        @test_throws ArgumentError solve_power_flow_continuation!(sys; min_step=0.1)
        @test_throws ArgumentError solve_power_flow_continuation!(sys; max_step=Inf)
        @test_throws ArgumentError solve_power_flow_continuation!(sys; abstol=NaN)
        @test_throws ArgumentError solve_power_flow_continuation!(sys; maxiters=-1)
        @test_throws ArgumentError solve_power_flow_continuation!(sys; maxsteps=0)
    end

    for filename in ("ieee14_fault_barq.xlsx", "ieee39_fault.xlsx")
        sys = build_system(load_data(joinpath(@__DIR__, "..", "cases", "Fault_Cases", filename)))
        sol = solve_power_flow_continuation!(sys; verbose=false)
        @test sol.converged
        @test norm(reference_powerflow(sol.u, sys), Inf) <= 1e-9
    end
end
