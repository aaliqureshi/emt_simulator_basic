#=
Run with: julia --project=. scripts/test_openipsl_fault.jl

Compare unmodified loads and both OpenIPSL low-voltage characteristics for the
IEEE 39-bus fault at bus 20, Xf = 0.015 p.u. Only the load law changes between
each paired run; line currents remain differential and no re-initialization is
performed. Each fault-on run attempts exactly the first post-event step.
=#
using Barq, MyDiffEq, ForwardDiff, LinearAlgebra, CSV, DataFrames, Test

function rhs_with_load!(du, u, p, t; load_options...)
    solve_generator!(du, u, p)
    solve_line!(du, u, p, t)
    solve_fault!(du, u, p, t)
    balance!(du, u, p; load_options...)
end

function first_step_residual(rhs!, u, previous, p, mass, h, method)
    f = similar(u)
    rhs!(f, u, p, h)
    r = mass * (u - previous) - h*f
    if method == :Trap
        fprevious = similar(u)
        rhs!(fprevious, previous, p, 0.0)
        r = mass * (u - previous) - (h/2)*(f + fprevious)
    end
    algebraic = findall(iszero, diag(mass))
    r[algebraic] = f[algebraic]
    return r
end

function main()
    output_dir = joinpath(@__DIR__, "..", "outputs", "openipsl_fault_comparison")
    mkpath(output_dir)
    models = load_data("cases/Fault_Cases/ieee39_fault.xlsx")
    models.fault.bus = [20]
    models.fault.x_fault[1] = 0.015
    sys = build_system(models)
    solve_power_flow!(sys)
    run_static_init!(sys)
    address = build_dynamic_address(sys)
    mass = build_mass_matrix(sys, address)
    initial = build_initial_conditions(sys, address)
    ppre = (address, models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, 0.0)
    ppost = (address, models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, 1.0)
    @assert all(>(0), diag(mass)[address["line_id"]])
    @assert all(>(0), diag(mass)[address["line_iq"]])

    # Compare the disabled-conversion path to the exact old implementation
    # preserved in bus.jl, including at voltages below both breakpoints.
    bus_source = read(joinpath(@__DIR__, "..", "src", "models", "bus.jl"), String)
    legacy = split(split(bus_source, "#= Previous balance! implementation, preserved for comparison.\n"; limit=2)[2], "=#"; limit=2)[1]
    legacy_module = Module(:LegacyBalanceForComparison)
    Base.include_string(legacy_module, legacy)
    old_balance! = Base.invokelatest(getfield, legacy_module, :balance!)
    factors = Barq.Models.BusModel._openipsl_load_factors
    @testset "OpenIPSL load integration" begin
        @test factors(0.175, 0.7, 1) == (0.125, 1.0)
        @test factors(0.525, 0.7, 1)[1] ≈ 0.875
        @test factors(1.0, 0.7, 1) == (1.0, 1.0)
        # Preserve upstream's exact-boundary fall-through, rather than silently
        # changing it in a comparison intended to reproduce that implementation.
        @test factors(0.35, 0.7, 1) == (1.0, 1.0)
        @test factors(0.0, 0.7, 1) == (1.0, 1.0)
        @test factors(0.7, 0.7, 2) == (1.0, 1.0)
        @test factors(0.2, 0.7, 2)[2] < 1.0
        for voltage_scale in (0.1, 0.4, 0.8, 1.0)
            u = copy(initial)
            u[address["balance_d"]] .*= voltage_scale
            u[address["balance_q"]] .*= voltage_scale
            old = zeros(length(u))
            new = zeros(length(u))
            Base.invokelatest(old_balance!, old, u, ppost)
            balance!(new, u, ppost; low_voltage=false)
            @test new ≈ old atol=1e-12 rtol=1e-12
        end
        for characteristic in (1, 2)
            u = copy(initial)
            u[address["balance_d"]] .*= 0.4
            u[address["balance_q"]] .*= 0.4
            residual = y -> begin
                du = similar(y)
                rhs_with_load!(du, y, ppost, 0; characteristic=characteristic, zip=(0.7, 0.1, 0.2))
                du
            end
            direction = collect(range(-0.5, 0.5; length=length(u)))
            jv = ForwardDiff.jacobian(residual, u) * direction
            epsilon = 1e-6
            fd = (residual(u + epsilon*direction) - residual(u - epsilon*direction))/(2epsilon)
            @test jv ≈ fd rtol=1e-5 atol=1e-6
        end
    end

    rows = NamedTuple[]
    traces = NamedTuple[]
    max_iter = 200
    for (mix_name, mix) in (("PQ", (0.0, 0.0, 1.0)), ("ZIP_70_10_20", (0.7, 0.1, 0.2)))
        # All three laws coincide at this normal-voltage pre-event solution.
        pre_rhs! = (du, u, p, t) -> rhs_with_load!(du, u, p, t; zip=mix, low_voltage=false)
        pre_prob = MyDiffEq.ODEProblem(pre_rhs!, copy(initial), (0.0, 0.001), ppre, mass)
        pre_sol = MyDiffEq.Solve(pre_prob, 0.0005; method=:Euler, always_new=true, adaptive=false)
        @assert pre_sol.retcode == :Success
        previous = copy(pre_sol.u[end])
        pre_voltages = hypot.(previous[address["balance_d"]], previous[address["balance_q"]])
        @assert minimum(pre_voltages) > 0.7

        for (law, characteristic, enabled) in (("unmodified", 1, false), ("OpenIPSL_1", 1, true), ("OpenIPSL_2", 2, true))
            rhs! = (du, u, p, t) -> rhs_with_load!(du, u, p, t; zip=mix, characteristic=characteristic, low_voltage=enabled)
            for method in (:Euler, :Trap), h in (5e-4, 5e-5, 5e-6)
                println("TEST: $mix_name, $law, $method, h=$h")
                prob = MyDiffEq.ODEProblem(rhs!, copy(previous), (0.0, h), ppost, mass)
                elapsed = @elapsed sol = MyDiffEq.Solve(prob, h; method=method, adaptive=false, always_new=true, max_iter=max_iter, tol=1e-7)
                candidate = sol.retcode == :Success ? sol.u[end] : sol.newton_log.u_final
                residual = first_step_residual(rhs!, candidate, previous, ppost, mass, h, method)
                voltages = hypot.(candidate[address["balance_d"]], candidate[address["balance_q"]])
                vmin = minimum(voltages)
                residual_inf = norm(residual, Inf)
                verified = sol.retcode == :Success && residual_inf <= 1e-6
                # A tiny power mismatch at nearly zero V need not imply a tiny
                # current mismatch. Check the equivalent KCL residual separately.
                power_mismatch = hypot.(residual[address["balance_d"]], residual[address["balance_q"]])
                current_mismatch = maximum(v > 0 ? mismatch/v : Inf for (mismatch, v) in zip(power_mismatch, voltages))
                bus_vd = copy(models.bus.vd)
                bus_vq = copy(models.bus.vq)
                bus_vd[sys.non_slack_buses] = candidate[address["balance_d"]]
                bus_vq[sys.non_slack_buses] = candidate[address["balance_q"]]
                load_multipliers = map(models.load.bus) do b
                    v = hypot(bus_vd[b], bus_vq[b])
                    kP, kI = enabled ? factors(v, 0.7, characteristic) : (1.0, 1.0)
                    mix[1]*(v/models.bus.v[b])^2 + kI*mix[2]*(v/models.bus.v[b]) + kP*mix[3]
                end
                push!(rows, (load=mix_name, law=law, method=string(method), step_s=h,
                    retcode=string(sol.retcode), residual_verified=verified,
                    newton_iterations=sol.newton_log.iters, residual_inf=residual_inf,
                    converged_within_30_iterations=verified && sol.newton_log.iters <= 30,
                    max_current_mismatch_pu=current_mismatch,
                    current_balance_verified=verified && current_mismatch <= 1e-6,
                    candidate_min_voltage_pu=vmin,
                    candidate_min_load_multiplier=minimum(load_multipliers), elapsed_s=elapsed))
                for (k, r) in enumerate(sol.newton_log.residual_norm)
                    push!(traces, (load=mix_name, law=law, method=string(method), step_s=h,
                        iteration=k, residual_2=r, correction_2=sol.newton_log.correction_norm[k]))
                end
                # Save incrementally so completed comparisons survive interruption.
                CSV.write(joinpath(output_dir, "summary.csv"), DataFrame(rows))
                CSV.write(joinpath(output_dir, "newton_traces.csv"), DataFrame(traces))
                println("RESULT: $(sol.retcode), verified=$verified, iterations=$(sol.newton_log.iters), residual_inf=$residual_inf, current_mismatch=$current_mismatch")
            end
        end
    end
    println("\nResults saved to $output_dir (elapsed_s includes first-call compilation).")
    show(stdout, MIME("text/plain"), DataFrame(rows))
    println()
end

main()
