using Pkg; Pkg.activate(".")

using Barq

using MyDiffEq

"""Run one non-adaptive simulation phase, returning the solution or the thrown exception."""
function _solve_phase(u0, address, mass_matrix, sys, lambda, t_end; dt=5e-4, method=:Euler)
    p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
    prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, t_end), p, mass_matrix)
    try
        return MyDiffEq.Solve(prob, dt,
                              method=method,
                              adaptive=false,
                              tstops=[],
                              always_new=true,
                              )
    catch e
        return e
    end
end

"""Return the retcode of a phase result, or the exception type if the phase threw."""
_retcode(result) = result isa Exception ? Symbol(nameof(typeof(result))) : result.retcode

"""Sweep fault bus and fault reactance to locate fault-on cases that fail with `:MaxIter`."""
function main()
    # data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
    data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
    # data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
    x_fault_sweep = range(0.01, 0.05, step=0.001)

    models = load_data(data_file)
    slack_bus = models.slack.bus[1]
    non_slack_buses = setdiff(models.bus.idx, slack_bus)

    maxiter_hits = Tuple{Int,Float64}[]
    other_fails = Tuple{Int,Float64,Symbol}[]
    test_count = 0
    success_count = 0

    for bus in non_slack_buses
        println("================================================")
        println("Testing fault for bus $bus")
        println("================================================")
        models = load_data(data_file)
        models.fault.bus = [bus]
        sys = build_system(models)
        solve_power_flow!(sys)
        run_static_init!(sys)
        address = build_dynamic_address(sys)
        mass_matrix = build_mass_matrix(sys, address)
        u0 = build_initial_conditions(sys, address)

        # pre-fault (lambda = 0 => fault removed) does not depend on x_fault, so solve it once per bus
        sol0 = _solve_phase(u0, address, mass_matrix, sys, 0.0, 0.01)
        if _retcode(sol0) != :Success
            println("Pre-fault failed for bus $bus (retcode = $(_retcode(sol0))) -- skipping this bus")
            continue
        end
        u1 = sol0.u[end]

        for x_fault in x_fault_sweep
            models.fault.x_fault[1] = x_fault
            sol1 = _solve_phase(u1, address, mass_matrix, sys, 1.0, 0.1)
            test_count += 1
            retcode = _retcode(sol1)
            if retcode == :Success
                success_count += 1
            elseif retcode == :MaxIter
                push!(maxiter_hits, (bus, x_fault))
                println(">>> MaxIter failure: bus = $bus, x_fault = $x_fault")
            else
                push!(other_fails, (bus, x_fault, retcode))
                println(">>> $retcode failure: bus = $bus, x_fault = $x_fault")
            end
        end
    end

    println("================================================")
    println("Test count: $test_count")
    println("Success count: $success_count")
    println("MaxIter failures: $(length(maxiter_hits))")
    for (bus, x_fault) in maxiter_hits
        println("    bus = $bus, x_fault = $x_fault")
    end
    println("Other failures: $(length(other_fails))")
    for (bus, x_fault, retcode) in other_fails
        println("    bus = $bus, x_fault = $x_fault, retcode = $retcode")
    end
    println("================================================")

    return maxiter_hits
end

main()
