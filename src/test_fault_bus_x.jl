using Pkg; Pkg.activate(".")

include("Models.jl"); using .Models
include("data_loader.jl"); using .DataLoader
include("utils.jl"); using .Utils
include("system.jl"); using .SystemModel
include("power_flow.jl"); using .PowerFlow
include("static_init.jl"); using .StaticInit
include("dynamic_sim.jl"); using .DynamicSim
include("homotopy.jl"); using .Homotopy

using MyDiffEq




function main()
    # data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
    data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
    # data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
    models = load_data(data_file)
    slack_bus = models.slack.bus[1]
    non_slack_buses = setdiff(models.bus.idx, slack_bus)
    method = :Trap
    test_count = 0
    success_count = 0

    for bus in non_slack_buses
        println("================================================")
        println("Testing fault for bus $bus")
        println("================================================")
        # fac = 0, 0.1, ..., 0.9 => line reactance at 100%, 90%, ..., 10% of original
        for fac in range(0.0, 0.1, length=2)
            println(fac)
            models = load_data(data_file)
            models.fault.bus = [bus]
            models.line.X[:] = (1-fac) * models.line.X[:]
            sys = build_system(models)
            solve_power_flow!(sys)
            run_static_init!(sys)
            address = build_dynamic_address(sys)
            mass_matrix = build_mass_matrix(sys, address)
            u0 = build_initial_conditions(sys, address)
            lambda = 0.0
            tstops = []
            p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
            prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, 0.1), p, mass_matrix)
            sol = MyDiffEq.Solve(prob, 5e-4, method=method, adaptive=false, tstops=tstops);
            test_count += 1
            if sol.retcode != :Success
                @show "!!!!!!!!!!!!!!Failed to solve fault-on for bus $bus!!!!!!!!!!!!!"
                continue
            else
                # println("-------Fault-on solved successfully-------")
                success_count += 1
            end
            
            u1 = sol.u[end]
            lambda = 1.0
            p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
            prob2 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u1, (0.0, 0.1), p, mass_matrix)
            sol2 = MyDiffEq.Solve(prob2, 5e-4, method=method, adaptive=false, tstops=tstops);
            test_count += 1
            if sol2.retcode != :Success
                @show "!!!!!!!!!!!!!!Failed to solve fault-off for bus $bus!!!!!!!!!!!!!"
                continue
            else
                # println("-------Fault-off solved successfully-------")
                success_count += 1
            end
        end


    end

    println("================================================")
    println("Test count: $test_count")
    println("Success count: $success_count")
    println("================================================")
end

main()