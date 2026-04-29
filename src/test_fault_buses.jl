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
    data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
    # data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
    # data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
    # data_file = "cases/Simple_Cases/wecc_full_slack.xlsx"
    models = load_data(data_file)
    slack_bus = models.slack.bus[1]
    non_slack_buses = setdiff(models.bus.idx, slack_bus)
    method = :Euler
    first_plot = true
    test_count = 0
    success_count = 0

    for bus in non_slack_buses
        local sol0, sol1
        local sol2
        println("================================================")
        println("Testing fault for bus $bus")
        println("================================================")
        models.fault.bus = [bus]
        models.fault.x_fault[1] = 0.015
        sys = build_system(models)
        solve_power_flow!(sys)
        run_static_init!(sys)
        address = build_dynamic_address(sys)
        mass_matrix = build_mass_matrix(sys, address)
        u0 = build_initial_conditions(sys, address)
        tstops = []
        
        lambda = 1.0
        p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
        prob0 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, 0.01), p, mass_matrix)
        try
            sol0 = MyDiffEq.Solve(prob0, 
                                  5e-4, 
                                  method=method, 
                                  adaptive=false, 
                                  tstops=tstops,
                                  always_new=true,
                                  )
        catch e
            println("Exception handling - steady state - bus = $bus !!!")
            continue
        end

        test_count += 1
        lambda = 0.0
        u1 = sol0.u[end]
        p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
        prob1 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u1, (0.0, 0.1), p, mass_matrix)
        try
            sol1 = MyDiffEq.Solve(prob1, 
                                  5e-4, 
                                  method=method,
                                  adaptive=false, 
                                  tstops=tstops,
                                  always_new=true,
                                  )
        catch e
            println("Exception handling - fault-on state - bus = $bus !!!")
            continue
        end
        if sol1.retcode != :Success
            @show "!!!!!!!!!!!!!!Failed to solve fault-on for bus $(bus)!!!!!!!!!!!!!"
            continue
        else
            success_count += 1
        end
        

        test_count += 1
        u2 = sol1.u[end]
        lambda = 1.0
        p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
        prob2 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u2, (0.0, 0.05), p, mass_matrix)
        try
            sol2 = MyDiffEq.Solve(prob2, 
                                  5e-4, 
                                  method=method, 
                                  adaptive=false, 
                                  tstops=tstops,
                                  always_new=true,
                                  )
        catch e
            println("Exception handling to solve fault-off for bus $(bus) !!!")
            continue
        end
        if sol2.retcode != :Success
            @show "!!!!!!!!!!!!!!Failed to solve fault-off for bus $(bus)!!!!!!!!!!!!!"
            continue
        else
            success_count += 1
        end
    end

    println("================================================")
    println("Test count: $test_count")
    println("Success count: $success_count")
    println("================================================")
end

main()