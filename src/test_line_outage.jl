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
    # data_file = "cases/Simple_Cases/wecc_full_slack.xlsx"
    # data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
    models = load_data(data_file)
    slack_bus = models.slack.bus[1]
    non_slack_buses = setdiff(models.bus.idx, slack_bus)
    method = :Euler
    first_plot = true
    test_count = 0
    success_count = 0

    lines = models.line.idx

    #scale load
    # models.load.p *= 2.0
    # models.load.q *= 2.0

    for line in lines
        local sol, sol2
        println("================================================")
        println("Testing fault for line $line")
        println("================================================")

        sys = build_system(models)
        solve_power_flow!(sys)
        run_static_init!(sys)
        address = build_dynamic_address(sys)
        mass_matrix = build_mass_matrix(sys, address)
        u0 = build_initial_conditions(sys, address)
        
        x = sys.models.line.X[line]
        r = sys.models.line.R[line]

        # sys.models.line.R[line] = 1e10
        # sys.models.line.X[line] = 1e10

        sys.models.line.R[line] *= 2.0
        sys.models.line.X[line] *= 2.0


        test_count += 1
        lambda = 1.0
        tstops = []
        p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
        prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, 0.1), p, mass_matrix)
        try
            sol = MyDiffEq.Solve(prob, 
                                5e-4, 
                                method=method, 
                                adaptive=false, 
                                tstops=tstops,
                                always_new=true,
                                )
        catch e
            print("****** Exception ******")
            continue
        end
        if sol.retcode != :Success
            println("!!!!!!!!!!!!!!Failed to solve fault-on for line $(line)!!!!!!!!!!!!!")
            continue
            # return sol
        else
            println("-------Fault-on solved successfully-------")
            success_count += 1
        end
        
        u1 = sol.u[end]
        # lambda = 1.0
        
        # revert line impedance
        sys.models.line.R[line] = r
        sys.models.line.X[line] = x

        test_count += 1
        p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
        prob2 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u1, (0.0, 0.1), p, mass_matrix)
        try
            sol2 = MyDiffEq.Solve(prob2, 
                                  5e-4,
                                  method=method,
                                  adaptive=false, 
                                  tstops=tstops,
                                  always_new=true,
                                  )    
        catch e
            println("****** Exception ******")
            continue
        end
        if sol2.retcode != :Success
            println("!!!!!!!!!!!!!!Failed to solve fault-off for line $(line)")
            continue
        else
            println("-------Fault-off solved successfully-------")
            success_count += 1
        end
    end

    println("================================================")
    println("Test count: $test_count")
    println("Success count: $success_count")
    println("================================================")
end

stats = main()