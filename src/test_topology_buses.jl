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

    # fault_lines = [1, 3, 4]
    lines = models.line.idx

    # models.load.p *= 1.0

    for line in lines
        println("================================================")
        println("Testing fault for line $line")
        println("================================================")
        # models.fault.bus = [bus]
        sys = build_system(models)
        solve_power_flow!(sys)
        run_static_init!(sys)
        address = build_dynamic_address(sys)
        mass_matrix = build_mass_matrix(sys, address)
        u0 = build_initial_conditions(sys, address)
        
        # sys.mo


        lambda = 0.0
        tstops = []
        p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
        prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, 0.1), p, mass_matrix)
        sol = MyDiffEq.Solve(prob, 5e-4, method=method, adaptive=false, tstops=tstops)
        test_count += 1
        if sol.retcode != :Success
            println("!!!!!!!!!!!!!!Failed to solve fault-on for line $(line)!!!!!!!!!!!!!")
            continue
            # return sol
        else
            println("-------Fault-on solved successfully-------")
            success_count += 1
        end
        
        u1 = sol.u[end]
        lambda = 1.0
        
        # change line impedance
        sys.models.line.R[fault_lines] .= 1e100
        sys.models.line.X[fault_lines] .= 1e100
        # sys.models.line.X[line] /= 2.0


        p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
        prob2 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u1, (0.0, 0.1), p, mass_matrix)
        sol2 = MyDiffEq.Solve(prob2, 5e-4, method=method, adaptive=false, tstops=tstops)
        test_count += 1
        if sol2.retcode != :Success
            println("!!!!!!!!!!!!!!Failed to solve fault-off for line $(line)")
            continue
        else
            println("-------Fault-off solved successfully-------")
            success_count += 1
        end

        # u2 = sol2.u[end]
        # v_post = complex.(u2[address["balance_d"]], u2[address["balance_q"]])
        # if first_plot
        #     plot(non_slack_buses, abs.(v_post), label="Fault Bus $bus", title="Post-fault voltage")
        #     first_plot = false
        # else
        #     plot!(non_slack_buses, abs.(v_post), label="Fault Bus $bus")
        # end
        # @show sol2.u[end][address["line_id"]][fault_lines]
        # @show sol2.u[end][address["balance_d"]]
    end

    println("================================================")
    println("Test count: $test_count")
    println("Success count: $success_count")
    println("================================================")
end

stats = main()