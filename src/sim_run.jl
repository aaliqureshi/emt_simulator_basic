using Pkg; Pkg.activate(".")

using Barq

using MyDiffEq, Plots
# using BenchmarkTools, Memoize

function initialize_simulation(data_file)
    models = load_data(data_file)
    sys = build_system(models)
    solve_power_flow!(sys)
    run_static_init!(sys)
    address = build_dynamic_address(sys)
    mass_matrix = build_mass_matrix(sys, address)
    u0 = build_initial_conditions(sys, address)
    return sys, address, mass_matrix, u0
end

function run_simulation(data_file, method, dt)
    sys, address, mass_matrix, u0 = initialize_simulation(data_file)
    
    # fault-on simulation
    lambda = 0.0
    tstops = []
    p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
    tstops = []
    prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, 0.1), p, mass_matrix)
    sol = MyDiffEq.Solve(prob, dt, method=method, adaptive=false, tstops=tstops)
    u1 = sol.u[end]

    # fault-off simulation
    lambda = 1.0
    p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
    prob2 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u1, (0.0, 0.1), p, mass_matrix)
    sol2 = MyDiffEq.Solve(prob2, dt, method=method, adaptive=false, tstops=tstops)
    return sol2
end


begin
    # data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
    # data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
    data_file = "cases/Simple_Cases/wecc_full_slack.xlsx"
    # data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
    
    method = :Euler
    dt = 5e-6
end


sol2 = run_simulation(data_file, method, dt)

function convergence_rate(sol)
    return [sol.newton_log.correction_norm[i] / sol.newton_log.correction_norm[i+1] for i in range(1, length(sol.newton_log.correction_norm)-1)]
end


# cr = convergence_rate(sol2)

# k = 15
# n = collect(2:k)
# plot!(n, cr[2:k]; xlabel = "iteration", ylabel = "convergence rate", yaxis = :log, label = "$dt")

# plot!(n, sol2.newton_log.correction_norm[2:k], xlabel = "iteration", ylabel = "correction norm", yaxis = :log, label = "$dt")


# save results
using CSV
file_name = "wecc_$(dt)"
CSV.write("$(file_name)_correction_norm.csv", (correction_norm = sol2.newton_log.correction_norm,))
