using Pkg; Pkg.activate(".")
using Revise

includet("Models.jl"); using .Models
includet("data_loader.jl"); using .DataLoader
includet("utils.jl"); using .Utils
includet("system.jl"); using .SystemModel
includet("power_flow.jl"); using .PowerFlow
includet("static_init.jl"); using .StaticInit
includet("dynamic_sim.jl"); using .DynamicSim
includet("SteadyStateAnalysis.jl"); using .SteadyStateAnalysis
includet("io/json.jl"); using .JsonRW

using MyDiffEq, Plots, Printf

# using BenchmarkTools
# 1. Load data
data_file = "cases/Fault_Cases/ieee14_fault_barq.xlsx"
# data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
# data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
# data_file = "cases/Fault_Cases/wecc_full.xlsx"
# data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
models = load_data(data_file)

# models.line.X[:] /= 2 
# models.line.R[:] .= 1e-5
# models.load.p[:]*=4.03

sys = build_system(models)

# models.bus.v[:] .= 1.0
# models.bus.theta[:] .= 0.0

# 2. Solve power flow
sol = solve_power_flow!(sys)
plot_pf_convergence(sol)
# sol.correction_norm
sol.residual_norm
diff(sol.residual_norm)

@show sys.models.bus.v
@show sys.models.bus.theta

# 3. Static initialization
run_static_init!(sys)

# 4. Dynamic simulation setup
address = build_dynamic_address(sys)
mass_matrix = build_mass_matrix(sys, address)
u0 = build_initial_conditions(sys, address)

# 4b. Small-signal analysis around the steady state (fault off, t=0)
# ssa = small_signal_analysis(solve_dynamic_sim!, sys, address, mass_matrix, u0;
#                             case_file=data_file, report_dir="reports", top_n=5)

lambda = 1.0
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
models.fault.t_fault=[1.0]
models.fault.t_clear=[1.1]
tstops = [sys.models.fault.t_fault[1], sys.models.fault.t_clear[1]] 
models.fault.x_fault .= 1e-4
# tstops = [1.0, 1.1]
# tstops = [2.0]
tstops = []

# 5. Pre-fault simulation
adaptive = false
always_new = true
dt = 50e-6
method = :Euler
NLsolver = :GD
prob0 = MyDiffEq.ODEProblem(solve_dynamic_sim!,
                           u0,
                           (0.0, 0.1),
                           p,
                           mass_matrix,
                           )


sol0 = MyDiffEq.Solve(prob0,
                     dt,
                     method=method,
                     adaptive=adaptive,
                     tstops=tstops,
                     always_new=always_new,
                     rtol=5e-3,
                     atol=1e-5,
                     NLsolver=NLsolver,
                     max_iter=10000,
                     alpha=1e-1,
                     )

states0 = MyDiffEq.get_states(sol0)


using Plots
lambda = 0.0; t=0
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
u1 = copy(sol0.u[end])

prob1 = MyDiffEq.ODEProblem(solve_dynamic_sim!,
                           u1,
                           (0.0, 0.1),
                           p,
                           mass_matrix,
                           )


sol1 = MyDiffEq.Solve(prob1,
                     dt,
                     method=method,
                     adaptive=adaptive,
                     tstops=tstops,
                     always_new=always_new,
                     rtol=5e-3,
                     atol=1e-5,
                     NLsolver=NLsolver,
                     max_iter=20000,
                     alpha=1e-13,
                     )

states1 = MyDiffEq.get_states(sol1)



mode = :gradient
# mode = :newton
gamma = -1e-1
# res1 = stabilize_algebraic!(u1, p, t, mass_matrix; mode = mode, tol = 1e-5, gamma = gamma)
res = stabilize_algebraic!(u1, p, t, mass_matrix; mode = mode, tol = 1e-8, gamma = 1e-8)

# @show res.converged, res.iters, res.g_norm

# plot!(res1.history, label="$gamma")
plot(res.history, label="$gamma")

# plot(res.history, yscale=:log10)

u2 = copy(res.u)
lambda = 0.0; t=0
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)

NLsolver=:GD

prob2 = MyDiffEq.ODEProblem(solve_dynamic_sim!,
                           u2,
                           (0.0, 0.1),
                           p,
                           mass_matrix,
                           )


sol2 = MyDiffEq.Solve(prob2,
                     dt,
                     method=method,
                     adaptive=adaptive,
                     tstops=tstops,
                     always_new=always_new,
                     rtol=5e-3,
                     atol=1e-5,
                     NLsolver=NLsolver,
                     max_iter=20000,
                     alpha=1e-12,
                     )

states2 = MyDiffEq.get_states(sol2)

plot(states2[11,:])

plot(1:res.iters, res.history[2:end], yscale=:log10, label = "", xlabel="iteration number", ylabel="log ||g||_2",
     title="Convergence behaviour of GD")