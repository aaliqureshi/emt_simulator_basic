using Pkg; Pkg.activate(".")

include("Models.jl"); using .Models
include("data_loader.jl"); using .DataLoader
include("utils.jl"); using .Utils
include("system.jl"); using .SystemModel
include("power_flow.jl"); using .PowerFlow
include("static_init.jl"); using .StaticInit
include("dynamic_sim.jl"); using .DynamicSim
include("homotopy.jl"); using .Homotopy

using MyDiffEq, Plots

# 1. Load data
# data_file = "cases/Simple_Cases/ieee39_fault.xlsx"
# data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
# data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
# data_file = "cases/Simple_Cases/wecc_full_slack.xlsx"
models = load_data(data_file)



# models.load.p *= 4.0

# build system
sys = build_system(models)



# 2. Solve power flow
solve_power_flow!(sys)


# 3. Static initialization
run_static_init!(sys)


# sys.models.line.X[4] *= 2.0
sys.models.fault.bus[1] = 3



# 4. Dynamic simulation setup
address = build_dynamic_address(sys)
mass_matrix = build_mass_matrix(sys, address)
u0 = build_initial_conditions(sys, address)



lambda = 0.72
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
tstops = []

# 5. Pre-fault simulation
prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, 0.1), p, mass_matrix)
sol = MyDiffEq.Solve(prob, 5e-4, method=:Euler, adaptive=false, tstops=tstops)
u1 = sol.u[end]




λ_target = 0.73
p_direct = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, λ_target)
p_base   = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)

u0_new = copy(u0)
# u0_damped = copy(u0)
u0_backtracking = copy(u0)
# u0_levenberg_marquardt = copy(u0)
u0_homotopy = copy(u0)
# u0_homotopy2 = copy(u0)
u0_homotopy_lm = copy(u0)
r1 = solve_newton!(u0_new,              p_direct, address; max_iter=1000)
# r2 = solve_damped_newton!(u0_damped,       p_direct, address; α=0.5, max_iter=1000)
r3 = solve_backtracking_newton!(u0_backtracking, p_direct, address; max_iter=1000)
# r4 = solve_levenberg_marquardt!(u0_levenberg_marquardt, p_direct, address; max_iter=1000)
r5 = solve_homotopy!(u0_homotopy, p_base, address; λ_target=λ_target, Δλ=0.01, max_iter = 25)
# r6 = solve_homotopy2!(u0_homotopy, p_base, address; λ_target=λ_target, Δλ=0.1)
r7 = solve_homotopy_lm!(u0_homotopy_lm, p_base, address; λ_target=λ_target, max_iter = 1000)

v0_homotopy = complex.(u0_homotopy[address["balance_d"]],u0_homotopy[address["balance_q"]])
# v0_trust = complex.(u0_levenberg_marquardt[address["balance_d"]],u0_levenberg_marquardt[address["balance_q"]])
v0_new = complex.(u0_new[address["balance_d"]],u0_new[address["balance_q"]])

  


