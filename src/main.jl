using Pkg; Pkg.activate(".")

include("Models.jl"); using .Models
include("data_loader.jl"); using .DataLoader
include("utils.jl"); using .Utils
include("system.jl"); using .SystemModel
include("power_flow.jl"); using .PowerFlow
include("static_init.jl"); using .StaticInit
include("dynamic_sim.jl"); using .DynamicSim
include("SteadyStateAnalysis.jl"); using .SteadyStateAnalysis

using MyDiffEq, Plots, Printf

# 1. Load data
# data_file = "cases/Fault_Cases/ieee14_fault_barq.xlsx"
# data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
# data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
models = load_data(data_file)
sys = build_system(models)

# 2. Solve power flow
solve_power_flow!(sys)

@show sys.models.bus.v
@show sys.models.bus.theta

# 3. Static initialization
run_static_init!(sys)

# 4. Dynamic simulation setup
address = build_dynamic_address(sys)
mass_matrix = build_mass_matrix(sys, address)
u0 = build_initial_conditions(sys, address)

# 4b. Small-signal analysis around the steady state (fault off, t=0)
ssa = small_signal_analysis(solve_dynamic_sim!, sys, address, mass_matrix, u0;
                            case_file=data_file, report_dir="reports", top_n=5)

lambda = 1.0
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
models.fault.t_fault=[100.]
models.fault.t_clear=[100.]
# tstops = [sys.models.fault.t_fault[1], sys.models.fault.t_clear[1], 2.0] 
# tstops = [1.0, 1.1]
# tstops = [2.0]
tstops = []

# 5. Pre-fault simulation
prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, 5.0), p, mass_matrix)
sol = MyDiffEq.Solve(prob, 5e-4, method=:Euler, adaptive=false, tstops=tstops)

states = MyDiffEq.get_states(sol)

# @show sol.iters[1:5]

using Plots
plot(sol.time[2:end], sol.iters)
plot(sol.time[2:end], sol.dt_hist)


plot(sol.time, states[address["delta"],:]')
plot(sol.time, states[address["omega"],:]')
plot(sol.time, states[address["line_id"],:]')
plot(sol.time, states[address["line_iq"],:]')
plot(sol.time, states[address["fault_id"],:]')
plot(sol.time, states[address["fault_iq"],:]')
plot(sol.time, states[address["balance_d"],:]')
plot(sol.time, states[address["balance_q"],:]')

v = complex.(states[address["balance_d"],:], states[address["balance_q"],:])
plot(sol.time, abs.(v'))

