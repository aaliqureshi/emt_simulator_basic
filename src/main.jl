using Pkg; Pkg.activate(".")

include("Models.jl"); using .Models
include("data_loader.jl"); using .DataLoader
include("utils.jl"); using .Utils
include("system.jl"); using .SystemModel
include("power_flow.jl"); using .PowerFlow
include("static_init.jl"); using .StaticInit
include("dynamic_sim.jl"); using .DynamicSim
include("SteadyStateAnalysis.jl"); using .SteadyStateAnalysis
include("io/json.jl"); using .JsonRW

using MyDiffEq, Plots, Printf

using BenchmarkTools
# 1. Load data
data_file = "cases/Fault_Cases/ieee14_fault_barq.xlsx"
# data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
# data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
# data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
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
# ssa = small_signal_analysis(solve_dynamic_sim!, sys, address, mass_matrix, u0;
                            # case_file=data_file, report_dir="reports", top_n=5)

lambda = 1.0
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
models.fault.t_fault=[1.0]
models.fault.t_clear=[1.1]
tstops = [sys.models.fault.t_fault[1], sys.models.fault.t_clear[1]] 
models.fault.x_fault .= 1e-4
# tstops = [1.0, 1.1]
# tstops = [2.0]
# tstops = []

# 5. Pre-fault simulation
adaptive = true
always_new = true
dt = 5e-4
method = :Euler
prob = MyDiffEq.ODEProblem(solve_dynamic_sim!,
                           u0,
                           (0.0, 5.0),
                           p,
                           mass_matrix,
                           )


@elapsed sol = MyDiffEq.Solve(prob,
                     dt,
                     method=method,
                     adaptive=adaptive,
                     tstops=tstops,
                     always_new=always_new,
                     rtol=5e-3,
                     atol=1e-5,
                     )

states = MyDiffEq.get_states(sol)

@show sol.sol_stats

# # 6. Save results
# case_base = splitext(basename(data_file))[1]
# save_results(joinpath("results_full", case_base * ".json"), sol, sys, address;
#              case_file=data_file, solver=method, adaptive=adaptive, dt=dt)

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