using Pkg; Pkg.activate(".")

include("Models.jl"); using .Models
include("data_loader.jl"); using .DataLoader
include("utils.jl"); using .Utils
include("system.jl"); using .SystemModel
include("power_flow.jl"); using .PowerFlow
include("static_init.jl"); using .StaticInit
include("dynamic_sim.jl"); using .DynamicSim

using MyDiffEq, Plots

# 1. Load data
data_file = "cases/Fault_Cases/ieee14_fault_barq.xlsx"
# data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
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
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)
tstops = [sys.models.fault.t_fault[1], sys.models.fault.t_clear[1]]

# 5. Pre-fault simulation
prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, 20.0), p, mass_matrix)
sol = MyDiffEq.Solve(prob, 50e-6, method=:Trap, adaptive=false, tstops=tstops)

states = MyDiffEq.get_states(sol)

using Plots
plot(sol.time[2:end], sol.iters)
plot(sol.time[2:end], sol.dt_hist)


plot(sol.time, states[address["delta"],:]')
plot(sol.time, states[address["omega"],:]')
