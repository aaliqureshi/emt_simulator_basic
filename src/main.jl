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
models = load_data(data_file)
sys = build_system(models)

# 2. Solve power flow
solve_power_flow!(sys)

@show sys.models.bus.v
@show sys.models.bus.theta

# 3. Static initialization
run_static_init!(sys)

@show sys.models.generator.delta
@show sys.models.generator.e_q_prime
@show sys.models.generator.i_d
@show sys.models.generator.i_q

# 4. Dynamic simulation setup
address = build_dynamic_address(sys)
mass_matrix = build_mass_matrix(sys, address)
u0 = build_initial_conditions(sys, address)
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)

# 5. Pre-fault simulation
prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, 10.0), p, mass_matrix)
sol_pre = MyDiffEq.Solve(prob, 50e-4, method=:Euler, adaptive=false)

states_pre = MyDiffEq.get_states(sol_pre)

# 6. Fault simulation
sys.models.fault.r_s = [1e-6]
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)

prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, sol_pre.u[end], (0.0, 0.1), p, mass_matrix)
sol_fault = MyDiffEq.Solve(prob, 50e-5, method=:Euler, adaptive=false)

states_fault = MyDiffEq.get_states(sol_fault)
