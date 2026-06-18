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
# data_file = "cases/Fault_Cases/ieee14_fault_barq.xlsx"
# data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
# data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
models = load_data(data_file)
sys = build_system(models)

# 2. Solve power flow
solve_power_flow!(sys)

# 3. Static initialization
run_static_init!(sys)

# 4. Dynamic simulation setup
address = build_dynamic_address(sys)
mass_matrix = build_mass_matrix(sys, address)
u0 = build_initial_conditions(sys, address)


sys.Y

z = inv(sys.Y)
z12 = z[1,2]
z13 = z[1,3]
z23 = z[2,3]