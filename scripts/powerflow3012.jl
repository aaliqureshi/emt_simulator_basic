using Pkg; Pkg.activate(".")
using Revise
using Barq
using MyDiffEq, Plots

# 1. Load data
# data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
data_file = "cases/Fault_Cases/case3012wp_barq.xlsx"

models = load_data(data_file)

# build system
sys = build_system(models)

# 2. Solve power flow using NR
sol_nr = solve_power_flow!(sys);
# fails with flat start; works with case guess (4 iters)


# 2. Solve power flow using NR
sol_hc = solve_power_flow_continuation!(sys);


using LinearAlgebra

@show norm((sol_nr.u - sol_hc.u), Inf)

@show norm((sol_nr.u - sol_hc.u), 2)