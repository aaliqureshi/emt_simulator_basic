using Pkg; Pkg.activate(".")
using Revise, LinearAlgebra

includet("Models.jl"); using .Models
includet("data_loader.jl"); using .DataLoader
includet("utils.jl"); using .Utils
includet("system.jl"); using .SystemModel
includet("power_flow.jl"); using .PowerFlow
includet("static_init.jl"); using .StaticInit
includet("dynamic_sim.jl"); using .DynamicSim
includet("SteadyStateAnalysis.jl"); using .SteadyStateAnalysis
includet("io/json.jl"); using .JsonRW
includet("../scripts/plot_pf_convergence.jl"); using .SolverPlot

using MyDiffEq, Plots, Printf

# using BenchmarkTools
# 1. Load data
# data_file = "cases/Fault_Cases/ieee14_fault_barq.xlsx"
# data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
# data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
# data_file = "cases/Fault_Cases/wecc_full.xlsx"
# data_file = "cases/Fault_Cases/kundur_full.xlsx"
# data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
# data_file = "cases/Fault_Cases/GBnetwork.xlsx"
# data_file = "cases/Fault_Cases/case118.xlsx"
data_file = "cases/Fault_Cases/case300.xlsx"
models = load_data(data_file)

# models.line.X[:] /= 2 
# models.line.R[:] .= 1e-5
# models.load.p[:] *=4.03

sys = build_system(models)

# models.bus.v[:] .= 1.0
# models.bus.theta[:] .= 0.0

# 2. Solve power flow
sol = solve_power_flow!(sys)

diff_vec = sol.u0 - sol.u_final
norm(diff_vec)

# label="GD (FE)"
# label="PGD (FE)"
# label="GD (BDF)"
label="PGD (BDF)"

# label="GD(BDF)"
# label="PGD(BDF)"
# label="PGD+GD (BDF)"
# label=""

# label="GD(FE)"
# label="PGD(FE)"
# label="PGD+GD (FE)"

# label="GD (FE) (fixed alpha)"
# label="GD (FE) (line search alpha)"

plt = plot_pf_convergence(sol, overlay=true, label=label)

# plt = plot_pf_convergence(sol, overlay=false, label=label)



# plt_eigen = plot_of_eigen_spectrum(sol)
#  plot(diff(sol.residual_norm[700:end]))

occursin("14", data_file) ? case="14" : nothing
occursin("SMIB", data_file) ? case="SMIB" : nothing
occursin("39", data_file) ? case="39" : nothing
occursin("wecc", data_file) ? case="wecc" : nothing


savefig(plt, "results/solver_convrg//convergence_residual_$case.pdf")


