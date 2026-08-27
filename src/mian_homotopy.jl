using Pkg; Pkg.activate(".")

using Barq

using MyDiffEq, Plots

# 1. Load data
# data_file = "cases/Simple_Cases/ieee39_fault.xlsx"
data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
# data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
# data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
# data_file = "cases/Simple_Cases/wecc_full_slack.xlsx"
models = load_data(data_file)

x_og = copy(models.line.X)

# models.line.R[:] = models.line.R .* 0.5
# models.line.X[:] = models.line.X .* 0.048
models.line.X[:] = models.line.X .* 0.72

# models.load.p[:] = models.load.p .* 1.1
# models.load.q[:] = models.load.q .* 1.1
# models.generator.p_m[:] = models.generator.p_m .* 2.0

@show maximum(models.line.R ./ models.line.X)

# build system
sys = build_system(models)

models.fault.bus = [3]


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



lambda = 0.0
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
tstops = []

# 5. Pre-fault simulation
prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, 0.1), p, mass_matrix)
sol = MyDiffEq.Solve(prob, 5e-3, method=:Euler, adaptive=false, tstops=tstops)

u1 = sol.u[end]

# u1 = sol.u[end]
# # u1 = copy(u0)
# lambda = 0.0
# p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
# tstops = []

# sys.models.line.X[:] = x_og .* 0.01

# # 5. post-fault simulation
lambda = 1.0
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
prob2 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u1, (0.0, 0.1), p, mass_matrix)
sol2 = MyDiffEq.Solve(prob2, 5e-4, method=:Euler, adaptive=false, tstops=tstops)

# u2 = sol2.u[end]

v1 = complex.(u1[address["balance_d"]],u1[address["balance_q"]])
v0 = complex.(u0[address["balance_d"]],u0[address["balance_q"]])

@show v1 - v0
# complex.(u1[address["balance_d"]],u1[address["balance_q"]])
# complex.(u2[address["balance_d"]],u2[address["balance_q"]])

plot(models.bus.idx[2:end],abs.(v1), label="v1")
plot!(models.bus.idx[2:end],abs.(v0), label="v0")

plot!(abs.(v1 - v0), label="$(maximum(models.line.R ./ models.line.X))")

"""
@ λ_target = 0.0989, and using u0 as guess, all fail but homotopy converges.
"""


# λ_target = 0.764
λ_target = 1.0
# sys.models.line.X .= sys.models.line.X .* 1.01
p_direct = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, λ_target)
p_base   = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)

u0_new = copy(u0)
u0_damped = copy(u0)
u0_backtracking = copy(u0)
u0_levenberg_marquardt = copy(u0)
u0_homotopy = copy(u0)
u0_homotopy2 = copy(u0)
u0_homotopy_lm = copy(u0)
u0_adapt_h = copy(u0)
r1 = solve_newton!(u0_new,              p_direct, address; max_iter=1000)
r2 = solve_damped_newton!(u0_damped,       p_direct, address; α=0.5, max_iter=1000)
r3 = solve_backtracking_newton!(u0_backtracking, p_direct, address; max_iter=1000)
r4 = solve_levenberg_marquardt!(u0_levenberg_marquardt, p_direct, address; max_iter=1000)
r5 = solve_homotopy!(u0_homotopy, p_base, address; λ_target=λ_target, Δλ=0.1)
# r6 = solve_homotopy2!(u0_homotopy, p_base, address; λ_target=λ_target, Δλ=0.1)
r7 = solve_homotopy_lm!(u0_homotopy_lm, p_base, address; λ_target=λ_target)
r8 = solve_adaptive_homotopy!(u0_adapt_h, p_base, address)

v0_homotopy = complex.(u0_homotopy[address["balance_d"]],u0_homotopy[address["balance_q"]])
v0_trust = complex.(u0_levenberg_marquardt[address["balance_d"]],u0_levenberg_marquardt[address["balance_q"]])
v0_new = complex.(u0_new[address["balance_d"]],u0_new[address["balance_q"]])
# v0_homotopy2 = complex.(u0_homotopy2[address["balance_d"]],u0_homotopy2[address["balance_q"]])
v0_homotopy_lm = complex.(u0_homotopy_lm[address["balance_d"]],u0_homotopy_lm[address["balance_q"]])


plot(models.bus.idx[2:end],abs.(v0_homotopy))
plot!(models.bus.idx[2:end],abs.(v0_trust))
plot!(models.bus.idx[2:end],abs.(v0_new))
# plot!(models.bus.idx[2:end],abs.(v0_homotopy2))
plot!(models.bus.idx[2:end],abs.(v0_homotopy_lm))

plot(models.bus.idx[2:end],abs.(v0_trust))
plot!(models.bus.idx[2:end],abs.(v0_homotopy_lm))

plot(models.bus.idx[2:end],angle.(v0_homotopy))
plot!(models.bus.idx[2:end],angle.(v0_trust))
plot!(models.bus.idx[2:end],angle.(v0_new))
plot!(models.bus.idx[2:end],angle.(v0_homotopy2))




# using Plots
# plot(r1.residuals[1:100], label="Newton",       lw=2, yscale=:log10)
# plot!(r2.residuals[1:100], label="Damped (α=0.5)", lw=2, yscale=:log10)
# plot!(r3.residuals[1:100], label="Backtracking",   lw=2, yscale=:log10)
# plot!(r4.residuals[1:end], label="Levenberg-Marquardt", lw=2, yscale=:log10)
# plot!(r5.residuals[1:end], label="Homotopy",       lw=2, ls=:dash, yscale=:log10)
# xlabel!("Newton iterations")
# ylabel!("‖g‖")



# # solve DAE with new results

# lambda = 0.0439
# p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
# tstops = []

# # 5. Pre-fault simulation
# prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0_homotopy, (0.0, 5e-4), p, mass_matrix)
# sol2 = MyDiffEq.Solve(prob, 5e-4, method=:Euler, adaptive=false, tstops=tstops)


# u2 = sol2.u[end]

# λ_target = 1.0
# p_direct = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, λ_target)
# p_base   = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)

# u2_new = copy(u2)
# u2_damped = copy(u2)
# u2_backtracking = copy(u2)
# u2_levenberg_marquardt = copy(u2)
# u2_homotopy = copy(u2)
# r1 = solve_newton!(u2_new,              p_direct, address; max_iter=1000)
# r2 = solve_damped_newton!(u2_damped,       p_direct, address; α=0.5, max_iter=1000)
# r3 = solve_backtracking_newton!(u2_backtracking, p_direct, address; max_iter=1000)
# r4 = solve_levenberg_marquardt!(u2_levenberg_marquardt, p_direct, address; max_iter=1000)
# r5 = solve_homotopy!(u2_homotopy,            p_base,   address; λ_target=λ_target, Δλ=0.0001)


# v0_homotopy = complex.(u0_homotopy[address["balance_d"]],u0_homotopy[address["balance_q"]])
# v0_trust = complex.(u0_levenberg_marquardt[address["balance_d"]],u0_levenberg_marquardt[address["balance_q"]])
# v0_new = complex.(u0_new[address["balance_d"]],u0_new[address["balance_q"]])












plot(collect(1:20), r5.correction_norm[3:23], yscale=:log10)




# homotopy with u1 as guess

λ_target = 1.0
p_direct = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, λ_target)
p_base   = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)

u1_new = copy(u1)
u1_damped = copy(u1)
u1_backtracking = copy(u1)
u1_levenberg_marquardt = copy(u1)
u1_homotopy = copy(u1)
# u1_homotopy2 = copy(u1)
u1_homotopy_lm = copy(u1)
r1 = solve_newton!(u1_new,              p_direct, address; max_iter=1000)
r2 = solve_damped_newton!(u1_damped,       p_direct, address; α=0.5, max_iter=1000)
r3 = solve_backtracking_newton!(u1_backtracking, p_direct, address; max_iter=1000)
r4 = solve_levenberg_marquardt!(u1_levenberg_marquardt, p_direct, address; max_iter=1000)
r5 = solve_homotopy!(u1_homotopy,            p_base,   address; λ_target=λ_target, Δλ=0.1, max_iter=100)
# r6 = solve_homotopy2!(u1_homotopy2,            p_base,   address; λ_target=λ_target, Δλ=0.01)
r7 = solve_homotopy_lm!(u1_homotopy_lm, p_base, address; λ_target=λ_target)
# using Plots
# plot(r1.residuals[1:100], label="Newton",       lw=2, yscale=:log10)
# plot!(r2.residuals[1:100], label="Damped (α=0.5)", lw=2, yscale=:log10)
# plot!(r3.residuals[1:100], label="Backtracking",   lw=2, yscale=:log10)
# plot!(r4.residuals[1:end], label="Levenberg-Marquardt", lw=2, yscale=:log10)
# plot!(r5.residuals[1:end], label="Homotopy",       lw=2, ls=:dash, yscale=:log10)
# xlabel!("Newton iterations")
# ylabel!("‖g‖")

v_homotopy = complex.(u1_homotopy[address["balance_d"]],u1_homotopy[address["balance_q"]])
v_trust = complex.(u1_levenberg_marquardt[address["balance_d"]],u1_levenberg_marquardt[address["balance_q"]])
v_new = complex.(u1_new[address["balance_d"]],u1_new[address["balance_q"]])
v_backtrack = complex.(u1_backtracking[address["balance_d"]],u1_backtracking[address["balance_q"]])
v_damped = complex.(u1_damped[address["balance_d"]],u1_damped[address["balance_q"]])
# v_homotopy2 = complex.(u1_homotopy2[address["balance_d"]],u1_homotopy2[address["balance_q"]])
v_homotopy_lm = complex.(u1_homotopy_lm[address["balance_d"]],u1_homotopy_lm[address["balance_q"]])

@show v_homotopy - v_new

@show v_new - v1

plot(models.bus.idx[2:end],abs.(v_homotopy))
plot!(models.bus.idx[2:end],abs.(v_trust))
plot!(models.bus.idx[2:end],abs.(v_new))
# plot!(models.bus.idx[2:end],abs.(v_homotopy2))
plot!(models.bus.idx[2:end],abs.(v_homotopy_lm))

plot(models.bus.idx[2:end],angle.(v_homotopy))
plot!(models.bus.idx[2:end],angle.(v_trust))
plot!(models.bus.idx[2:end],angle.(v_new))

using Plots


