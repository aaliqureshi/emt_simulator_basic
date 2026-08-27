using Pkg; Pkg.activate(".")

using Barq

using MyDiffEq, Plots

# 1. Load data
# data_file = "cases/Simple_Cases/ieee39_fault.xlsx"
# data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
# data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
# data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx" 
data_file = "cases/Simple_Cases/wecc_full_slack.xlsx"
models = load_data(data_file)

# models.load.p[:] = models.load.p .* 4.0

# build system
sys = build_system(models)



# 2. Solve power flow
solve_power_flow!(sys)

# 3. Static initialization
run_static_init!(sys)

# 4. Dynamic simulation setup
address = build_dynamic_address(sys)
mass_matrix = build_mass_matrix(sys, address)
u0 = build_initial_conditions(sys, address)

models.fault.bus[1] = 20
models.fault.x_fault[1] = 0.015

# sys.models.line.X[line] = 1e10
# sys.models.line.R[line] = 1e10

lambda = 0.0
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
tstops = []
dt1 = 5e-4
prob1 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, 0.1), p, mass_matrix)
sol1 = MyDiffEq.Solve(prob1, dt1, method=:Euler, adaptive=false, tstops=tstops)
u1 = sol1.u[end]

u0 = build_initial_conditions(sys, address)
dt2= 5e-4
# dt2 = dt1/2
prob2 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, 0.1), p, mass_matrix)
sol2 = MyDiffEq.Solve(prob2, dt2, method=:Euler, adaptive=false, tstops=tstops)
# u1 = sol.u[end]



λ_target = 1.0


p_direct = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, λ_target)
p_base   = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)

# u0_new = copy(u0)
u0_new = copy(sol2.u[end])
u0_damped = copy(u0)
u0_backtracking = copy(u0)
u0_levenberg_marquardt = copy(u0)
u0_homotopy = copy(u0)
u0_homotopy2 = copy(u0)
u0_homotopy_lm = copy(u0)
r1 = solve_newton!(u0_new,              p_direct, address; max_iter=1000)
r2 = solve_damped_newton!(u0_damped,       p_direct, address; α=0.5, max_iter=1000)
r3 = solve_backtracking_newton!(u0_backtracking, p_direct, address; max_iter=1000)
r4 = solve_levenberg_marquardt!(u0_levenberg_marquardt, p_direct, address; max_iter=1000)
r5 = solve_homotopy!(u0_homotopy, p_base, address; λ_target=λ_target, Δλ=0.1)
# r6 = solve_homotopy2!(u0_homotopy, p_base, address; λ_target=λ_target, Δλ=0.1)
r7 = solve_homotopy_lm!(u0_homotopy_lm, p_base, address; λ_target=λ_target)




# u1_new = copy(r1.u)
# lambda = 1.0
# p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
# tstops = []
# prob2 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u1_new, (0.0, 0.1), p, mass_matrix)
# sol2 = MyDiffEq.Solve(prob2, 5e-4, method=:Euler, adaptive=false, tstops=tstops)
# states = get_states(sol2)
# plot(states[address["omega"],:])
# plot(sol2.newton_log.correction_norm)

norm_diverge1=sol1.newton_log.correction_norm
norm_diverge2=sol2.newton_log.correction_norm
norm_reinit=r1.correction_norm

# plot(1:length(norm_diverge), norm_diverge,
#      yscale = :log10,
#      label = "Standard Newton",
#      marker = :circle,
#      linewidth = 2)

# plot!(1:length(norm_reinit), norm_reinit,
#       label = "After reinitialization",
#       marker = :square,
#       linewidth = 2)

k = 10
g = 10
plot(sol1.newton_log.correction_norm[1:k], yscale = :log10, label="$(dt1)", marker = :circle, lw = 2.5, ms=5)
plot!(sol2.newton_log.correction_norm[1:g], yscale = :log10, label="$(dt2)", marker = :square, lw=2.5, ms = 5)
plot!(norm_reinit, yscale=:log10, marker = :diamond, lw = 3.0, ms=5, label="reinit", ylabel="‖Δxₖ‖", xlabel="Newton iteration", legend=:topleft)

savefig("14bus_3fault.pdf")

p1 = plot(collect(1:k), [norm_diverge1[1:k] norm_diverge2[1:k]], marker = [:circle :square])
p2 = plot(norm_reinit, marker = :square)

plot(p1, p2, layout=(2, 1), yscale=:log10, lw = 3.0, xlabel="Newton Iteration")

savefig("14bus_3fault_separate.pdf")



# u1_new = sol2.u[end]
u1_new = copy(u0_new)
lambda = 1.0
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
tstops = []
prob5 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u1_new, (0.0, 1.0), p, mass_matrix)
sol5 = MyDiffEq.Solve(prob5, 5e-4, method=:Euler, adaptive=false, tstops=tstops)
states = get_states(sol5)
plot(states[address["omega"],:])
# plot(sol2.newton_log.correction_norm)

# plot(states[address["delta"], :])

# plot(states[1,:])

v = complex.(states[address["balance_d"], :], states[address["balance_q"], :])

v_mag = abs.(v)
v_ang = angle.(v)

plot!(v_mag[1,:])

# t=sol5.time
# v_sin = @. v_mag[4,:] * cos(2*pi*60*t + v_ang[4,:])

# plot(sol5.time,v_sin)

# plot(sol5.time,v_sin