using Pkg; Pkg.activate(".")

using Barq

using MyDiffEq, Plots

# 1. Load data
# data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
# data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
# data_file = "cases/Simple_Cases/wecc_full_slack.xlsx"
models = load_data(data_file)

# build system
sys = build_system(models)

models.fault.bus = [16]

# 2. Solve power flow
solve_power_flow!(sys)

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
sol = MyDiffEq.Solve(prob, 5e-4, method=:Euler, adaptive=false, tstops=tstops)

u1 = sol.u[end]

# # 5. post-fault simulation
lambda = 1.0
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
prob2 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u1, (0.0, 0.1), p, mass_matrix)
sol2 = MyDiffEq.Solve(prob2, 5e-4, method=:Euler, adaptive=false, tstops=tstops)
u2 = sol2.u[end]

v0 = complex.(u0[address["balance_d"]], u0[address["balance_q"]])
v1 = complex.(u1[address["balance_d"]], u1[address["balance_q"]])
v2 = complex.(u2[address["balance_d"]], u2[address["balance_q"]])

# @show v1 - v0

plot(models.bus.idx[2:end],abs.(v1), label="v1")
plot!(models.bus.idx[2:end],abs.(v0), label="v0")
plot!(models.bus.idx[2:end], abs.(v2), label = "v2")


# λ_target = 0.764
λ_target = 1.0
# sys.models.line.X .= sys.models.line.X .* 1.01
p_direct = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, λ_target)
p_base   = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)

u0_new = copy(u1)
u0_homotopy = copy(u1)
u0_adapt_h = copy(u1)

r1 = solve_newton!(u0_new, p_direct, address; max_iter=1000)
r2 = solve_homotopy!(u0_homotopy, p_base, address; λ_target=λ_target, Δλ=0.1)
r3 = solve_adaptive_homotopy!(u0_adapt_h, p_base, address; λ_target=λ_target)

v0_homotopy = complex.(u0_homotopy[address["balance_d"]],u0_homotopy[address["balance_q"]])
v0_new = complex.(u0_new[address["balance_d"]],u0_new[address["balance_q"]])
v0_homotopy2 = complex.(u0_adapt_h[address["balance_d"]],u0_adapt_h[address["balance_q"]])


plot(models.bus.idx[2:end],abs.(v0_homotopy))
plot!(models.bus.idx[2:end],abs.(v0_new))
plot!(models.bus.idx[2:end],abs.(v0_homotopy2))


plot(models.bus.idx[2:end],angle.(v0_homotopy))
plot!(models.bus.idx[2:end],angle.(v0_new))
plot!(models.bus.idx[2:end],angle.(v0_homotopy2))




# # 5. post-fault simulation
lambda = 1.0
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
u3 = copy(u0_new)
# u3= copy(u2)
# u3= copy(u1)
prob3 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u3, (0.0, 0.1), p, mass_matrix)
sol3 = MyDiffEq.Solve(prob3, 5e-4, method=:Euler, adaptive=false, tstops=tstops)
u3 = sol3.u[end]

using OrdinaryDiffEq
prob00 = OrdinaryDiffEq.ODEFunction(solve_dynamic_sim!, mass_matrix=mass_matrix)
prob01 = OrdinaryDiffEq.ODEProblem(prob00, u3, (0.0, 0.1), p)
sol00 = OrdinaryDiffEq.solve(prob01, 
                            # ImplicitEuler(nlsolve=NLNewton(;always_new=true)),
                            ImplicitEuler(nlsolve=NLNewton(;always_new=false)),
                            adaptive=true,
                            dt=5e-4)

states3 = get_states(sol3)

v0 = complex.(u0[address["balance_d"]], u0[address["balance_q"]])
v1 = complex.(u1[address["balance_d"]], u1[address["balance_q"]])
v2 = complex.(u2[address["balance_d"]], u2[address["balance_q"]])
v3 = complex.(u3[address["balance_d"]], u3[address["balance_q"]])

v3_td_c = complex.(states3[address["balance_d"], :], states3[address["balance_q"]])
v3_td_mag = abs.(v3_td_c)

# @show v1 - v0

plot(models.bus.idx[2:end],abs.(v1), label="v1")
plot!(models.bus.idx[2:end],abs.(v0), label="v0")
plot!(models.bus.idx[2:end], abs.(v2), label = "v2")
plot!(models.bus.idx[2:end], abs.(v3), label = "v3")
plot(sol3.time, v3_td_mag[2,:])


lambda = 0.0
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
u4 = copy(u3)
# u3= copy(u2)
# u3= copy(u1)
prob4 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u4, (0.0, 0.1), p, mass_matrix)
sol4 = MyDiffEq.Solve(prob4, 5e-4, method=:Euler, adaptive=false, tstops=tstops)

states1 = get_states(sol)
states2 = get_states(sol3)
states3 = get_states(sol4)

states = hcat(states1, states2, states3)

plot(states[address["omega"],:][6,:])

v3_td = complex.(states[address["balance_d"], :], states[address["balance_q"]])
v3_td_mag = abs.(v3_td)

t1 = sol.time
t2 = sol3.time
t3 = sol4.time
t2 = @. t2+5e-4 + t1[end] 
t3 = @. t3+5e-4 + t2[end]

t = vcat(t1,t2,t3)

plot(diff(t))

plot(t, v3_td_mag[16,:])