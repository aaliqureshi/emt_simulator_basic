using Pkg; Pkg.activate(".")
using Revise
using Barq
using MyDiffEq, Plots

# 1. Load data
# data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
# data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
# data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
# data_file = "cases/Simple_Cases/wecc_full_slack.xlsx"
# data_file = "cases/Fault_Cases/case2383wp_gc.xlsx"
data_file = "cases/Fault_Cases/case118_gc.xlsx"
# data_file = "cases/Fault_Cases/case3012wp_barq.xlsx"
models = load_data(data_file)

# build system
sys = build_system(models)

# 1, 2, 5, 6, 9, 10, 11, 13, 14, 17, 19, 22
# 38,63,64,68,71,81
# models.fault.bus = [24]
# models.fault.bus = [3]
# models.fault.bus = [24]
# models.fault.x_fault[1] = 0.015
# models.fault.bus = [6]
# models.fault.x_fault[1] = 0.002

# models.fault.bus = [1396]
# models.fault.x_fault[1] = 0.01

# case 118 fault
models.fault.bus=[48]
models.fault.x_fault[1]=0.001
# models.fault.x_fault[1]=0.03

# models.load.p[:] .*= 1.42

# models.load.p[end-553] = 4.0

# 2. Solve power flow
solve_power_flow!(sys);

# 3. Static initialization
run_static_init!(sys)

# 4. Dynamic simulation setup
address = build_dynamic_address(sys);
mass_matrix = build_mass_matrix(sys, address)
u0 = build_initial_conditions(sys, address);


function run_simulation(u0, lambda, dt; 
                        method=:Euler, 
                        always_new=true, 
                        tstops=[], 
                        adaptive=false, 
                        t_end=0.1,
    )
    p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
    prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, t_end), p, mass_matrix)
    sol = MyDiffEq.Solve(prob, dt, method=method, adaptive=adaptive, tstops=tstops, always_new=always_new)
    return sol
end

# Pre-fault simulation
lambda=0.0
dt0=5e-4
# method=:Trap
method=:Euler
sol_pf=run_simulation(u0, lambda, dt0, t_end=2*dt0, method=method);
u1 = sol_pf.u[end]


vd_fault_idx = address["balance_d"][models.fault.bus[1]]
vq_fault_idx = address["balance_q"][models.fault.bus[1]]

# fault-on simulation 3
dt = 5e-4
lambda=1.0

# ux = run_simulation(u1, lambda, dt_list[iter], t_end=2*dt_list[iter],method=method)
ux = run_simulation(u1, lambda, 
                    dt, t_end=0.1,
                    method=method)
# ux = run_simulation(u1, lambda, 
#                     dt, t_end=dt,
#                     method=method)



# # re-init
#
λ_target = 1.0
# lambda=1.0
always_new=true
p_direct = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, λ_target)
p_base   = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)

u0_new = copy(u1)
u0_homotopy = copy(u1)
u0_adapt_h = copy(u1)
u0_resid_h = copy(u1)

solver_stats = true
r1 = solve_newton!(u0_new, p_direct, address; max_iter=100, always_new=always_new, solver_stats=solver_stats)
r2 = solve_homotopy!(u0_homotopy, 
                    p_base, address; λ_target=λ_target, 
                    Δλ=0.01, always_new=always_new,
                    vd_idx=vd_fault_idx,
                    vq_idx=vq_fault_idx,
                    )
r3 = solve_adaptive_homotopy!(u0_adapt_h, 
                              p_base, 
                              address; 
                              λ_target=λ_target, 
                              always_new=always_new,
                              vd_idx=vd_fault_idx,
                              vq_idx=vq_fault_idx,
                              )

# =#

# 5. post re-init simulation
# lambda = 0.0
# u3 = copy(ux.u[end])
# u3= copy(u0_homotopy)
# u3= copy(u0_new)
u3 = copy(u1)
# u3[address["line_iq"][end]+1:end] .= u0_new[address["line_iq"][end]+1: end]
dt_post = 5e-4
# method=:Trap
sol_fault = run_simulation(u3, lambda, dt_post, t_end=0.1, method=method)
# sol_fault = run_simulation(u3, lambda, dt_post, t_end=dt_post, method=method)




# # re-init
#
λ_target = 0.0
# lambda=1.0
always_new=true
p_direct = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, λ_target)
p_base   = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)

u0_new = copy(sol_fault.u[end])
u0_homotopy = copy(sol_fault.u[end])
u0_adapt_h = copy(sol_fault.u[end])
u0_resid_h = copy(sol_fault.u[end])

solver_stats = true
r1 = solve_newton!(u0_new, p_direct, address; max_iter=100, always_new=always_new, solver_stats=solver_stats)
r2 = solve_homotopy!(u0_homotopy, 
                    p_base, address; 
                    λ_start=1.0,
                    λ_target=λ_target, 
                    Δλ=0.01, always_new=always_new,
                    vd_idx=vd_fault_idx,
                    vq_idx=vq_fault_idx,
                    )
r3 = solve_adaptive_homotopy!(u0_adapt_h, 
                              p_base, 
                              address; 
                              λ_start=1.0,
                              λ_target=λ_target, 
                              always_new=always_new,
                              vd_idx=vd_fault_idx,
                              vq_idx=vq_fault_idx,
                              )

# =#


u4_flat=copy(sol_fault.u[end])
u4_flat[address["balance_d"]] .= 1.0
u4_flat[address["balance_q"]] .= 0.0 

# 5. post re-init simulation
lambda = 0.0
# u4 = copy(ux.u[end])
# u4= copy(u0_new)
u4= copy(sol_fault.u[end])
# u4 = copy(u0_homotopy)
# u4[address["line_iq"][end]+1:end] .= u0_new[address["line_iq"][end]+1: end]
dt_post = 1e-4
# method=:Trap
sol_post = run_simulation(u4, lambda, dt_post, t_end=0.2, method=method)
# sol_post = run_simulation(u4, lambda, dt_post, t_end=dt_post, method=method)
# sol_post = run_simulation(u4_flat, lambda, dt_post, t_end=0.5, method=method)





id_pre = [u[address["line_id"]] for u in sol_pf.u]
iq_pre = [u[address["line_iq"]] for u in sol_pf.u]
id_f = [u[address["line_id"]] for u in sol_fault.u]
iq_f = [u[address["line_iq"]] for u in sol_fault.u]
id_post = [u[address["line_id"]] for u in sol_post.u]
iq_post = [u[address["line_iq"]] for u in sol_post.u]

i_pre = stack(id_pre) + 1im*stack(iq_pre)
i_f =  stack(id_f) + 1im*stack(iq_f)
i_post = stack(id_post) + 1im*stack(iq_post)

i = hcat(i_pre, i_f, i_post)
i_mag = @. abs(i)
plot(i_mag[1:end,1:300]', label="")


vd_pre = [u[address["balance_d"]] for u in sol_pf.u]
vq_pre = [u[address["balance_q"]] for u in sol_pf.u]
vd_f = [u[address["balance_d"]] for u in sol_fault.u]
vq_f = [u[address["balance_q"]] for u in sol_fault.u]
vd_post = [u[address["balance_d"]] for u in sol_post.u]
vq_post = [u[address["balance_q"]] for u in sol_post.u]

v_pre = stack(vd_pre) + 1im*stack(vq_pre)
v_f =  stack(vd_f) + 1im*stack(vq_f)
v_post = stack(vd_post) + 1im*stack(vq_post)

v = hcat(v_pre, v_f, v_post)
v_mag = @. abs(v)

plot(v_mag[1:end,:]')

plot(v_mag[2,:])


norm(sol_fault_np - sol_post.u[end], Inf)
norm(sol_fault_np - sol_post.u[end], 2)

norm(sol_fault_np[address["delta"]] - sol_post.u[end][address["delta"]], Inf)
norm(sol_fault_np[address["omega"]] - sol_post.u[end][address["omega"]], Inf)
norm(sol_fault_np[address["balance_d"]] - sol_post.u[end][address["balance_d"]], Inf)
norm(sol_fault_np[address["balance_q"]] - sol_post.u[end][address["balance_q"]], Inf)




omega_pre = [u[address["omega"]] for u in sol_pf.u]
omega_f = [u[address["omega"]] for u in sol_fault.u]
omega_post = [u[address["omega"]] for u in sol_post.u]

omega_pre = stack(omega_pre)
omega_f =  stack(omega_f) 
omega_post = stack(omega_post)

omega = hcat(omega_pre, omega_f, omega_post)
# v_mag = @. abs(v)
plot(omega[1:end,:]')

# omega_pqbrak = copy(omega)
# v_pqbrak = copy(v_mag)

idx = 15
plot(omega_pqbrak[idx,:], label="pqbrak")
plot!(omega[idx,:], label="zip")

plot(v_pqbrak[idx,:], label="pqbrak")
plot!(v_mag[idx,:], label="zip")
