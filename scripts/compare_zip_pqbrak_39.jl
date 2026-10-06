using Pkg; Pkg.activate(".")
using Revise
using Barq
using MyDiffEq, Plots

# 1. Load data
data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
# data_file = "cases/Fault_Cases/case118_gc.xlsx"
models = load_data(data_file)

# build system
sys = build_system(models)

models.fault.bus = [20]
# models.fault.x_fault[1] = 0.013
models.fault.x_fault[1] = 0.009

# case 118 fault
# models.fault.bus=[12]
# models.fault.x_fault[1]=0.009
# models.fault.x_fault[1]=0.009

models.load.p[:] *= 1.42


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
# vd_fault_idx = address["fault_id"][end]
# vq_fault_idx = address["fault_iq"][end]


# re-init
λ_target = 1.0
# lambda=1.0
always_new=true
p_direct = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, λ_target)
p_base   = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)

u0_new = copy(u1)
u0_homotopy = copy(u1)
u0_adapt_h = copy(u1)

solver_stats = true
r1 = solve_newton!(u0_new, p_direct, address; max_iter=100, always_new=always_new, solver_stats=solver_stats)



# 5. post re-init simulation
lambda = 1.0
# u3 = copy(u0_new)
u3= copy(u1)
# u3= copy(u1)
dt_post = 5e-4
# method=:Trap
# sol_post = run_simulation(u3, lambda, dt_post, t_end=dt_post, method=method)
# sol_post = run_simulation(u3, lambda, dt_post, t_end=0.1, method=method)
sol_post = run_simulation(u3, lambda, dt_post, t_end=0.2, method=method)
# sol_post = run_simulation(u3, lambda, dt_post, t_end=0.0095, method=method)


# 5. remove fault simulation
lambda = 0.0
u4 = sol_post.u[end]
dt_post = 5e-3
sol_ss = run_simulation(u4, lambda, dt_post, t_end=10.0, method=method, always_new=true)



using Plots

omega_pre = [u[address["omega"]] for u in sol_pf.u]
omega_f = [u[address["omega"]] for u in sol_post.u]
omega_s = [u[address["omega"]] for u in sol_ss.u]

omega = Matrix(undef, 9, length(omega_pre)+length(omega_f)+length(omega_s))

col_idx=1
for data in vcat(omega_pre, omega_f, omega_s)
    omega[:,col_idx] = data
    col_idx+=1
end

# plot(omega[5,:])

# omega22=copy(omega)

# plot(omega[6,:], label="constant ZIP load")
# plot!(omega22[6,:], label="PQBRAL-type conversion")

vd_pre = [u[vd_fault_idx] for u in sol_pf.u]
vd_f = [u[vd_fault_idx] for u in sol_post.u]
vd_s = [u[vd_fault_idx] for u in sol_ss.u]

vq_pre = [u[vq_fault_idx] for u in sol_pf.u]
vq_f = [u[vq_fault_idx] for u in sol_post.u]
vq_s = [u[vq_fault_idx] for u in sol_ss.u]

v_pre = @. abs(vd_pre + 1im*vq_pre)
v_f = @. abs(vd_f + 1im*vq_f)
v_s = @. abs(vd_s + 1im*vq_s)

V = vcat(v_pre, v_f, v_s)

# plot(V)

# V_pqbrak = copy(V)
# omega_pqbrak = copy(omega)

# plot(V, label="ZIP")
# plot!(V_pqbrak, label="PQBRAK")

# x = 1:length(V)
# plot(x, omega[6,:], label="ZIP")
# plot!(x, omega_pqbrak[6,:], label="PQBRAK")
# plot!(twinx(), x,  V)

# plt = plot()

# for idx in collect(1:9)
#     plot!(plt, x, omega[idx, :])
# end

# for idx in collect(1:9)
#     plot!(plt, x, omega_pqbrak[idx, :], linestyle=:dash)
# end

# display(plt)

using JSON
# fname_v="v_zip.json"
# fname_w="omega_zip.json"
fname_v="v_pqbrak.json"
fname_w="omega_pqbrak.json"

open(fname_v, "w") do io
    JSON.print(io, V)
end

open(fname_w, "w") do io
    JSON.print(io, omega)
end

omega_zip = Matrix(undef, 9, length(omega_pre)+length(omega_f)+length(omega_s))
omega_pqbrak = Matrix(undef, 9, length(omega_pre)+length(omega_f)+length(omega_s))

# col_idx=1
# for data in w_zip
#     omega_zip[:,col_idx] = data
#     col_idx+=1
# end


v_zip = open("v_zip.json", "r") do io
    JSON.parse(io)
end

w_zip = open("omega_zip.json", "r") do io
    JSON.parse(io)
end

v_brak = open("v_pqbrak.json", "r") do io
    JSON.parse(io)
end

w_brak = open("omega_pqbrak.json", "r") do io
    JSON.parse(io)
end

plot(v_zip, label = "ZIP")
plot!(v_brak, label="PQBRAK")

col_idx=1
for (data_zip, data_brak) in zip(w_zip, w_brak)
    omega_zip[:,col_idx] = data_zip
    omega_pqbrak[:,col_idx] = data_brak
    col_idx+=1
end


idx = 5
plot(omega_zip[idx, :], label="ZIP")
plot!(omega_pqbrak[idx, :], label="PQBRAK")

# plot(omega_zip[idx, :] - omega_pqbrak[idx, :])