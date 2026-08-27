using Pkg; Pkg.activate(".")

using Barq

using MyDiffEq, Plots, Printf

# 1. Load data
# data_file = "cases/Fault_Cases/ieee14_fault_barq.xlsx"
# data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
# data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
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

lambda = 1.0
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
# tstops = [sys.models.fault.t_fault[1], sys.models.fault.t_clear[1], 2.0] 
# tstops = [1.0, 1.1]
# tstops = [2.0]
tstops = []

# 5. Pre-fault simulation
prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, 1.0), p, mass_matrix)
sol = MyDiffEq.Solve(prob, 5e-4, method=:Euler, adaptive=false, tstops=tstops)

states = MyDiffEq.get_states(sol)

@show sol.iters[1:5]

using Plots
plot(sol.time[2:end], sol.iters)
plot(sol.time[2:end], sol.dt_hist)


plot(sol.time, states[address["delta"],:]')
plot(sol.time, states[address["omega"],:]')
plot(sol.time, states[address["line_id"],:]')
plot(sol.time, states[address["line_iq"],:]')
plot(sol.time, states[address["fault_id"],:]')
plot(sol.time, states[address["fault_iq"],:]')
plot(sol.time, states[address["balance_d"],:]')
plot(sol.time, states[address["balance_q"],:]')

v = complex.(states[address["balance_d"],:], states[address["balance_q"],:])
plot(sol.time, abs.(v'))

nlog = sol.newton_log
nlog.time
nlog.iters
nlog.residual_norm
nlog.correction_norm
nlog.u_prev
nlog.u_final

# plot(2:nlog.iters, nlog.residual_norm[2:end])
plot(2:nlog.iters, nlog.correction_norm[2:end])



# perform homotopy on algebraic variables only
lambda_target = 1.0
lambda_step = 0.01
num_iters = []
u_homotopy = copy(u0)
for lambda in 0.0:lambda_step:lambda_target
    p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
    converged, iters = solve_algebraic!(u_homotopy, p, mass_matrix, address)
    println("λ=$lambda, converged=$converged, iters=$iters, R: $((1e10)^(1-lambda) * (models.fault.x_fault[1])^(lambda))")
    push!(num_iters, iters)
    if !converged
        break
    end
    # sleep(0.5)
end

plot(lambda_step:lambda_step:lambda_target, num_iters[2:end])
sum_iters = [ sum(num_iters[2:i]) for i in 2:length(num_iters) ]
plot!(lambda_step:lambda_step:lambda_target, sum_iters)

v_homotopy = complex.(u_homotopy[address["balance_d"]], u_homotopy[address["balance_q"]])
println("Bus voltage magnitudes at λ=0.803:")
for (i, bus) in enumerate(sys.non_slack_buses)
    println("  Bus $bus: |V| = $(abs(v_homotopy[i]))")
end

n_mech = length(address["delta"]) + length(address["omega"])
u0[n_mech+1:end] = u_homotopy[n_mech+1:end]


lambda = 0.8
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, 1.0), p, mass_matrix)
sol = MyDiffEq.Solve(prob, 5e-4, method=:Euler, adaptive=false, tstops=tstops)


# perform pseudo-arclength homotopy
# u_homotopy = copy(u0)                                                                                                
# p_base = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)
                                                                                                                     
# path_λ, path_u, success = pseudo_arclength!(u_homotopy, p_base, mass_matrix, address;                              
#                                              λ_start=0.0, λ_target=1.0, Δs=0.01,
#                                              Δλ_boot=0.01)

# # Plot the nose curve (voltage at a bus vs λ)
# using Plots
# n_mech = length(address["delta"]) + length(address["omega"])
# v14_mag = [abs(complex(pu[address["balance_d"][end] - n_mech],
#                         pu[address["balance_q"][end] - n_mech])) for pu in path_u]
# plot(path_λ, v14_mag, xlabel="λ", ylabel="|V₁₄|", marker=:circle)



# Step 1: Natural parameter continuation (fast, up to near the turning point)                                        
u_homotopy = copy(u0)
u_history = [copy(u0)]
λ_history = [0.0]
for lambda in 0.001:0.001:1.0
    p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
    converged, iters = solve_algebraic!(u_homotopy, p, mass_matrix, address)
    if !converged
        println("Natural parameter failed at λ=$lambda")
        break
    end
    push!(u_history, copy(u_homotopy))
    push!(λ_history, lambda)
end

# Step 2: Pseudo-arclength from the last two converged points
p_base = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)
u_prev = u_history[end-1]
u_last = u_history[end]
λ_prev = λ_history[end-1]
λ_last = λ_history[end]

path_λ, path_u, success = pseudo_arclength!(u_prev, λ_prev, u_last, λ_last,
                                             p_base, mass_matrix, address;
                                             λ_target=1.0)
                                            
  # Pick λ_target slightly below where homotopy fails                                                                  
#   λ_target = 0.80   
λ_target = 0.785
p_direct = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, λ_target)
p_base   = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)

u0_new = copy(u0)
u0_damped = copy(u0)
u0_backtracking = copy(u0)
u0_levenberg_marquardt = copy(u0)
u0_homotopy = copy(u0)
r1 = solve_newton!(u0_new,              p_direct, address; max_iter=1000)
r2 = solve_damped_newton!(u0_damped,       p_direct, address; α=0.5, max_iter=1000)
r3 = solve_backtracking_newton!(u0_backtracking, p_direct, address; max_iter=1000)
r4 = solve_levenberg_marquardt!(u0_levenberg_marquardt, p_direct, address; max_iter=1000)
r5 = solve_homotopy!(u0_homotopy,            p_base,   address; λ_target=λ_target, Δλ=0.001)

using Plots
plot(r1.residuals[1:100], label="Newton",       lw=2, yscale=:log10)
plot!(r2.residuals[1:100], label="Damped (α=0.5)", lw=2, yscale=:log10)
plot!(r3.residuals[1:100], label="Backtracking",   lw=2, yscale=:log10)
plot!(r4.residuals[1:end], label="Levenberg-Marquardt", lw=2, yscale=:log10)
plot!(r5.residuals[1:end], label="Homotopy",       lw=2, ls=:dash, yscale=:log10)
xlabel!("Newton iterations")
ylabel!("‖g‖")

u1_new = u0_new[]



