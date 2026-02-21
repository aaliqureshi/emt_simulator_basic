using Pkg; Pkg.activate(".")

include("Models.jl"); using .Models
include("data_loader.jl"); using .DataLoader
include("utils.jl"); using .Utils
include("system.jl"); using .SystemModel
include("power_flow.jl"); using .PowerFlow
include("static_init.jl"); using .StaticInit
include("dynamic_sim.jl"); using .DynamicSim
include("homotopy.jl"); using .Homotopy

using MyDiffEq, Plots

# 1. Load data
# data_file = "cases/Fault_Cases/ieee14_fault_barq.xlsx"
data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
# data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
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
u_homotopy = copy(u0)
for lambda in 0.0:0.01:1.0
    p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
    converged, iters = solve_algebraic!(u_homotopy, p, mass_matrix, address)
    println("λ=$lambda, converged=$converged, iters=$iters, R: $((1e10)^(1-lambda) * (models.fault.r_fault[1])^(lambda))")
    if !converged
        break
    end
    # sleep(0.5)
end

v_homotopy = complex.(u_homotopy[address["balance_d"]], u_homotopy[address["balance_q"]])
println("Bus voltage magnitudes at λ=0.803:")
for (i, bus) in enumerate(sys.non_slack_buses)
    println("  Bus $bus: |V| = $(abs(v_homotopy[i]))")
end




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