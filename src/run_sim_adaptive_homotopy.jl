using Pkg; Pkg.activate(".")

include("Models.jl"); using .Models
include("data_loader.jl"); using .DataLoader
include("utils.jl"); using .Utils
include("system.jl"); using .SystemModel
include("power_flow.jl"); using .PowerFlow
include("static_init.jl"); using .StaticInit
include("dynamic_sim.jl"); using .DynamicSim
# include("homotopy.jl"); using .Homotopy

using MyDiffEq, Plots

# 1. Load data
# data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
# data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
# data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
data_file = "cases/Simple_Cases/wecc_full_slack.xlsx"
models = load_data(data_file)

# build system
sys = build_system(models)

models.fault.bus = [6]
models.fault.x_fault[1] = 0.0009

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
prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, 0.4), p, mass_matrix)
sol = MyDiffEq.Solve(prob, 5e-4, method=:Euler, adaptive=false, tstops=tstops, always_new=false)

u1 = sol.u[end]

# # 5. post-fault simulation
lambda = 1.0
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
prob2 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u1, (0.0, 0.1), p, mass_matrix)
sol2 = MyDiffEq.Solve(prob2, 5e-4, method=:Euler, adaptive=false, tstops=tstops, always_new=true)
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
always_new=true
# sys.models.line.X .= sys.models.line.X .* 1.01
p_direct = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, λ_target)
p_base   = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)

u0_new = copy(u1)
u0_homotopy = copy(u1)
u0_adapt_h = copy(u1)

r1 = solve_newton!(u0_new, p_direct, address; max_iter=100, always_new=always_new)
r2 = solve_homotopy!(u0_homotopy, p_base, address; λ_target=λ_target, Δλ=0.1, always_new=always_new)
r3 = solve_adaptive_homotopy!(u0_adapt_h, p_base, address; λ_target=λ_target, always_new=always_new)

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
always_new=true
prob3 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u3, (0.0, 5e-4), p, mass_matrix)
sol3 = MyDiffEq.Solve(prob3, 5e-4, method=:Euler, adaptive=false, tstops=tstops, always_new=always_new)
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





# λ_target = 0.764
λ_target = 0.0
# sys.models.line.X .= sys.models.line.X .* 1.01
p_direct = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, λ_target)
p_base   = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)

u0_new = copy(u3)
u0_homotopy = copy(u3)
u0_adapt_h = copy(u3)

r1 = solve_newton!(u0_new, p_direct, address; max_iter=1000, always_new=true)
r2 = solve_homotopy!(u0_homotopy, p_base, address; λ_target=λ_target, Δλ=0.1, always_new=true)
r3 = solve_adaptive_homotopy!(u0_adapt_h, p_base, address; λ_target=λ_target,always_new=true)

# v0_homotopy = complex.(u0_homotopy[address["balance_d"]],u0_homotopy[address["balance_q"]])
# v0_new = complex.(u0_new[address["balance_d"]],u0_new[address["balance_q"]])
# v0_homotopy2 = complex.(u0_adapt_h[address["balance_d"]],u0_adapt_h[address["balance_q"]])





lambda = 0.0
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
u4 = copy(u3)
# u4= copy(u2)
# u4= copy(u1)
# u4=copy(u0_new)
always_new=true
prob4 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u4, (0.0, 10e-4), p, mass_matrix)
sol4 = MyDiffEq.Solve(prob4, 5e-4, method=:Euler, adaptive=false, tstops=tstops,always_new=always_new)

states1 = get_states(sol)
states2 = get_states(sol3)
states3 = get_states(sol4)

states = hcat(states1, states2, states3)

plot(states[address["omega"],:][6,:])

v3_td = complex.(states[address["balance_d"], :], states[address["balance_q"], :])
v3_td_mag = abs.(v3_td)
v3_td_ang = angle.(v3_td)


t1 = sol.time
t2 = sol3.time
t3 = sol4.time
t2 = @. t2+5e-4 + t1[end] 
t3 = @. t3+5e-4 + t2[end]

t = vcat(t1,t2,t3)

# plot(diff(t))

plot(t, v3_td_mag[16,:])
plot(t, v3_td_mag[20,:])

v_sin = @. v3_td_mag[16,:]*cos(2*pi*60*t + v3_td_ang[16, :])

s = 500
e = 1500
plot(t[s:e], v_sin[s:e])

v_m = 1.0
sig = @. v_m * cos(2*pi*60*t)
plot(t[1:100], sig[1:100])



begin
    e = 10
    lw=1.5
    ms = 4
    tol_line = 1e-6
    # n1 = min(length(sol2.newton_log.residual_norm), 50)
    n1 = min(length(r1.residuals), 50)
    p1 = plot(collect(1:e), sol2.newton_log.residual_norm[1:e], 
         yscale=:log10, 
         lw=lw, 
         label="No initialization", 
         marker = :diamond,
         markersize=ms,
         markerstrokewidth=0,)
    # e2=length(r1.residuals)
    e2=length(sol3.newton_log.residual_norm)
    plot!(p1, collect(1:e2), 
        #   r1.residuals[1:e2],
          sol3.newton_log.residual_norm[1:e2],
          yscale=:log10, 
          xlabel = "Newton itertaion",
          ylabel = "Residual Norm ||f(z)||",
          lw=lw,
          label = "with initialization",
          marker=:square,
          markersize=ms,
          markerstrokewidth=0)
    # hline!([1e-6], label = "convergnce value", legend=(0.75, 0.25),
    hline!(p1, [tol_line], color=:black, linestyle=:dot, lw=1.0, label="")
    annotate!(n1 * 0.55, tol_line * 3, text(L"\epsilon = 10^{-6}", 8, :left))

    plot!(p1,
    xlabel="Newton iterations",
    # ylabel=L"\| g(\mathbf{z}) \|_2",
    ylabel="||g(x)||",
    # legend=:bottomright,
    legend=(0.65,0.20),
    legendfontsize=7,
    tickfontsize=8,
    guidefontsize=9,
    framestyle=:box,
    grid=true, gridalpha=0.2,
    size=(400, 300),
    dpi=300,
    margin=3Plots.mm)

display(p1)
# savefig(p1, "39bus20_newton_convergence_comparison.pdf")
    end

begin
        e = 10
        lw=1.5
        ms = 4
        tol_line = 1e-6
        n1 = min(length(r1.residuals), 50)
        # p1 = plot(collect(1:e), sol2.newton_log.residual_norm[1:e], 
        #      yscale=:log10, 
        #      lw=lw, 
        #      label="No initialization", 
        #      marker = :diamond,
        #      markersize=ms,
        #      markerstrokewidth=0,)
        e2=length(r1.residuals)
        # e2=length(sol3.newton_log.residual_norm)
        p1 = plot(collect(1:e2), 
              r1.residuals[1:e2],
            #   sol3.newton_log.residual_norm[1:e2],
              yscale=:log10, 
              xlabel = "Newton itertaion",
              ylabel = "Residual Norm ||f(z)||",
              lw=lw,
              label = "",
              marker=:square,
              markersize=ms,
              markerstrokewidth=0)
        # hline!([1e-6], label = "convergnce value", legend=(0.75, 0.25),
        hline!(p1, [tol_line], color=:black, linestyle=:dot, lw=1.0, label="")
        annotate!(n1 * 0.55, tol_line * 3, text(L"\epsilon = 10^{-6}", 8, :left))
    
        plot!(p1,
        xlabel="Newton iterations",
        # ylabel=L"\| g(\mathbf{z}) \|_2",
        ylabel="||g(x)||",
        # legend=:bottomright,
        legend=(0.65,0.20),
        legendfontsize=7,
        tickfontsize=8,
        guidefontsize=9,
        framestyle=:box,
        grid=true, gridalpha=0.2,
        size=(400, 300),
        dpi=300,
        margin=3Plots.mm)
    
    display(p1)
    # savefig(p1, "39bus20_newton_convergence_comparison.pdf")
        end

## ── IEEE figure: Newton convergence comparison ──
#
# Run both solvers from the SAME pre-event state u1, targeting the post-event λ.
# This isolates the reinitialization effect from any other differences.

# "Without reinitialization": plain Newton from pre-event state
u_no_reinit = copy(u1)
r_no_reinit = solve_newton!(u_no_reinit, p_direct, address; max_iter=100)

# "With reinitialization": homotopy continuation, then final Newton
u_reinit = copy(u1)
r_reinit = solve_adaptive_homotopy!(u_reinit, p_base, address; λ_target=λ_target)

# ── Panel 1: Newton convergence at t_event ──
using LaTeXStrings
begin
    tol_line = 1e-6
    lw = 1.5
    ms = 4

    n1 = min(length(r_no_reinit.residuals), 50)
    n2 = length(r_reinit.iter_hist)

    # Cumulative Newton iterations for homotopy (each λ-step costs some iterations)
    cum_iters = cumsum(r_reinit.iter_hist)

    p1 = plot(1:n1, r_no_reinit.residuals[1:n1],
         yscale=:log10, lw=lw, color=:red, linestyle=:dash,
         marker=:diamond, markersize=ms, markerstrokewidth=0,
         label="Without reinitialization")

    # For the homotopy, we don't have per-iteration residuals for each λ-step,
    # but we can show the final residual at each λ-step vs cumulative iterations.
    # A simpler and cleaner approach: just show the final Newton solve residuals
    # after homotopy has found consistent initial conditions.
    u_final_newton = copy(u_reinit)
    r_final = solve_newton!(u_final_newton, p_direct, address; max_iter=100)
    n3 = length(r_final.residuals)
    offset = r_reinit.total_newton_iters  # homotopy cost

    plot!(p1, offset .+ (1:n3), r_final.residuals[1:n3],
         lw=lw, color=:blue,
         marker=:circle, markersize=ms, markerstrokewidth=0,
         label="With reinitialization")

    # Mark the homotopy cost region
    vspan!(p1, [0, offset], fillalpha=0.08, color=:blue, label="Homotopy phase")

    hline!(p1, [tol_line], color=:black, linestyle=:dot, lw=1.0, label="")
    annotate!(p1, n1 * 0.55, tol_line * 3, text(L"\epsilon = 10^{-6}", 8, :left))

    plot!(p1,
         xlabel="Cumulative Newton iterations",
         ylabel=L"\| g(\mathbf{z}) \|_\infty",
         legend=:topright,
         legendfontsize=7,
         tickfontsize=8,
         guidefontsize=9,
         framestyle=:box,
         grid=true, gridalpha=0.2,
         size=(400, 300),
         dpi=300,
         margin=3Plots.mm)

    display(p1)
    # savefig(p1, "newton_convergence_comparison.pdf")
end