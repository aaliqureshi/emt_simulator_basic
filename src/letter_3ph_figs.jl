using Pkg; Pkg.activate(".")
using Revise

includet("Models.jl"); using .Models
includet("data_loader.jl"); using .DataLoader
includet("utils.jl"); using .Utils
includet("system.jl"); using .SystemModel
includet("power_flow.jl"); using .PowerFlow
includet("static_init.jl"); using .StaticInit
includet("dynamic_sim.jl"); using .DynamicSim
# includet("homotopy.jl"); using .Homotopy

using MyDiffEq, Plots

# 1. Load data
# data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
# data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
# data_file = "cases/Simple_Cases/wecc_full_slack.xlsx"
models = load_data(data_file)

# build system
sys = build_system(models)

models.fault.bus = [20]
models.fault.x_fault[1] = 0.015

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
always_new=true
prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, 0.1), p, mass_matrix)
sol = MyDiffEq.Solve(prob, 5e-4, method=:Euler, adaptive=false, tstops=tstops, always_new=always_new)
u1 = sol.u[end]

# # 5. post-fault simulation
# lambda = 1.0
lambda = 1.0
always_new=true
dt1 = 5e-4
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
prob2 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u1, (0.0, 0.1), p, mass_matrix)
sol2 = MyDiffEq.Solve(prob2, dt1, method=:Euler, adaptive=false, tstops=tstops, always_new=always_new)
u2 = sol2.u[end]

# # 5. post-fault simulation with smaller dt
lambda = 1.0
always_true=true
dt2=5e-5
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
prob99 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u1, (0.0, 0.1), p, mass_matrix)
sol99 = MyDiffEq.Solve(prob99, dt2, method=:Euler, adaptive=false, tstops=tstops, always_new=always_new)
u99 = sol99.u[end]



# re-init
λ_target = 1.0
# lambda=1.0
always_new=true
p_direct = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, λ_target)
p_base   = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)

u0_new = copy(u1)
u0_homotopy = copy(u1)
u0_adapt_h = copy(u1)

r1 = solve_newton!(u0_new, p_direct, address; max_iter=100, always_new=always_new)
r2 = solve_homotopy!(u0_homotopy, p_base, address; λ_target=λ_target, Δλ=0.1, always_new=always_new)
r3 = solve_adaptive_homotopy!(u0_adapt_h, p_base, address; λ_target=λ_target, always_new=always_new)




# # 5. post re-init simulation
lambda = 0.0
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
u3 = copy(u0_new)
# u3= copy(u2)
# u3= copy(u1)
always_new=true
prob3 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u3, (0.0, 0.1), p, mass_matrix)
sol3 = MyDiffEq.Solve(prob3, 5e-4, method=:Euler, adaptive=false, tstops=tstops, always_new=always_new)
u3 = sol3.u[end]



# Optional - re-init post-fault
λ_target = 0.0
p_direct = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, λ_target)
p_base   = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses)

u1_new = copy(u3)
u1_homotopy = copy(u3)
u1_adapt_h = copy(u3)

r4 = solve_newton!(u1_new, p_direct, address; max_iter=1000, always_new=true)
r5 = solve_homotopy!(u1_homotopy, p_base, address; λ_target=λ_target, Δλ=0.1, always_new=true)
r6 = solve_adaptive_homotopy!(u1_adapt_h, p_base, address; λ_target=λ_target,always_new=true)


# post-fault integration

lambda = 0.0
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
u4 = copy(u3)
# u4= copy(u2)
# u4= copy(u1)
# u4=copy(u1_new)
always_new=true
prob4 = MyDiffEq.ODEProblem(solve_dynamic_sim!, u4, (0.0, 10e-4), p, mass_matrix)
sol4 = MyDiffEq.Solve(prob4, 5e-4, method=:Euler, adaptive=false, tstops=tstops,always_new=always_new)



using LaTeXStrings

# begin
#     lw = 1.5
#     ms = 4
#     tol_line = 1e-6

#     # ── Left panel (a): with vs without initialization ──
#     e = 10
#     pa = plot(1:e, sol2.newton_log.residual_norm[1:e],
#          yscale=:log10, lw=lw, color=:black, linestyle=:dash,
#          marker=:diamond, markersize=ms, markerstrokewidth=0,
#          label="No reinitialization (dt = 500 μs)")
        
#     plot!(1:e, sol99.newton_log.residual_norm[1:e],
#          yscale=:log10, lw=lw, color=:red, linestyle=:dash,
#          marker=:diamond, markersize=ms, markerstrokewidth=0,
#          label="No reinitialization (dt = 50 μs)")

#     e2 = length(sol3.newton_log.residual_norm)
#     plot!(pa, 1:e2, sol3.newton_log.residual_norm[1:e2],
#           yscale=:log10, lw=lw, color=:blue,
#           marker=:square, markersize=ms, markerstrokewidth=0,
#           label="With reinitialization")

#     hline!(pa, [tol_line], color=:black, linestyle=:dot, lw=1.0, label="")
#     plot!(pa,
#          xlabel="Newton iterations", ylabel=L"\| F(\mathbf{z}) \|_2",
#          title="", titlelocation=:center,
#          legend=(0.42, 0.40), legendfontsize=6,
#          tickfontsize=8, guidefontsize=9,
#          framestyle=:box, grid=true, gridalpha=0.2)

#     # ── Right panel (b): homotopy reinitialization detail ──
#     e3 = length(r1.residuals)
#     pb = plot(1:e3, r1.residuals[1:e3],
#          yscale=:log10, lw=lw, color=:orangered,
#          marker=:circle, markersize=ms, markerstrokewidth=0,
#          label="")

#     hline!(pb, [tol_line], color=:black, linestyle=:dot, lw=1.0, label="")
#     plot!(pb,
#          xlabel="Newton iterations", ylabel=L"\| g(\mathbf{z}) \|_2",
#          title="", titlelocation=:center,
#          legendfontsize=6, tickfontsize=8, guidefontsize=9,
#          framestyle=:box, grid=true, gridalpha=0.2)

#     # ── Combine side-by-side ──
#     fig = plot(pa, pb, layout=(1, 2), size=(700, 300), dpi=300,
#               margin=4Plots.mm, bottom_margin=5Plots.mm)
#     display(fig)
#     # savefig(fig, "39bus16_newton_convergence_comparison_complete.pdf")
# end

# begin
#      # lw = 3.5
#      # ms = 7
#      tol_line = 1e-6

# default(
#     fontfamily = "Computer Modern",
#     linewidth = 3.5,
#     markersize = 9,
#     legendfontsize = 13,
#     guidefontsize = 14,
#     tickfontsize = 14,
#     titlefontsize = 14,
#     dpi = 300,
#     size = (1200, 500),
#     # margin = 6Plots.m
#     )
 
#      # ── Left panel (a): with vs without initialization ──
#      e = 10
#      pa = plot(1:e, 
#                sol2.newton_log.residual_norm[1:e],
#                yscale=:log10, 
#                # lw=lw, 
#                color=:blue, 
#                linestyle=:dash,
#                marker=:diamond, 
#                # markersize=ms,
#                markerstrokewidth=0,
#                label="No re-initialization (h = 500μs)",
#                )
         
#      plot!(1:e, 
#            sol99.newton_log.residual_norm[1:e],
#           yscale=:log10, 
#           # lw=lw, 
#           color=:red,
#           linestyle=:dash,
#           marker=:diamond,
#           # markersize=ms, 
#           markerstrokewidth=0,
#           label="No re-initialization (h = 50μs)",
#           )
 
#      e2 = length(sol3.newton_log.residual_norm)

#      plot!(pa, 
#           1:e2, 
#           sol3.newton_log.residual_norm[1:e2],
#           yscale=:log10, 
#           # lw=lw, 
#           color=:teal,
#           marker=:square, 
#           # markersize=ms, 
#           markerstrokewidth=0,
#           label="With re-initialization",
#           )
 
#      hline!(pa, 
#           [tol_line], 
#           color=:black, 
#           linestyle=:dot, 
#           lw=1.6, 
#           label="",
#           )

#      plot!(pa,
#           xlabel="Newton iterations", 
#           ylabel=L"\| F(\mathbf{z}) \|_2",
#           title="", 
#           # titlelocation=:center,
#           # legend=(0.42, 0.40),
#           framestyle=:box,
#           grid=true, 
#           # gridalpha=0.2
#           )
 
#      # ── Right panel (b): homotopy reinitialization detail ──
#      e3 = length(r1.residuals)
#      # pb = plot(1:e3, r1.residuals[1:e3],
#      #      yscale=:log10, lw=lw, color=:orange,
#      #      marker=:circle, markersize=ms, markerstrokewidth=0,
#      #      label="")
 
#      # hline!(pb, [tol_line], color=:black, linestyle=:dot, lw=1.6, label="")
#      # plot!(pb,
#      #      xlabel="Newton iterations", ylabel=L"\| g(\mathbf{z}) \|_2",
#      #      title="", titlelocation=:center,
#      #      framestyle=:box, grid=true,
#      #      # gridalpha=0.2
#      #      )

#      lambdas = r3.λ_hist[2:end]
#      iterations = r3.iter_hist[1:end]
#      pb = plot(
#           lambdas,
#           iterations,
#           marker = :circle,
#       #     markersize = 5,
#       #     linewidth = 2,
#           xlabel = "Homotopy parameter, λ",
#           ylabel = "Newton iterations",
#           xticks = 0:0.2:1.0,
#           yticks = 1:1:7,
#           legend = false,
#           grid = true,
#           framestyle = :box,
#           # size = (700, 400)
#       )
 
#      # ── Combine side-by-side ──
#      fig = plot(pa, 
#                pb, 
#                layout=(1, 2), 
#                size=(1300, 580), 
#                dpi=300,
#                # margin=4Plots.mm, bottom_margin=5Plots.mm
#                )
#      display(fig)
#      # savefig(fig, "39bus16_newton_convergence_comparison_complete.pdf")
#  end

# plot(r3.λ_hist[2:end], r3.iter_hist)

# plot(r3.λ_hist[2:end])

# using Plots

# # Assuming you have collected these arrays during your homotopy loop
# lambdas = r3.λ_hist[3:end]
# iterations = r3.iter_hist[2:end]

# step_indices = 1:length(lambdas)
# # lambdas = [...] 
# # iterations = [...]

# p1 = plot(step_indices, lambdas, 
#     marker=:circle, markersize=3, 
#     label="λ Path", ylabel="λ", color=:blue, lw=2)

# p2 = bar(step_indices, iterations, 
#     label="Newton Iterations", ylabel="Iterations", 
#     xlabel="Continuation Step Index", color=:orange, bar_width=0.5)

# # Combine into a stacked layout
# plot(p1, p2, layout=(2,1), size=(600, 500), link=:x)

# using Plots

# # Mock data
# step_indices = 1:length(lambdas)
# # lambdas = [...]
# # iterations = [...]

# # 1. Plot the primary axis (λ vs Step)
# p = plot(step_indices, lambdas, 
#     label="λ Path", 
#     ylabel="Continuation Parameter (λ)", 
#     xlabel="Continuation Step Index", 
#     color=:blue, lw=2, 
#     yguidefontcolor=:blue, ytickfontcolor=:blue, # Color-match the axis
#     legend=:topleft)

# # 2. Create the twin axis overlay
# p_twin = twinx(p)

# # 3. Plot the secondary data (Iterations vs Step) as a bar chart to avoid line-crossing
# bar!(p_twin, step_indices, iterations, 
#     label="Newton Iterations", 
#     ylabel="Iterations", 
#     color=:orange, alpha=0.6, bar_width=0.5,
#     yguidefontcolor=:orange, ytickfontcolor=:orange, # Color-match the axis
#     legend=:bottomright)

# # Display the combined plot
# display(p)

# plt = plot(lambdas, 
#      iterations, 
#      marker=:circle, 
#      markersize=3,
#      )
# lambdas2 = @. round(lambdas, digits=2)
# xticks!(plt, lambdas2)

# begin
#      default(
#     fontfamily = "Computer Modern",
#     linewidth = 3.5,
#     markersize = 9,
#     legendfontsize = 13,
#     guidefontsize = 14,
#     tickfontsize = 14,
#     titlefontsize = 14,
#     dpi = 300,
#     size = (1200, 500),
#     # margin = 6Plots.m
#     )

#     plt = plot(
#     lambdas,
#     iterations,
#     marker = :circle,
# #     markersize = 5,
# #     linewidth = 2,
#     xlabel = "Homotopy parameter, λ",
#     ylabel = "Newton iterations",
#     xticks = 0:0.2:1.0,
#     yticks = 1:1:7,
#     legend = false,
#     grid = true,
#     framestyle = :box,
#     size = (700, 400)
# )


# end

# plt = plot(
#     lambdas,
#     iterations,
#     marker = :circle,
# #     markersize = 5,
# #     linewidth = 2,
#     xlabel = "Homotopy parameter, λ",
#     ylabel = "Newton iterations",
#     xticks = 0:0.2:1.0,
#     yticks = 1:1:7,
#     legend = false,
#     grid = true,
#     framestyle = :box,
#     size = (700, 400)
# )



using Plots
using LaTeXStrings

begin
    tol_line = 1e-6

    default(
        fontfamily = "Computer Modern",
        linewidth = 3.5,
        markersize = 7,
        markerstrokewidth = 0,
        legendfontsize = 13,
        guidefontsize = 14,
        tickfontsize = 14,
        titlefontsize = 14,
        dpi = 300,
        size = (1000, 380)
    )

    # -------------------------------
    # Left panel (a): Newton residual convergence
    # -------------------------------
    e = 10
    e2 = length(sol3.newton_log.residual_norm)
    xmax_a = max(e, e2)

    pa = plot(
        1:e,
        sol2.newton_log.residual_norm[1:e],
        yscale = :log10,
        color = :blue,
        linestyle = :dash,
        marker = :diamond,
        label = "No re-init. (h = 500 μs)",
    )

    plot!(
        pa,
        1:e,
        sol99.newton_log.residual_norm[1:e],
        color = :red,
        linestyle = :dash,
        marker = :diamond,
        label = "No re-init. (h = 50 μs)",
    )

    plot!(
        pa,
        1:e2,
        sol3.newton_log.residual_norm[1:e2],
        color = :teal,
        linestyle = :solid,
        marker = :square,
        label = "With re-init.",
    )

    hline!(
        pa,
        [tol_line],
        color = :black,
        linestyle = :dot,
        linewidth = 2.5,
        label = "Convergence tol. (1e-6)",
    )

    plot!(
        pa,
        xlabel = "Newton iteration",
        ylabel = L"\|F(\mathbf{z})\|_2",
        framestyle = :box,
        grid = true,
        gridalpha = 0.18,
        minorgrid = false,
     #    legend = :topright,
        legend = (0.42,0.42),
        title = "(a)",
        titlelocation = :left,
        xlims = (0.7, xmax_a + 0.3),   # small padding on both sides
        left_margin = 6Plots.mm,
        bottom_margin = 5Plots.mm,
    )

    # -------------------------------
    # Right panel (b): Homotopy continuation effort
    # -------------------------------
    lambdas = r3.λ_hist[1:end]

    iterations = prepend!(r3.iter_hist[1:end], 1)

    pb = plot(
        lambdas,
        iterations,
        color = :purple,
        linestyle = :solid,
     #    linewidth = 2.5,
        marker = :circle,
        xlabel = L"\lambda",
        ylabel = "Newton iterations",
        xticks = 0:0.2:1.0,
        yticks = minimum(iterations):1:maximum(iterations),
        xlims = (-0.02, 1.03),   # padding so λ=1 marker is fully visible
        ylims = (0.5, maximum(iterations) + 0.5),
        framestyle = :box,
        grid = true,
        gridalpha = 0.18,
        legend = false,
        title = "(b)",
        titlelocation = :left,
        left_margin = 4Plots.mm,
        bottom_margin = 5Plots.mm,
    )

#     vline!(
#         pb,
#         [1.0],
#         color = :black,
#         linestyle = :dot,
#         linewidth = 1.2,
#         label = "",
#     )

    # -------------------------------
    # Combine panels
    # -------------------------------
    fig = plot(
        pa,
        pb,
        layout = @layout([a b]),
        size = (1000, 380),
        dpi = 300,
    )

    display(fig)

    # savefig(fig, "homotopy_newton_results.pdf")
    # savefig(fig, "homotopy_newton_results.png")
end



##########################################################
# save each plot separately


using Plots
using LaTeXStrings

begin
    tol_line = 1e-6
#     tosave = false
     tosave = true
    size = (520, 380)

    default(
        fontfamily = "Computer Modern",
        linewidth = 3.0,
        markersize = 6,
        markerstrokewidth = 0,
        legendfontsize = 13,
        guidefontsize = 14,
        tickfontsize = 14,
        titlefontsize = 13,
        dpi = 300,
    )

    # -------------------------------
    # Figure 1: Newton residual convergence
    # -------------------------------
    e = 10
    e2 = length(sol3.newton_log.residual_norm)
    xmax_a = max(e, e2)

    pa = plot(
        1:e,
        sol2.newton_log.residual_norm[1:e],
        yscale = :log10,
        color = :blue,
        linestyle = :dash,
        marker = :diamond,
        label = "No re-init. (h = 500 μs)",
     #    size = (560, 380),
     #    size = (520, 380),
          size = size,
    )

    plot!(
        pa,
        1:e,
        sol99.newton_log.residual_norm[1:e],
        color = :red,
        linestyle = :dash,
        marker = :diamond,
        label = "No re-init. (h = 50 μs)",
    )

    plot!(
        pa,
        1:e2,
        sol3.newton_log.residual_norm[1:e2],
        color = :teal,
        linestyle = :solid,
        marker = :square,
        label = "With re-init.",
    )

    hline!(
        pa,
        [tol_line],
        color = :black,
        linestyle = :dot,
        linewidth = 2.0,
        label = "Convergence tol.",
    )

    plot!(
        pa,
        xlabel = "NR iteration",
        ylabel = L"\|R(\mathbf{z})\|_2",
     #    ylabel = "|R(z)|_2",
        framestyle = :box,
        grid = true,
        gridalpha = 0.1,
        minorgrid = false,
     #    legend = :topright,
        legend = (0.46, 0.45),
        xlims = (0.7, xmax_a + 0.3),
        left_margin = 1Plots.mm,
        right_margin = 1Plots.mm,
        bottom_margin = 2Plots.mm,
        top_margin = 2Plots.mm,
        title = "",
    )

    display(pa)
    tosave && savefig(pa, "newton_residual_convergence2.pdf")
#     savefig(pa, "newton_residual_convergence.png")

    # -------------------------------
    # Figure 2: Homotopy continuation effort
    # -------------------------------
    lambdas = r3.λ_hist
    iterations = vcat(1, r3.iter_hist)

#     pb = plot(
#         lambdas,
#         iterations,
#         color = :purple,
#         linestyle = :solid,
#         marker = :circle,
#         xlabel = L"\lambda",
#         ylabel = "NR iterations",
#         xticks = 0:0.2:1.0,
#         yticks = minimum(iterations):1:maximum(iterations),
#         xlims = (-0.02, 1.03),
#         ylims = (0.5, maximum(iterations) + 0.5),
#         framestyle = :box,
#         grid = true,
#         gridalpha = 0.1,
#         minorgrid = false,
#         legend = false,
#         left_margin = 2Plots.mm,
#         right_margin = 2Plots.mm,
#         bottom_margin = 2Plots.mm,
#         top_margin = 2Plots.mm,
#         #size = (520, 380),
#           size = size,
#         title = "",
#     )

pb = plot(
    lambdas,
    iterations,
    seriestype = :sticks,
    color = :purple,
    alpha = 0.8,
    linewidth = 2.2,
    xlabel = L"\lambda",
    ylabel = "Newton iterations",
    xticks = 0:0.2:1.0,
    yticks = 0:1:maximum(iterations),
    xlims = (-0.02, 1.02),
    ylims = (0.0, maximum(iterations) + 0.5),
    framestyle = :box,
    grid = true,
    gridalpha = 0.18,
    legend = false,
#     title = "(b)",
    size = size,
    titlelocation = :left,
     left_margin = 2Plots.mm,
     right_margin = 2Plots.mm,
     bottom_margin = 2Plots.mm,
     top_margin = 2Plots.mm,
)

scatter!(
    pb,
    lambdas,
    iterations,
    color = :purple,
    marker = :circle,
    markersize = 7,
    label = "",
)

    display(pb)
    tosave && savefig(pb, "homotopy_continuation_effort2.pdf")
#     savefig(pb, "homotopy_continuation_effort.png")
end



###########################
# one plot

##########################################################
# save each plot separately


using Plots
using LaTeXStrings

begin
    tol_line = 1e-6
#     tosave = false
     tosave = true
    size = (1200, 500)

    default(
        fontfamily = "Computer Modern",
        linewidth = 3.5,
        markersize = 9,
        markerstrokewidth = 0,
        legendfontsize = 15,
        guidefontsize = 15,
        tickfontsize = 15,
        titlefontsize = 15,
        dpi = 300,
        size = size,
        # margin = 6Plots.mm
        left_margin = 6Plots.mm,
        right_margin = 2Plots.mm,
        bottom_margin = 7Plots.mm,
        top_margin = 2Plots.mm,
    )

    # -------------------------------
    # Figure 1: Newton residual convergence
    # -------------------------------
    e = 10
    e2 = length(sol3.newton_log.residual_norm)
    xmax_a = max(e, e2)

    pa = plot(
        1:e,
        sol2.newton_log.residual_norm[1:e],
        yscale = :log10,
        color = :blue,
        linestyle = :dash,
        marker = :diamond,
        label = "No re-init. (h = 500 μs)",
     #    size = (560, 380),
     #    size = (520, 380),
          # size = size,
    )

    plot!(
        pa,
        1:e,
        sol99.newton_log.residual_norm[1:e],
        color = :red,
        linestyle = :dash,
        marker = :diamond,
        label = "No re-init. (h = 50 μs)",
    )

    plot!(
        pa,
        1:e2,
        sol3.newton_log.residual_norm[1:e2],
        color = :teal,
        linestyle = :solid,
        marker = :square,
        label = "With re-init.",
    )

    hline!(
        pa,
        [tol_line],
        color = :black,
        linestyle = :dot,
        linewidth = 2.0,
        label = "Convergence tol.",
    )

    plot!(
        pa,
        xlabel = "NR iteration",
        ylabel = L"\|R(\mathbf{z})\|_2",
     #    ylabel = "|R(z)|_2",
        framestyle = :box,
        grid = true,
        gridalpha = 0.1,
        minorgrid = false,
     #    legend = :topright,
        legend = (0.5, 0.45),
        xlims = (0.7, xmax_a + 0.3),
     #    left_margin = 1Plots.mm,
     #    right_margin = 1Plots.mm,
     #    bottom_margin = 2Plots.mm,
     #    top_margin = 2Plots.mm,
        title = "(a)",
    )

#     display(pa)
#     tosave && savefig(pa, "newton_residual_convergence2.pdf")
#     savefig(pa, "newton_residual_convergence.png")

    # -------------------------------
    # Figure 2: Homotopy continuation effort
    # -------------------------------
    lambdas = r3.λ_hist
    iterations = vcat(1, r3.iter_hist)

#     pb = plot(
#         lambdas,
#         iterations,
#         color = :purple,
#         linestyle = :solid,
#         marker = :circle,
#         xlabel = L"\lambda",
#         ylabel = "NR iterations",
#         xticks = 0:0.2:1.0,
#         yticks = minimum(iterations):1:maximum(iterations),
#         xlims = (-0.02, 1.03),
#         ylims = (0.5, maximum(iterations) + 0.5),
#         framestyle = :box,
#         grid = true,
#         gridalpha = 0.1,
#         minorgrid = false,
#         legend = false,
#         left_margin = 2Plots.mm,
#         right_margin = 2Plots.mm,
#         bottom_margin = 2Plots.mm,
#         top_margin = 2Plots.mm,
#         #size = (520, 380),
#           size = size,
#         title = "",
#     )

pb = plot(
    lambdas,
    iterations,
    seriestype = :sticks,
    color = :purple,
    alpha = 0.8,
#     linewidth = 2.2,
    xlabel = L"\lambda",
    ylabel = "Newton iterations",
    xticks = 0:0.2:1.0,
    yticks = 0:1:maximum(iterations),
    xlims = (-0.02, 1.02),
    ylims = (0.0, maximum(iterations) + 0.5),
    framestyle = :box,
    grid = true,
    gridalpha = 0.1,
    legend = false,
    title = "(b)",
    size = size,
#     titlelocation = :left,
     # left_margin = 2Plots.mm,
     # right_margin = 2Plots.mm,
     # bottom_margin = 2Plots.mm,
     # top_margin = 2Plots.mm,
)

scatter!(
    pb,
    lambdas,
    iterations,
    color = :purple,
    marker = :circle,
#     markersize = 7,
    label = "",
)

plt = plot(
    pa, pb,
    layout = (1, 2),
    size = (1300, 550)
)

    display(plt)
    tosave && savefig(plt, "homotopy_continuation_effort_combined.pdf")
#     savefig(pb, "homotopy_continuation_effort.png")
end