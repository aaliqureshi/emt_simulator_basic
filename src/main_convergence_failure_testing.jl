using Pkg; Pkg.activate(".")
using Revise
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

models.fault.bus = [20]
models.fault.x_fault[1] = 0.015
# models.fault.x_fault[1] = 0.038

# models.fault.bus = [16]
# models.fault.x_fault[1] = 0.015

# models.fault.bus=[3]
# models.fault.x_fault[1] = 0.00015
# models.load.p .*= 4.03

# 2. Solve power flow
solve_power_flow!(sys)

# 3. Static initialization
run_static_init!(sys)

# 4. Dynamic simulation setup
address = build_dynamic_address(sys)
mass_matrix = build_mass_matrix(sys, address)
u0 = build_initial_conditions(sys, address)


function run_simulation(u0, lambda, dt; always_new=true, tstops=[], adaptive=false, t_end=0.1)
    p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
    prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, t_end), p, mass_matrix)
    sol = MyDiffEq.Solve(prob, dt, method=:Euler, adaptive=adaptive, tstops=tstops, always_new=always_new)
    return sol
end

# Pre-fault simulation
lambda=0.0
dt0=5e-4
sol_pf=run_simulation(u0, lambda, dt0, t_end=0.01)
u1 = sol_pf.u[end]

# fault-on simulation 3
# dt_list = [5e-4, 5e-5, 5e-6]
dt_list = [5e-4, 5e-5]
# dt_list=[5e-4]
lambda = 1.0
sol_list = []
num_fails = 0
for iter in eachindex(dt_list)
    ux = run_simulation(u1, lambda, dt_list[iter], t_end=0.01)
    if ux.retcode == :MaxIter
        num_fails+=1
    end
    push!(sol_list, ux)
    # println("For dt=$(dt_list[iter]), retcode = $(ux.retcode)")
end
@show num_fails

# fault-on flat initialization
u2=copy(u1)
u2[address["balance_d"]] .= 1.0
u2[address["balance_q"]] .= 0.0 
dt_flat = 5e-4
sol_flat = run_simulation(u2, lambda, dt_flat)


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



# 5. post re-init simulation
lambda = 0.0
# u3 = copy(u0_new)
u3= copy(u0_homotopy)
# u3= copy(u1)
dt_post = 5e-4
sol_post = run_simulation(u3, lambda, dt_post, t_end=dt_post)


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
    e2 = length(sol_post.newton_log.residual_norm)
    xmax_a = max(e, e2)



    pa = plot(
        1:e,
        sol_list[1].newton_log.residual_norm[1:e],
        yscale = :log10,
        color = :blue,
        linestyle = :dash,
        marker = :diamond,
        # label = "No re-init. (h = 500 μs)",
        label = "No re-init. (h = $(trunc(Int, dt_list[1]/1e-6)) μs)"
    )

    # color_list = [:red, ]

    for iter in collect(2:length(sol_list))
        plot!(
            pa,
            1:e,
            sol_list[iter].newton_log.residual_norm[1:e],
            # color = :red,
            palette = :Dark2,
            linestyle = :dash,
            marker = :diamond,
            # label = "No re-init. (h = 50 μs)",
            label = "No re-init. (h = $(trunc(Int, dt_list[iter]/1e-6)) μs)",
        )
    end

    plot!(
        pa,
        1:e,
        sol_flat.newton_log.residual_norm[1:e],
        # color = :teal,
        palette = :Dark2,
        linestyle = :dash,
        marker = :diamond,
        label = "Flat re-init (h = $(trunc(Int, dt_flat/1e-6)) μs).",
    )

    plot!(
        pa,
        1:e2,
        sol_post.newton_log.residual_norm[1:e2],
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

begin
    tol_line = 1e-6
    tosave = false
    #  tosave = true
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
    e2 = length(sol_post.newton_log.residual_norm)
    xmax_a = max(e, e2)

    pa = plot(
        1:e,
        sol_list[1].newton_log.residual_norm[1:e],
        yscale = :log10,
        color = :blue,
        linestyle = :dash,
        marker = :diamond,
        # label = "No re-init. (h = 500 μs)",
        label = "No re-init. (h = $(trunc(Int, dt_list[1]/1e-6)) μs)",
     #    size = (560, 380),
     #    size = (520, 380),
          size = size,
    )

    for iter in collect(2:length(sol_list))
        plot!(
            pa,
            1:e,
            sol_list[iter].newton_log.residual_norm[1:e],
            # color = :red,
            palette = :Dark2,
            linestyle = :dash,
            marker = :diamond,
            # label = "No re-init. (h = 50 μs)",
            label = "No re-init. (h = $(trunc(Int, dt_list[iter]/1e-6)) μs)",
        )
    end

    plot!(
        pa,
        1:e,
        sol_flat.newton_log.residual_norm[1:e],
        # color = :teal,
        palette = :Dark2,
        linestyle = :dash,
        marker = :diamond,
        label = "Flat re-init (h = $(trunc(Int, dt_flat/1e-6)) μs).",
    )

    plot!(
        pa,
        1:e2,
        sol_post.newton_log.residual_norm[1:e2],
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
        label = "Convergence tol. (1e-6)",
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
    tosave = false
    #  tosave = true
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
    e2 = length(sol_post.newton_log.residual_norm)
    xmax_a = max(e, e2)

    pa = plot(
        1:e,
        sol_list[1].newton_log.residual_norm[1:e],
        yscale = :log10,
        color = :blue,
        linestyle = :dash,
        marker = :diamond,
        # label = "No re-init. (h = 500 μs)",
        label = "No re-init. (h = $(trunc(Int, dt_list[1]/1e-6)) μs)",
     #    size = (560, 380),
     #    size = (520, 380),
          # size = size,
    )

    for iter in collect(2:length(sol_list))
        plot!(
            pa,
            1:e,
            sol_list[iter].newton_log.residual_norm[1:e],
            # color = :red,
            palette = :Dark2,
            linestyle = :dash,
            marker = :diamond,
            # label = "No re-init. (h = 50 μs)",
            label = "No re-init. (h = $(trunc(Int, dt_list[iter]/1e-6)) μs)",
        )
    end

    plot!(
        pa,
        1:e,
        sol_flat.newton_log.residual_norm[1:e],
        # color = :teal,
        palette = :Dark2,
        linestyle = :dash,
        marker = :diamond,
        label = "Flat re-init (h = $(trunc(Int, dt_flat/1e-6)) μs).",
    )

    plot!(
        pa,
        1:e2,
        sol_post.newton_log.residual_norm[1:e2],
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
    folder = "/Users/aali27/Work/repos/emt_simulator_basic/figures/review1/"
    tosave && savefig(plt, folder*"homotopy_continuation_effort_combined.pdf")
end