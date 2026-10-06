using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Revise
using Plots, LaTeXStrings
using Barq
includet("SimulateDivergence.jl")


begin
    size = (1200, 500)
    default(
    fontfamily = "Computer Modern",
    linewidth = 3.5,
    markersize = 8,
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
end

# run complete simulation
begin
    case = 39
    # case = 118
    method = :Trap
    # method = :Euler
    # dt_list = [5e-4, 5e-5]
    dt_list = [5e-4, 5e-5, 5e-6]
    sim_data = Simulate.run_simulation(case_bus=case, method=method, dt_list=dt_list);
end

begin
    # case = 39
    case = 118
    method = :Trap
    # method = :Euler
    # dt_list = [5e-4, 5e-5]
    sim_data_trap = Simulate.run_simulation(case_bus=case, method=method, dt_list=dt_list);
end

## plot for case 39 and case 118 convergence
begin
    const C_BLUE   = "#0072B2"
    const C_ORANGE = "#D55E00"
    const C_GREEN  = "#009E73"
    const C_SKY    = "#56B4E9"
    const C_PURPLE = "#CC79A7"
    const C_BLACK  = "#333333"

    tol_line = 1e-6
    lw = 3.5
    tosave = false
    # tosave = true

    # idx_div = length(sim_data.sim_list[1].newton_log.residual)
    idx_div = 10
    idx_reinit = length(sim_data.sol_post.newton_log.residual_norm)

    xmax_a = max(idx_div, idx_reinit)

    pa = plot(1:idx_div,
              sim_data.sol_list[1].newton_log.residual_norm[1:idx_div],
              yscale=:log10,
            #   color=:blue,
              color=C_BLUE,
              linestyle=:solid,
              marker=:circle,
              label="No re-init (h = $(trunc(Int, dt_list[1]/1e-6)) μs)",
              linewidth=lw,
              markerstrokecolor = :black,
              markerstrokewidth = 1.0,
              )
    
    no_reinit_color = [C_ORANGE, C_PURPLE]
    no_reinit_marker = [:diamond, :utriangle]
    no_reinit_style = [:dash, :dashdot]
    
    for iter in collect(2:length(dt_list))
        plot!(
            pa,
            1:idx_div,
            sim_data.sol_list[iter].newton_log.residual_norm[1:idx_div],
            # color = :red,
            yscale=:log10,
            # palette = :Dark2,
            color=no_reinit_color[iter-1],
            # linestyle = :dash,
            linestyle=no_reinit_style[iter-1],
            # marker = :diamond,
            marker=no_reinit_marker[iter-1],
            # label = "No re-init. (h = 50 μs)",
            label = "No re-init. (h = $(trunc(Int, dt_list[iter]/1e-6)) μs)",
            linewidth=lw,
            markerstrokecolor = :black,
            markerstrokewidth = 1.0,
        )
    end

    plot!(
        pa,
        1:idx_div,
        sim_data_trap.sol_flat.newton_log.residual_norm[1:idx_div],
        # color = :teal,
        # palette = :Dark2,
        color = C_SKY,
        # linestyle = :dash,
        linestyle=:dashdotdot,
        # marker = :diamond,
        marker=:dtriangle,
        label = "Flat start re-init. (h = $(trunc(Int, sim_data.sol_flat.dt/1e-6)) μs)",
        linewidth=lw,
        markerstrokecolor = :black,
        markerstrokewidth = 1.0,
    )

    plot!(
        pa,
        1:idx_reinit,
        # sim_data.sol_post.newton_log.residual_norm[1:idx_reinit],
        sim_data.sol_post.newton_log.residual_norm,
        # color = :teal,
        color=C_GREEN,
        linestyle = :solid,
        marker = :square,
        label = "With re-init. (h = $(trunc(Int, sim_data.sol_post.dt/1e-6)) μs)",
        linewidth=lw,
        markerstrokecolor = :black,
        markerstrokewidth = 1.0,
    )

    hline!(
        pa,
        [tol_line],
        # color = :black,
        color=C_BLACK,
        linestyle = :dot,
        linewidth = 1.5,
        label = "Convergence tolerance",
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
        # legend = :bottomleft,
        # legend = (0.5, 0.45),
        legend=false,
        xlims = (0.7, xmax_a + 0.3),
     #    left_margin = 1Plots.mm,
     #    right_margin = 1Plots.mm,
     #    bottom_margin = 2Plots.mm,
     #    top_margin = 2Plots.mm,
        title = "(a) IEEE 39-bus",
    )
    # display(pa)
end

begin
    # -------------------------------
    # Figure 2: Case 118
    # -------------------------------

    pb = plot(1:idx_div,
              sim_data_trap.sol_list[1].newton_log.residual_norm[1:idx_div],
              yscale=:log10,
            #   color=:blue,
            #   linestyle=:dash,
            #   marker=:diamond,
              color=C_BLUE,
              linestyle=:solid,
              marker=:circle,
              label="No re-init (h = $(trunc(Int, dt_list[1]/1e-6)) μs)",
              linewidth=lw,
              markerstrokecolor = :black,
              markerstrokewidth = 1.0,
              )
    
    for iter in collect(2:length(dt_list))
        plot!(
            pb,
            1:idx_div,
            sim_data_trap.sol_list[iter].newton_log.residual_norm[1:idx_div],
            # color = :red,
            yscale=:log10,
            # palette = :Dark2,
            # linestyle = :dash,
            # marker = :diamond,
            color=no_reinit_color[iter-1],
            linestyle=no_reinit_style[iter-1],
            marker=no_reinit_marker[iter-1],
            # label = "No re-init. (h = 50 μs)",
            label = "No re-init., (h = $(trunc(Int, dt_list[iter]/1e-6)) μs)",
            linewidth=lw,
            markerstrokecolor = :black,
            markerstrokewidth = 1.0,
        )
    end

    plot!(
        pb,
        1:idx_div,
        sim_data_trap.sol_flat.newton_log.residual_norm[1:idx_div],
        # color = :teal,
        # palette = :Dark2,
        # linestyle = :dash,
        # marker = :diamond,
        color = C_SKY,
        linestyle=:dashdotdot,
        marker=:dtriangle,
        label = "Flat start re-init. (h = $(trunc(Int, sim_data_trap.sol_flat.dt/1e-6)) μs)",
        linewidth=lw,
        markerstrokecolor = :black,
        markerstrokewidth = 1.0,
    )

    plot!(
        pb,
        1:length(sim_data_trap.sol_post.newton_log.residual_norm),
        # sim_data_trap.sol_post.newton_log.residual_norm[1:idx_reinit],
        sim_data_trap.sol_post.newton_log.residual_norm,
        # color = :teal,
        # linestyle = :solid,
        # marker = :square,
        color=C_GREEN,
        linestyle = :solid,
        marker = :square,
        label = "HC re-init. (h = $(trunc(Int, sim_data_trap.sol_post.dt/1e-6)) μs)",
        linewidth=lw,
        markerstrokecolor = :black,
        markerstrokewidth = 1.0,
    )

    hline!(
        pb,
        [tol_line],
        # color = :black,
        color=C_BLACK,
        linestyle = :dot,
        linewidth = 1.5,
        label = "Convergence tol.",
    )

    plot!(
        pb,
        xlabel = "NR iteration",
        ylabel = L"\|R(\mathbf{z})\|_2",
     #    ylabel = "|R(z)|_2",
        framestyle = :box,
        grid = true,
        gridalpha = 0.1,
        minorgrid = false,
        legend = :bottomright,
        # legend = (0.5, 0.45),
        # legend=false,
        xlims = (0.7, xmax_a + 0.3),
     #    left_margin = 1Plots.mm,
     #    right_margin = 1Plots.mm,
     #    bottom_margin = 2Plots.mm,
     #    top_margin = 2Plots.mm,
        title = "(b) IEEE 118-bus",
    )
end

begin
    plt = plot(
        pa, pb,
        layout = (1, 2),
        # size = (1300, 550),
        size=(1300,400),
        # linewidth=3.75,
        tickfontsize = 16,  # Bump up fonts so they are legible
        guidefontsize = 16,
        legendfontsize=13,
        left_margin = 8Plots.mm,
        right_margin = 1Plots.mm,
        bottom_margin = 8Plots.mm,
        top_margin = 2Plots.mm,
        grid=false,
        # gridalpha=0.05,
        )
    
    display(plt)
    folder = "/Users/aali27/Work/repos/emt_simulator_basic/figures/review1/"
    tosave && savefig(plt, folder*"nr_divergence_new.pdf")
end


## plot for continuation effort
# lambdas_39 = sim_data.reinit.r3.λ_hist
# lambdas_118 = sim_data_trap.reinit.r3.λ_hist
# iterations_39 = vcat(1, sim_data.reinit.r3.iter_hist)
# iterations_118 = vcat(1, sim_data_trap.reinit.r3.iter_hist)


c39  = "#0072B2"
c118 = "#D55E00"
lw = 3.5
# ## plot voltages
begin
    # tosave = true
    tosave = false
    
    bus_118 = 12
    bus_39 = 24

    sol39 = sim_data
    sol118 = sim_data_trap

    vd_118_pre = [u[sol118.address["balance_d"]][bus_118] for u in sol118.sol_pf.u]
    vq_118_pre = [u[sol118.address["balance_q"]][bus_118] for u in sol118.sol_pf.u]
    vd_39_pre = [u[sol39.address["balance_d"]][bus_39] for u in sol39.sol_pf.u]
    vq_39_pre = [u[sol39.address["balance_q"]][bus_39] for u in sol39.sol_pf.u]

    vd_118_post = [u[sol118.address["balance_d"]][bus_118] for u in sol118.sol_post.u]
    vq_118_post = [u[sol118.address["balance_q"]][bus_118] for u in sol118.sol_post.u]
    vd_39_post = [u[sol39.address["balance_d"]][bus_39] for u in sol39.sol_post.u]
    vq_39_post = [u[sol39.address["balance_q"]][bus_39] for u in sol39.sol_post.u]

    v_39_pre = @. abs(vd_39_pre + 1im*vq_39_pre)
    v_118_pre = @. abs(vd_118_pre + 1im*vq_118_pre)

    v_39_post = @. abs(vd_39_post + 1im*vq_39_post)
    v_118_post = @. abs(vd_118_post + 1im*vq_118_post)

    v_39 = vcat(v_39_pre, v_39_post)
    v_118 = vcat(v_118_pre, v_118_post)

    t_pre39 = sol39.sol_pf.time
    t_post39 = @. t_pre39[end] + sol39.sol_post.time
    t39 = vcat(t_pre39, t_post39)

    t_pre118 = sol118.sol_pf.time
    t_post118 = @. t_pre118[end] + sol118.sol_post.time
    t118 = vcat(t_pre118, t_post118)
    # dt = 5e-4

    t39 = t39/1e-3
    t118 = t118/1e-3

    idx_end = 16


    
    p1 = plot(
        t39[1:idx_end],
        v_39[1:idx_end],
        # color = :blue,
        color=c39,
        linestyle=:solid,
        label = "39-bus",
        linewidth=lw,
    )
    plot!(p1,
        t118[1:idx_end],
        v_118[1:idx_end],
        # color = :orange,
        color=c118,
        linestyle=:dash,
        label = "118-bus",
        linewidth=lw,
    )

    vline!(
        p1,
        [5.0],
        color = :black,
        linestyle = :dot,
        linewidth = 1.5,
        label = "",
    )
    # annotate!(5.2, 0.95, text("Fault", 8))


    plot!(
        p1,
        xlabel = "Time (ms)",
        ylabel = L"|V_f| \ \mathrm{(p.u.)}",
     #    ylabel = "|R(z)|_2",
        framestyle = :box,
        grid = true,
        gridalpha = 0.1,
        minorgrid = false,
        legend = :topright,
        # legend = (0.5, 0.45),
        # legend=false,
        # xlims = (0.7, xmax_a + 0.3),
        ylims = (0.0, 1.2),
        xticks = 0.0:5.0:30.0,
     #    left_margin = 1Plots.mm,
     #    right_margin = 1Plots.mm,
     #    bottom_margin = 2Plots.mm,
     #    top_margin = 2Plots.mm,
        title = "(b) Time-domain voltage",
    )
    display(p1)

    folder = "/Users/aali27/Work/repos/emt_simulator_basic/figures/review1/"
    tosave && savefig(p1, folder*"voltage_plot.pdf")

end

## plot continuation path of voltages
# begin
#     # tosave = true
#     tosave = false
    
#     t_fine = 0.0:0.01:1.0
#     t_adap_39 = sim_data.reinit.r3.λ_hist
#     t_adap_118 = sim_data_trap.reinit.r3.λ_hist

#     v_fine_39 = @. abs(sim_data.reinit.r2.vd_hist + 1im*sim_data.reinit.r2.vq_hist)
#     v_fine_118 = @. abs(sim_data_trap.reinit.r2.vd_hist + 1im*sim_data_trap.reinit.r2.vq_hist)

#     v_adap_39 = @. abs(sim_data.reinit.r3.vd_hist + 1im*sim_data.reinit.r3.vq_hist)
#     v_adap_118 = @. abs(sim_data_trap.reinit.r3.vd_hist + 1im*sim_data_trap.reinit.r3.vq_hist)
    
    
#     p1 = plot(
#         t_fine,
#         v_fine_39,
#         # label = "IEEE 39-bus",
#         label=""
#     )
#     scatter!(p1,
#         t_adap_39,
#         v_adap_39,
#         label=""
#     )
#     plot!(
#         p1,
#         xlabel = "Continuation Parameter (λ)",
#         ylabel = "|V| (p.u.)",
#      #    ylabel = "|R(z)|_2",
#         framestyle = :box,
#         grid = true,
#         gridalpha = 0.1,
#         minorgrid = false,
#         legend = :bottomleft,
#         # ylims = (0.0, 1.1),
#         title = "(a) IEEE 39-bus",
#     )

#     p2 = plot(
#         t_fine,
#         v_fine_118,
#         # label = "IEEE 39-bus",
#         label=""
#     )
#     scatter!(p2,
#         t_adap_118,
#         v_adap_118,
#         label=""
#     )
#     plot!(
#         p2,
#         xlabel = "Continuation Parameter (λ)",
#         ylabel = "|V| (p.u.)",
#      #    ylabel = "|R(z)|_2",
#         framestyle = :box,
#         grid = true,
#         gridalpha = 0.1,
#         minorgrid = false,
#         legend = :bottomleft,
#         # ylims = (0.0, 1.1),
#         title = "(b) IEEE 118-bus",
#     )
#     # display(p1)
#     plt = plot(
#             p1, p2,
#             layout = (1, 2),
#             size = (1300, 550),
#             )

#     display(plt)

#     folder = "/Users/aali27/Work/repos/emt_simulator_basic/figures/review1/"
#     tosave && savefig(plt, folder*"voltage_plot.pdf")

# end

begin
    # tosave = true
    tosave = false
    
    t_fine = 0.0:0.01:1.0
    t_adap_39 = sim_data.reinit.r3.λ_hist
    t_adap_118 = sim_data_trap.reinit.r3.λ_hist

    v_fine_39 = @. abs(sim_data.reinit.r2.vd_hist + 1im*sim_data.reinit.r2.vq_hist)
    v_fine_118 = @. abs(sim_data_trap.reinit.r2.vd_hist + 1im*sim_data_trap.reinit.r2.vq_hist)

    v_adap_39 = @. abs(sim_data.reinit.r3.vd_hist + 1im*sim_data.reinit.r3.vq_hist)
    v_adap_118 = @. abs(sim_data_trap.reinit.r3.vd_hist + 1im*sim_data_trap.reinit.r3.vq_hist)

    common_args = (
        xlabel = "Continuation parameter (λ)",
        ylabel = L"|V_f| \ \mathrm{(p.u.)}",
        # xlims = (-0.0, 1.05),
        ylims = (0.0, 1.2),
        xticks = 0.0:0.25:1.0,
        # yticks = 0.0:0.25:1.0,
        framestyle = :box,
        grid = true,
        gridalpha = 0.1,
        minorgrid = false,
        linewidth = 3.5,
    )
    
    
    p2 = plot(
        t_fine,
        v_fine_39;
        label = "39-bus",
        common_args...,
        # label="",
        # title = "(a) IEEE 39-bus",
        color=c39,
        linewidth=lw,
        # common_args...
        # yscale=:log10,
    )
    scatter!(p2,
        t_adap_39,
        v_adap_39,
        label="",
        markersize = 8.0,
        # markerstrokewidth = 0.5,
        color=c39,
        # yscale=:log10,
        markerstrokecolor = :black,
        markerstrokewidth = 1.0,
    )

    # plot!(
    #     p1,
    #     xlabel = "Continuation Parameter (λ)",
    #     ylabel = "|V| (p.u.)",
    #  #    ylabel = "|R(z)|_2",
    #     framestyle = :box,
    #     grid = true,
    #     gridalpha = 0.1,
    #     minorgrid = false,
    #     legend = :bottomleft,
    #     # ylims = (0.0, 1.1),
    #     title = "(a) IEEE 39-bus",
    # )

    # p2 = plot(
    #     t_fine,
    #     v_fine_118,
    #     # label = "IEEE 39-bus",
    #     label = "Fine-step continuation",
    #     title = "(b) IEEE 118-bus";
    #     common_args...
    # )
    plot!(
        p2,
        t_fine,
        v_fine_118,
        # label = "IEEE 39-bus",
        label = "118-bus";
        # title = "(b) IEEE 118-bus";
        common_args...,
        color=c118,
        linestyle=:dash,
        linewidth=lw,
        # yscale=:log10,
    )
    # scatter!(p2,
    #     t_adap_118,
    #     v_adap_118,
    #     # label="",
    #     label = "Adaptive HC points",
    #     markersize = 9.0,
    #     markerstrokewidth = 0.5,
    # )
    scatter!(p2,
        t_adap_118,
        v_adap_118,
        # label="",
        label = "",
        marker=:diamond,
        markersize = 8.0,
        # markerstrokewidth = 0.5,
        color=c118,
        title="(a) Continuation path",
        markerstrokecolor = :black,
        markerstrokewidth = 1.0,
        )
    # plot!(
    #     p2,
    #     legend = :bottomleft,
    # )
    # display(p1)
    # plt = plot(
    #         p1, p2,
    #         layout = (1, 2),
    #         size = (1300, 550),
    #         )

    # display(plt)
    # display(p2)

    folder = "/Users/aali27/Work/repos/emt_simulator_basic/figures/review1/"
    tosave && savefig(p2, folder*"voltage_continuation.pdf")

end

# combine figures here
begin
    # tosave=false
    tosave=true
    plt2 = plot(
        p2, p1,
        layout=(1,2),
        size=(1300,400),
        # linewidth=3.75,
        tickfontsize = 16,  # Bump up fonts so they are legible
        guidefontsize = 16,
        legendfontsize=14,
        left_margin = 6Plots.mm,
        right_margin = 2Plots.mm,
        bottom_margin = 9Plots.mm,
        top_margin = 2Plots.mm,
        grid=false,
        )
    folder = "/Users/aali27/Work/repos/emt_simulator_basic/figures/review1/"
    tosave && savefig(plt2, folder*"voltage_combined.pdf")
end