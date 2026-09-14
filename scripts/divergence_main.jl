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
    tol_line = 1e-6
    tosave = false
    # tosave = true

    # idx_div = length(sim_data.sim_list[1].newton_log.residual)
    idx_div = 10
    idx_reinit = length(sim_data.sol_post.newton_log.residual_norm)

    xmax_a = max(idx_div, idx_reinit)

    pa = plot(1:idx_div,
              sim_data.sol_list[1].newton_log.residual_norm[1:idx_div],
              yscale=:log10,
              color=:blue,
              linestyle=:dash,
              marker=:diamond,
              label="No re-init (h = $(trunc(Int, dt_list[1]/1e-6)) μs)",
              )
    
    for iter in collect(2:length(dt_list))
        plot!(
            pa,
            1:idx_div,
            sim_data.sol_list[iter].newton_log.residual_norm[1:idx_div],
            # color = :red,
            yscale=:log10,
            palette = :Dark2,
            linestyle = :dash,
            marker = :diamond,
            # label = "No re-init. (h = 50 μs)",
            label = "No re-init. (h = $(trunc(Int, dt_list[iter]/1e-6)) μs)",
        )
    end

    plot!(
        pa,
        1:idx_div,
        sim_data_trap.sol_flat.newton_log.residual_norm[1:idx_div],
        # color = :teal,
        palette = :Dark2,
        linestyle = :dash,
        marker = :diamond,
        label = "Flat start re-init. (h = $(trunc(Int, sim_data.sol_flat.dt/1e-6)) μs)",
    )

    plot!(
        pa,
        1:idx_reinit,
        # sim_data.sol_post.newton_log.residual_norm[1:idx_reinit],
        sim_data.sol_post.newton_log.residual_norm,
        color = :teal,
        linestyle = :solid,
        marker = :square,
        label = "With re-init. (h = $(trunc(Int, sim_data.sol_post.dt/1e-6)) μs)",
    )

    hline!(
        pa,
        [tol_line],
        color = :black,
        linestyle = :dot,
        linewidth = 2.0,
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
              color=:blue,
              linestyle=:dash,
              marker=:diamond,
              label="No re-init (h = $(trunc(Int, dt_list[1]/1e-6)) μs)",
              )
    
    for iter in collect(2:length(dt_list))
        plot!(
            pb,
            1:idx_div,
            sim_data_trap.sol_list[iter].newton_log.residual_norm[1:idx_div],
            # color = :red,
            yscale=:log10,
            palette = :Dark2,
            linestyle = :dash,
            marker = :diamond,
            # label = "No re-init. (h = 50 μs)",
            label = "No re-init. (h = $(trunc(Int, dt_list[iter]/1e-6)) μs)",
        )
    end

    plot!(
        pb,
        1:idx_div,
        sim_data_trap.sol_flat.newton_log.residual_norm[1:idx_div],
        # color = :teal,
        palette = :Dark2,
        linestyle = :dash,
        marker = :diamond,
        label = "Flat start re-init. (h = $(trunc(Int, sim_data_trap.sol_flat.dt/1e-6)) μs)",
    )

    plot!(
        pb,
        1:length(sim_data_trap.sol_post.newton_log.residual_norm),
        # sim_data_trap.sol_post.newton_log.residual_norm[1:idx_reinit],
        sim_data_trap.sol_post.newton_log.residual_norm,
        color = :teal,
        linestyle = :solid,
        marker = :square,
        label = "With re-init. (h = $(trunc(Int, sim_data_trap.sol_post.dt/1e-6)) μs)",
    )

    hline!(
        pb,
        [tol_line],
        color = :black,
        linestyle = :dot,
        linewidth = 2.0,
        label = "Convergence tolerance",
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
        size = (1300, 550),
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



# ## plot voltages
# begin
#     # tosave = true
#     tosave = false
    
#     bus_118 = 12
#     bus_39 = 24

#     vd_118_pre = [u[sim_data_trap.address["balance_d"]][bus_118] for u in sim_data_trap.sol_pf.u]
#     vq_118_pre = [u[sim_data_trap.address["balance_q"]][bus_118] for u in sim_data_trap.sol_pf.u]
#     vd_39_pre = [u[sim_data.address["balance_d"]][bus_39] for u in sim_data.sol_pf.u]
#     vq_39_pre = [u[sim_data.address["balance_q"]][bus_39] for u in sim_data.sol_pf.u]

#     vd_118_post = [u[sim_data_trap.address["balance_d"]][bus_118] for u in sim_data_trap.sol_post.u]
#     vq_118_post = [u[sim_data_trap.address["balance_q"]][bus_118] for u in sim_data_trap.sol_post.u]
#     vd_39_post = [u[sim_data.address["balance_d"]][bus_39] for u in sim_data.sol_post.u]
#     vq_39_post = [u[sim_data.address["balance_q"]][bus_39] for u in sim_data.sol_post.u]

#     v_39_pre = @. abs(vd_39_pre + 1im*vq_39_pre)
#     v_118_pre = @. abs(vd_118_pre + 1im*vq_118_pre)

#     v_39_post = @. abs(vd_39_post + 1im*vq_39_post)
#     v_118_post = @. abs(vd_118_post + 1im*vq_118_post)

#     v_39 = vcat(v_39_pre, v_39_post)
#     v_118 = vcat(v_118_pre, v_118_post)

#     t_pre = sim_data.sol_pf.time
#     t_post = @. t_pre[end] + sim_data.sol_post.time
#     t = vcat(t_pre, t_post)
#     dt = 5e-4
    
#     p1 = plot(
#         t,
#         v_39,
#         # color = :blue,
#         label = "IEEE 39-bus",
#     )
#     plot!(p1,
#         t,
#         v_118,
#         # color = :orange,
#         label = "IEEE 118-bus",
#     )
#     plot!(
#         p1,
#         xlabel = "Time (sec)",
#         ylabel = "V (p.u.)",
#      #    ylabel = "|R(z)|_2",
#         framestyle = :box,
#         grid = true,
#         gridalpha = 0.1,
#         minorgrid = false,
#         legend = :topright,
#         # legend = (0.5, 0.45),
#         # legend=false,
#         # xlims = (0.7, xmax_a + 0.3),
#         ylims = (0.0, 1.1),
#      #    left_margin = 1Plots.mm,
#      #    right_margin = 1Plots.mm,
#      #    bottom_margin = 2Plots.mm,
#      #    top_margin = 2Plots.mm,
#         # title = "(a)",
#     )
#     display(p1)

#     folder = "/Users/aali27/Work/repos/emt_simulator_basic/figures/review1/"
#     tosave && savefig(plt, folder*"voltage_plot.pdf")

# end

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
        ylabel = "|V| (p.u.)",
        xlims = (-0.0, 1.05),
        ylims = (0.0, 1.10),
        # xticks = 0.0:0.25:1.0,
        yticks = 0.0:0.25:1.0,
        framestyle = :box,
        grid = true,
        gridalpha = 0.1,
        minorgrid = false,
        linewidth = 2.5,
    )
    
    
    p1 = plot(
        t_fine,
        v_fine_39;
        # label = "IEEE 39-bus",
        common_args...,
        label="",
        title = "(a) IEEE 39-bus",
        # common_args...
    )
    scatter!(p1,
        t_adap_39,
        v_adap_39,
        label="",
        markersize = 7.5,
        markerstrokewidth = 0.5,
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

    p2 = plot(
        t_fine,
        v_fine_118,
        # label = "IEEE 39-bus",
        label = "Fine-step continuation",
        title = "(b) IEEE 118-bus";
        common_args...
    )
    scatter!(p2,
        t_adap_118,
        v_adap_118,
        # label="",
        label = "Adaptive HC points",
        markersize = 7.5,
        markerstrokewidth = 0.5,
    )
    plot!(
        p2,
        legend = :bottomleft,
    )
    # display(p1)
    plt = plot(
            p1, p2,
            layout = (1, 2),
            size = (1300, 550),
            )

    display(plt)

    folder = "/Users/aali27/Work/repos/emt_simulator_basic/figures/review1/"
    tosave && savefig(plt, folder*"voltage_continuation.pdf")

end
