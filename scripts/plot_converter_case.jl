# Figures for the GFM converter re-initialization case study.
#
#   julia --project=. scripts/plot_converter_case.jl [output_directory]
#
# Everything plotted here is an actual state variable of the 7-state DAE. Nothing
# is eliminated, reduced or projected onto a derived coordinate:
#
#   (a) converter current plane (i_cd, i_cq). The current limit
#       i_cd^2 + i_cq^2 = i_max^2 is one of the six algebraic equations and
#       involves only these two variables, so it is drawn exactly as a circle.
#   (b) algebraic residual against Newton iteration.
#   (c) PCC voltage against time.
#   (d) converter current magnitude against time.
#
# The physical story the figure tells:
#   pre-event the converter is voltage-forming and draws |i_c| = i_L = 0.32 p.u.
#   The fault would require |i_c| = i_L + g_f = 1.439 p.u. > i_max = 1.43 p.u.,
#   so the current limit binds and the model must switch to the limited mode.
#   The inherited state is then inconsistent with the new algebraic equations.

ENV["GKSwstype"] = "100"

using LinearAlgebra, Printf, ForwardDiff
using Plots, LaTeXStrings
using Barq, MyDiffEq

include(joinpath(@__DIR__, "converter_convex_homotopy.jl"))

const OUT   = isempty(ARGS) ? joinpath(@__DIR__, "..", "figures", "review1") : abspath(ARGS[1])
const TEND  = 0.01
const DLAM  = 0.1
const H_OK  = 1e-5      # step size used for the re-initialized run
const H_BAD = 1e-4      # step size at which the un-re-initialized run converges to the wrong branch

const C_PRE  = "#0072BD"
const C_WRONG= "#D95319"
const C_HC   = "#009E73"
const C_SOL  = "#CC3399"
const C_PATH = "#7F7F7F"
const C_LIM  = "#333333"

style!() = default(fontfamily = "Computer Modern", linewidth = 3.0, markersize = 6,
                   markerstrokewidth = 0, legendfontsize = 10, guidefontsize = 13,
                   tickfontsize = 11, titlefontsize = 13, framestyle = :box,
                   grid = false, dpi = 300)

"""Converter current components implied by (v_d, v_q) at homotopy level λ.

The line equations and the PCC current balance are linear in
(i_ld, i_lq, i_cd, i_cq) once the voltage is fixed, so they invert in closed
form. This is used only to turn the stored continuation path into current
coordinates; it agrees with the solver state to machine precision.
"""
function converter_current_components(v_d, v_q, p, λ)
    v_slack, r_line, x_line, i_L, g_f = p[2], p[3], p[4], p[9], p[10]
    z2 = r_line^2 + x_line^2
    i_ld = (r_line*(v_d - v_slack) + x_line*v_q)/z2
    i_lq = (-x_line*(v_d - v_slack) + r_line*v_q)/z2
    i_Ld, i_Lq = _load_currents(v_d, v_q, i_L)
    return i_ld + i_Ld + λ*g_f*v_d, i_lq + i_Lq + λ*g_f*v_q
end

"""Direct Newton on the post-event algebraic block, recording every iterate.

Same iteration as `solve_newton!`; reproduced here only so the iterates can be
plotted.
"""
function newton_iterates(u_start, p; max_iter = 40, tol = 1e-9)
    u = copy(u_start)
    hist = [copy(u)]
    for _ in 1:max_iter
        g = zeros(7); converter_homotopy!(g, u, post_event(p), 0.0)
        r = g[2:end]
        all(isfinite, r) || break
        norm(r, Inf) <= tol && break
        J = ForwardDiff.jacobian(y -> begin
                ut = Vector{eltype(y)}(u); ut[2:7] .= y
                dut = similar(ut); converter_homotopy!(dut, ut, post_event(p), 0.0)
                dut[2:7]
            end, u[2:7])
        Δ = try J \ (-r) catch; break end
        all(isfinite, Δ) || break
        u[2:7] .+= Δ
        push!(hist, copy(u))
        norm(u) > 1e8 && break
    end
    return hist
end

simulate(u_init, p, h) =
    MyDiffEq.Solve(MyDiffEq.ODEProblem(converter_homotopy!, copy(u_init), (0.0, TEND),
                                       post_event(p), begin
                                           m = zeros(7,7); m[EQ,EQ] = 1.0; m
                                       end),
                   h; method = :Euler, adaptive = false, always_new = true)

# ── Panels ─────────────────────────────────────────────────────────────────────

"""(a) Converter current plane: the limit circle, the demand that violates it, and
the continuation path that walks onto it."""
function panel_current_plane(p, u_pre, u_hc, hc, nr_hist)
    i_max = p[11]
    θ = range(0, 2π; length = 800)
    pl = plot(title = "(a) Converter current plane", xlims = (0.05, 1.80), ylims = (-0.4, 0.82),
              xlabel = L"i_{cd}\ \mathrm{(p.u.)}", ylabel = L"i_{cq}\ \mathrm{(p.u.)}",
              legend = :topleft, legendfontsize = 9)
    plot!(pl, i_max .* cos.(θ), i_max .* sin.(θ); color = C_LIM, linewidth = 2.2,
          label = L"i_{cd}^2 + i_{cq}^2 = i_{max}^2")

    xs = Float64[]; ys = Float64[]
    for (i, λ) in enumerate(range(0, 1; length = length(hc.vd_hist)))
        a, b = converter_current_components(hc.vd_hist[i], hc.vq_hist[i], p, λ)
        push!(xs, a); push!(ys, b)
    end
    plot!(pl, xs, ys; color = C_PATH, linewidth = 2.2, linestyle = :dash, label = "Continuation path")
    scatter!(pl, xs, ys; color = C_HC, markersize = 6, label = "")
    offs = [(0.0,0.0), (-0.055,-0.105), (-0.125,0.105), (0.135,0.105), (0.165,0.010), (0.0,0.0)]
    # for (i, λ) in enumerate(range(0, 1; length = length(xs)))
    #     (i == 1 || i == length(xs)) && continue
    #     annotate!(pl, xs[i] + offs[i][1], ys[i] + offs[i][2],
    #               text(@sprintf("\u03bb = %.1f", λ), 9, C_HC))
    # end

    # direct Newton leaves the admissible region on its first update
    # far = maximum(hypot(u[ICD], u[ICQ]) for u in nr_hist)
    # plot!(pl, [u_pre[ICD], u_pre[ICD]], [u_pre[ICQ], -0.34]; color = C_WRONG, linewidth = 2.6,
    #       arrow = :closed, label = "Direct Newton iterates")
    # annotate!(pl, 0.40, -0.30,
    #           text(@sprintf("escapes to %.0f p.u.", far), 9, C_WRONG, :left))

    scatter!(pl, [u_pre[ICD]], [u_pre[ICQ]]; marker = :utriangle, markersize = 11,
             color = C_PRE, label = L"\mathrm{Pre\!-\!event\ (GFM)},\ |i_c| = 0.32")
    # scatter!(pl, [p[9] + p[10]], [0.0]; marker = :xcross, markersize = 10,
    #          markerstrokewidth = 3, color = C_WRONG, label = L"\mathrm{GFM\ demand},\ |i_c| = 1.439")
    scatter!(pl, [u_hc[ICD]], [u_hc[ICQ]]; marker = :star5, markersize = 13,
             color = C_SOL, label = L"\mathrm{After\ re\!-\!init},\ |i_c| = 1.430")
    return pl
end

"""Inset: the GFM demand sits 0.63% outside the circle — too small to see unzoomed."""
function panel_zoom(p, u_hc)
    i_max = p[11]; demand = p[9] + p[10]
    θ = range(-0.16, 0.16; length = 400)
    pl = plot(title = "(b) Detail at the limit", xlabel = L"i_{cd}\ \mathrm{(p.u.)}",
              ylabel = L"i_{cq}\ \mathrm{(p.u.)}", legend = :bottomleft, legendfontsize = 9,
              xlims = (1.4175, 1.4475), ylims = (-0.088, 0.062))
    plot!(pl, i_max .* cos.(θ), i_max .* sin.(θ); color = C_LIM, linewidth = 2.4,
          label = L"|i_c| = i_{max} = 1.43")
    scatter!(pl, [demand], [0.0]; marker = :xcross, markersize = 11, markerstrokewidth = 3,
             color = C_WRONG, label = "GFM demand = 1.439")
    scatter!(pl, [u_hc[ICD]], [u_hc[ICQ]]; marker = :star5, markersize = 13, color = C_SOL,
             label = "Re-initialized point")
    annotate!(pl, 1.4388, 0.022, text("outside the circle\nby 0.63%", 9, C_WRONG))
    return pl
end

"""(c) Residual histories."""
function panel_residual(nr, hc)
    bounds = Float64[]; acc = 0
    for k in 1:length(hc.stage_iters)-1
        acc += hc.stage_iters[k]; push!(bounds, acc + 0.5)
    end
    pl = plot(title = "(c) Re-initialization", xlabel = "Newton iteration",
              ylabel = L"\|G(\mathbf{z})\|_\infty", yscale = :log10,
              legend = :bottomleft, 
              ylims = (1e-1, 1e4),
              )
    # for b in bounds
    #     vline!(pl, [b]; color = :gray75, linewidth = 1.0, linestyle = :dot, label = "")
    # end
    # plot!(pl, 1:length(nr.residuals), max.(nr.residuals, 1e-10); color = C_WRONG,
    #       linestyle = :dash, marker = :diamond, markersize = 4, label = "Direct DAE Integration")
    plot!(pl, 1:length(nr.residuals[1:10]), nr.residuals[1:10]; color = C_WRONG, yscale=:log10,
          linestyle = :dash, marker = :diamond, markersize = 9, label = "Direct DAE Integration")
    # plot!(pl, 1:length(hc.residuals), max.(hc.residuals, 1e-10); color = C_HC,
    #       marker = :square, markersize = 4, label = "Homotopy continuation")
    # plot!(pl, 1:length(nr.residuals[1:30]), nr.correction_norm[1:30]; color = C_WRONG,
    # linestyle = :dash, marker = :diamond, markersize = 4, label = "Direct DAE Integration")
    return pl
end

"""Time-domain response: correct restart vs the branch reached without re-initialization."""
function panel_time(p, u_pre, sol_hc, sol_bad, quantity, title)
    f, ylab = quantity === :V  ? (voltage, L"V\ \mathrm{(p.u.)}") :
              quantity === :ic ? (current, L"|i_c|\ \mathrm{(p.u.)}") :
                                 (u -> u[ICQ], L"i_{cq}\ \mathrm{(p.u.)}")
    pl = plot(title = title, xlabel = "Time (ms)", ylabel = ylab,
              legend = quantity === :ic ? :right : :best, legendfontsize = 9)
    plot!(pl, [-2.0, 0.0], [f(u_pre), f(u_pre)]; color = C_PRE, label = "Pre-event (GFM)")
    if quantity === :ic
        hline!(pl, [p[11]]; color = C_LIM, linestyle = :dashdot, linewidth = 2.0, label = L"i_{max}")
        scatter!(pl, [0.0], [p[9]+p[10]]; marker = :xcross, markersize = 9,
                 markerstrokewidth = 3, color = C_WRONG, label = "GFM demand (infeasible)")
    end
    plot!(pl, 1e3 .* sol_bad.time, [f(u) for u in sol_bad.u]; color = C_WRONG,
          linestyle = :dash, label = "No re-init, h = 100 \u03bcs")
    plot!(pl, 1e3 .* sol_hc.time, [f(u) for u in sol_hc.u]; color = C_HC,
          label = "Re-initialized, h = 10 \u03bcs")
    vline!(pl, [0.0]; color = :gray60, linewidth = 1.2, linestyle = :dot, label = "")
    return pl
end

# ── Driver ─────────────────────────────────────────────────────────────────────

function main_plot()
    mkpath(OUT); gr(); style!()
    p = P
    @printf("parameters: i_L = %.2f, g_f = %.3f, i_max = %.2f\n", p[9], p[10], p[11])
    @printf("GFM demand under fault = i_L + g_f = %.4f  >  i_max = %.2f  (violated by %.2f%%)\n",
            p[9]+p[10], p[11], 100*((p[9]+p[10])/p[11] - 1))

    u_pre = pre_event_state(p)
    u_nr = copy(u_pre)
    nr = solve_newton!(u_nr, post_event(p), ADDRESS; tol = 1e-9, max_iter = 100,
                       always_new = true, model! = converter_homotopy!)
    u_hc = copy(u_pre)
    hc = solve_homotopy!(u_hc, p, ADDRESS; tol = 1e-9, max_iter = 100, Δλ = DLAM,
                         vd_idx = VD, vq_idx = VQ, always_new = true,
                         model! = converter_homotopy!)
    @assert hc.converged "continuation failed"
    @assert !nr.converged "direct Newton unexpectedly converged"
    nr_hist = newton_iterates(u_pre, p)

    sol_hc  = simulate(u_hc,  p, H_OK)
    sol_bad = simulate(u_pre, p, H_BAD)
    @assert sol_hc.retcode == :Success
    @printf("re-initialized run   : %s, V(10 ms) = %.6f\n", sol_hc.retcode, voltage(sol_hc.u[end]))
    if sol_bad.retcode == :Success
        @printf("no re-init, h=100 us : Success but V(10 ms) = %.6f  <-- different trajectory\n",
                voltage(sol_bad.u[end]))
    else
        @printf("no re-init, h=100 us : %s\n", sol_bad.retcode)
    end
    for h in (1e-5, 1e-6)
        s = simulate(u_pre, p, h)
        @printf("no re-init, h=%3d us : %s\n", round(Int, h*1e6), s.retcode)
    end

    pa = panel_current_plane(p, u_pre, u_hc, hc, nr_hist)
    pb = panel_zoom(p, u_hc)
    pc = panel_residual(nr, hc)
    pd = panel_time(p, u_pre, sol_hc, sol_bad, :V,   "(d) PCC voltage")
    pe = panel_time(p, u_pre, sol_hc, sol_bad, :ic,  "(e) Converter current magnitude")
    pf = panel_time(p, u_pre, sol_hc, sol_bad, :icq, "(f) Converter q-axis current")
    annotate!(pd, 5.2, 1.155, text("no re-init at h = 10 \u03bcs or 1 \u03bcs:\nintegrator cannot start", 9, C_WRONG))

    fig = plot(pa, pb, pc, pd, pe, pf; layout = (2, 3),
               size = (1850, 1000), left_margin = 10Plots.mm, bottom_margin = 10Plots.mm,
               top_margin = 5Plots.mm, right_margin = 5Plots.mm)
    savefig(fig, joinpath(OUT, "converter_case.pdf"))
    savefig(fig, joinpath(OUT, "converter_case.png"))

    for (nm, pl, w, h) in (("converter_current_plane", pa, 760, 660),
                           ("converter_limit_detail",  pb, 700, 520),
                           ("converter_residual",      pc, 700, 520),
                           ("converter_voltage",       pd, 780, 520),
                           ("converter_current_time",  pe, 780, 520),
                           ("converter_icq_time",      pf, 780, 520))
        plot!(pl; size = (w, h), left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)
        savefig(pl, joinpath(OUT, nm * ".pdf"))
        savefig(pl, joinpath(OUT, nm * ".png"))
    end
    println("\nWrote figures to $(normpath(OUT))")
end

main_plot()
