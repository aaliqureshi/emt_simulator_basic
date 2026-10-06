# Individual panels for the converter case study, each written as its own PDF so
# they can be assembled into a 2x2 with LaTeX subfigures.
#
#   julia --project=. scripts/plot_converter_panels.jl [output_directory]
#
# Produces, in <out>/panels/:
#   panel_b_current_plane.pdf   converter current plane, limit circle, continuation path
#   panel_b_zoom.pdf            standalone version of the inset (in case the inset is
#                               too small once the 2x2 is assembled)
#   panel_c_divergence.pdf      Newton divergence in BOTH roles, h = 10 us
#   panel_d_continuity.pdf      |i_c| jumps, e_q is continuous, over -100..500 ms
#
# Panel (a) is the circuit schematic, drawn separately.
#
# All panels use the nine-state filtered model so the schematic may legitimately
# show a coupling reactance. The filter block is triangular, so the numbers are
# identical to the seven-state model; that is asserted before anything is drawn.

ENV["GKSwstype"] = "100"

using LinearAlgebra, Printf, ForwardDiff
using Plots, LaTeXStrings
using Barq, MyDiffEq

include(joinpath(@__DIR__, "converter_filtered.jl"))

const OUT = joinpath(isempty(ARGS) ? joinpath(@__DIR__, "..", "figures", "review1") :
                     abspath(ARGS[1]), "panels")

const H_DIVERGE = 1e-4    # step size at which the integrator's Newton diverges
const H_DIVERGE2 = 1e-5     # step size at which the integrator's Newton diverges
const T_PRE     = 0.1       # pre-fault simulation length (s)
const T_POST    = 0.5       # post-event simulation length (s)
const H_SIM     = 1e-4      # step size for the long-horizon runs
const DLAM      = 0.2

# const C_PRE   = "#0072BD"
# const C_WRONG = "#D95319"
# const C_HC    = "#009E73"
# const C_SOL   = "#CC3399"
# const C_PATH  = "#7F7F7F"
# const C_LIM   = "#333333"

const C_PRE   = "#0072B2"   # blue
const C_WRONG = "#D55E00"   # vermillion
const C_HC    = "#009E73"   # bluish green
const C_SOL   = "#CC79A7"   # reddish purple
const C_PATH  = "#7F7F7F"   # gray
const C_LIM   = "#222222"   # near black

style!() = default(fontfamily = "Computer Modern", linewidth = 3.5, markersize = 9,
                   markerstrokewidth = 0, legendfontsize = 15, guidefontsize = 15,
                   tickfontsize = 16, titlefontsize = 15, framestyle = :box,
                   grid = false, dpi = 300,
                   left_margin = 6Plots.mm,
                   right_margin = 2Plots.mm,
                   bottom_margin = 7Plots.mm,
                   top_margin = 2Plots.mm,
                   )

save(pl, name; w, h) = begin
    plot!(pl; size = (w, h), left_margin = 9Plots.mm, bottom_margin = 8Plots.mm,
          right_margin = 6Plots.mm, top_margin = 4Plots.mm)
    savefig(pl, joinpath(OUT, name * ".pdf"))
    savefig(pl, joinpath(OUT, name * ".png"))
end

# ── residual traces ────────────────────────────────────────────────────────────

"""Algebraic residual ‖g‖_inf at a state."""
alg_res(u, p) = (du = zeros(eltype(u), 9); converter_filtered!(du, u, p, 0.0);
                 norm(du[2:end], Inf))

"""
Newton on the implicit-Euler step equations, recording ‖g‖_inf at every iterate.

This reproduces what the integrator does on its first step: it solves
  r(u) = [u_eq - u_eq_prev - h*f_eq(u) ; g(u)] = 0.
Rows 2..9 of r ARE g, so plotting ‖g‖_inf puts this trace and the re-initialization
trace on one axis with a single meaning.
"""
function step_newton_trace(u_prev, p, h; max_iter = 60)
    u = copy(u_prev)
    trace = Float64[]
    for _ in 1:max_iter
        du = zeros(9); converter_filtered!(du, u, p, 0.0)
        r = copy(du); r[EQ_F] = u[EQ_F] - u_prev[EQ_F] - h*du[EQ_F]
        push!(trace, norm(r[2:end], Inf))
        all(isfinite, r) || break
        norm(r, Inf) < 1e-10 && break
        J = ForwardDiff.jacobian(y -> begin
                d = zeros(eltype(y), 9); converter_filtered!(d, y, p, 0.0)
                rr = copy(d); rr[EQ_F] = y[EQ_F] - u_prev[EQ_F] - h*d[EQ_F]; rr
            end, u)
        Δ = try J \ (-r) catch; break end
        all(isfinite, Δ) || break
        u += Δ
        norm(u) > 1e10 && (push!(trace, norm(alg_res(u, p), Inf)); break)
    end
    return trace
end

"""Direct Newton on the post-event algebraic block, recording ‖g‖_inf and the iterates."""
function reinit_newton_trace(u0, p; max_iter = 60)
    u = copy(u0)
    trace = Float64[]; iters = [copy(u)]
    for _ in 1:max_iter
        du = zeros(9); converter_filtered!(du, u, p, 0.0)
        r = du[2:end]
        push!(trace, norm(r, Inf))
        all(isfinite, r) || break
        norm(r, Inf) < 1e-10 && break
        J = ForwardDiff.jacobian(y -> begin
                ut = Vector{eltype(y)}(u); ut[2:9] .= y
                d = similar(ut); converter_filtered!(d, ut, p, 0.0); d[2:9]
            end, u[2:9])
        Δ = try J \ (-r) catch; break end
        all(isfinite, Δ) || break
        u[2:9] .+= Δ
        push!(iters, copy(u))
        norm(u) > 1e10 && break
    end
    return trace, iters
end

# ── panels ─────────────────────────────────────────────────────────────────────

"""Converter current components along the continuation, in true current coordinates."""
function path_currents(hc, p)
    xs = Float64[]; ys = Float64[]
    v_slack, r_line, x_line, i_L, g_f = p[2], p[3], p[4], p[9], p[10]
    z2 = r_line^2 + x_line^2
    for (i, λ) in enumerate(range(0, 1; length = length(hc.vd_hist)))
        v_d, v_q = hc.vd_hist[i], hc.vq_hist[i]
        i_ld = (r_line*(v_d - v_slack) + x_line*v_q)/z2
        i_lq = (-x_line*(v_d - v_slack) + r_line*v_q)/z2
        i_Ld, i_Lq = _load_currents(v_d, v_q, i_L)
        push!(xs, i_ld + i_Ld + λ*g_f*v_d)
        push!(ys, i_lq + i_Lq + λ*g_f*v_q)
    end
    return xs, ys
end

function panel_current_plane(p, u_pre, u_hc, hc, nr_iters; inset = false)
    # legendfontsize=13
    i_max = p[11]; demand = p[9] + p[10]
    θ = range(0, 2π; length = 800)
    pl = plot(xlims = (0.2, 1.6),
             xticks = collect(0.2:0.4:1.4),
            # xticks=[0.2, 0.],
            ylims = (-0.52, 0.75),
              xlabel = L"i_{cd}\ \mathrm{(p.u.)}", ylabel = L"i_{cq}\ \mathrm{(p.u.)}",
            #   legend = :bottomleft,
            legend=(0.32, 0.37),
            #    legendfontsize = legendfontsize, 
            #   aspect_ratio = :equal,
            # title="(b)",
              )
    plot!(pl, i_max .* cos.(θ), i_max .* sin.(θ); 
        color = C_LIM, 
        # linewidth = 3.5,
        label = L"i_{cd}^2 + i_{cq}^2 = i_{max}^2",
        # label="",
        # label = L"i_{max}^2",
        )

    xs, ys = path_currents(hc, p)
    plot!(pl, xs, ys; color = C_PATH,
    #  linewidth = 3.5,
      linestyle = :dash,
          label = "Continuation path",
        #   label="",
          )
    xs_mark = [xs[2], xs[5]]
    ys_mark = [ys[2], ys[5]]
    scatter!(pl, xs_mark, ys_mark; color = C_HC, markersize = 10, label = "")
    offs = [(0.0,0.0), (-0.12,0.12), (-0.08,0.15), (0.135,0.03), (-0.15,-0.1), (0.0,0.0)]
    for (i, λ) in enumerate(range(0, 1; length = length(xs)))
        (i == 1 || i == length(xs) ) && continue
    end
    annotate!(pl, xs[2] + offs[2][1], ys[2] + offs[2][2],
              text(@sprintf("λ=0.2"), 13, C_HC))
    # annotate!(pl, xs[3] + offs[3][1], ys[3] + offs[3][2],
    #           text(@sprintf("λ = 0.4"), 9, C_HC))
    # annotate!(pl, xs[4] + offs[4][1], ys[4] + offs[4][2],
    #           text(@sprintf("λ = 0.6"), 9, C_HC))
    annotate!(pl, xs[5] + offs[5][1], ys[5] + offs[5][2],
              text(@sprintf("λ=0.8"), 13, C_HC))

    # quote the excursion after a fixed, small number of updates: the maximum over a
    # divergent sequence depends only on where the iteration is cut off.
    nshow = min(4, length(nr_iters))
    far = maximum(hypot(nr_iters[k][ICD_F], nr_iters[k][ICQ_F]) for k in 1:nshow)


    scatter!(pl, [u_pre[ICD_F]], [u_pre[ICQ_F]]; marker = :utriangle, markersize = 12,
             color = C_PRE,
            # label = L"\mathrm{Pre\!-\!event}\ |i_c| = 0.32",
            # label="Pre-event solution |i_c| = 0.32",
            # label="",
            label="Pre-event solution",
            )
    # scatter!(pl, [demand], [0.0]; marker = :xcross, markersize = 10, markerstrokewidth = 3,
            #  color = C_WRONG, label = L"\mathrm{GFM\ demand},\ |i_c| = 1.439")
    scatter!(pl, [u_hc[ICD_F]], [u_hc[ICQ_F]]; marker = :star5, markersize = 14,
             color = C_SOL, 
            #  label = L"\mathrm{After\ re\!-\!init}\ |i_c| = 1.430",
            #  label="",
             label="Post-event solution",
             )

    if inset
        plot!(pl; inset = (1, bbox(0.50, 0.04, 0.44, 0.32, :bottom, :left)), subplot = 2)
        sp = pl[2]
        θz = range(-0.10, 0.10; length = 300)
        plot!(sp, i_max .* cos.(θz), i_max .* sin.(θz); color = C_LIM, linewidth = 2.0,
              label = "", xlims = (1.4235, 1.4445), ylims = (-0.070, 0.035),
              xticks = [1.43, 1.44], yticks = [-0.05, 0.0],
              tickfontsize = 8, framestyle = :box, grid = false)
        scatter!(sp, [demand], [0.0]; marker = :xcross, markersize = 8,
                 markerstrokewidth = 2.5, color = C_WRONG, label = "")
        scatter!(sp, [u_hc[ICD_F]], [u_hc[ICQ_F]]; marker = :star5, markersize = 10,
                 color = C_SOL, label = "")
        annotate!(sp, 1.4340, 0.023, text("0.63% outside", 7, C_WRONG),)
    end
    return pl
end

function panel_zoom(p, u_hc)
    i_max = p[11]; demand = p[9] + p[10]
    θ = range(-0.16, 0.16; length = 400)
    pl = plot(xlabel = L"i_{cd}\ \mathrm{(p.u.)}", ylabel = L"i_{cq}\ \mathrm{(p.u.)}",
              legend = :bottomleft, legendfontsize = 9,
              xlims = (1.4175, 1.4475), ylims = (-0.088, 0.062))
    plot!(pl, i_max .* cos.(θ), i_max .* sin.(θ); color = C_LIM, linewidth = 2.4,
          label = L"|i_c| = i_{max} = 1.43")
    scatter!(pl, [demand], [0.0]; marker = :xcross, markersize = 11, markerstrokewidth = 3,
             color = C_WRONG, label = L"\mathrm{GFM\ demand} = 1.439")
    scatter!(pl, [u_hc[ICD_F]], [u_hc[ICQ_F]]; marker = :star5, markersize = 13,
             color = C_SOL, label = "Re-initialized point")
    annotate!(pl, 1.4388, 0.022, text("outside the circle\nby 0.63%", 9, C_WRONG))
    return pl
end

function panel_divergence(step_trace, step_trace2, reinit_trace; nshow = 10)
    a = step_trace[1:min(nshow, end)]
    a2 = step_trace2[1:min(nshow, end)]
    b = reinit_trace[1:min(nshow, end)]
    # legendfontsize = 13
    pl = plot(xlabel = "Newton iteration", 
             ylabel = L"\|G(\mathbf{z})\|_\infty",
              yscale = :log10, legend = :bottomright,
            #    legendfontsize = 9,
            #   ylims = (5e-1, 1e4), 
            #   yticks = 10.0 .^ (0:6.5),
                yticks = [10^1, 10^3, 10^6],
              xlims = (0.7, 10.3),
              xticks=[1, 5, 10],
              )
    plot!(pl, 1:length(a), a; color = C_WRONG, linestyle = :dash,
          marker = :diamond, markersize = 6,
          yscale=:log10,
        #   label = L"\mathrm{DAE\ step,\ no\ re\!-\!init}\ (h = 10\ \mu s)",
        #   label="DAE step, no re-init",
          label="DAE step, (h=100 μs)",
        #   legendfontsize = legendfontsize,
            markerstrokecolor = :black,
              markerstrokewidth = 1.0,
          )
    plot!(pl, 1:length(a2), a2; color = C_SOL, linestyle = :dashdot,
          marker = :utriangle, markersize = 6,
          yscale=:log10,
        #   label = L"\mathrm{DAE\ step,\ no\ re\!-\!init}\ (h = 10\ \mu s)",
        #   label="DAE step, no re-init",
        label="DAE step, (h=10 μs)",
        #   legendfontsize = legendfontsize,
            markerstrokecolor = :black,
              markerstrokewidth = 1.0,
          )
    plot!(pl, 1:length(b), b; color = C_PRE, linestyle = :solid,
          marker = :circle, markersize = 6,
          yscale=:log10,
        #   label = L"\mathrm{Direct\ Newton\ re\!-\!init}",
          label="Direct NR re-init",
        #   legendfontsize = legendfontsize,
        # title="(a)"
        markerstrokecolor = :black,
        markerstrokewidth = 1.0,
          )
    return pl
end

function panel_continuity(t_pre, u_pre_series, t_post, u_post_series, p)
    ic_pre  = [current_f(u) for u in u_pre_series]
    ic_post = [current_f(u) for u in u_post_series]
    eq_pre  = [u[EQ_F] for u in u_pre_series]
    eq_post = [u[EQ_F] for u in u_post_series]
    # legendfontsize = 13

    pl = plot(xlabel = "Time (ms)", 
            ylabel = L"|i_c|\ \mathrm{(p.u.)}",
              legend = :right,
            #    legendfontsize = legendfontsize,
               yguidefontcolor = C_PRE,
            #   xlims = (-1e3*T_PRE, 1e3*T_POST),
                ylims = (0.15, 1.62),
                # xlims = (-1.0,5.0),
                # xticks=collect(0:1:6), 
              framestyle = :box,
            #   title="(c)",
              )
    t_post_added = @. t_pre[end] + t_post
    t_full = vcat(t_pre, t_post_added)
    # @show t_pre
    t_pre_ms = 1e3 .* (t_pre .- t_pre[end])
    t_post_ms = 1e3 .* t_post
    # t_full = 1e3 .* t_full 
    t_full = vcat(t_pre_ms, t_post_ms)
    # idx = findfirst(t-> t == 400, t_full)
    ic = vcat(ic_pre, ic_post)
    eq = vcat(eq_pre, eq_post)
    idx_end = 2000
    t_full = t_full[1:end-idx_end]
    ic = ic[1:end-idx_end]
    eq = eq[1:end-idx_end]
    hline!(pl, [p[11]]; color = C_LIM, linestyle = :dashdot, 
    # linewidth = 1.8,
     label = L"i_{max}")
    # plot!(pl, 1e3 .* (t_pre .- T_PRE), ic_pre; color = C_PRE, label = L"|i_c|\ \mathrm{(left\ axis)}")
    plot!(pl, t_full, ic; color = C_PRE, label = L"|i_c|\ \mathrm{(left\ axis)}")


    # single legend: a dummy series on the main axis stands in for the e_q curve
    plot!(pl, [NaN], [NaN]; color = C_HC, 
    # linewidth = 3.0,
          label = L"x\ \mathrm{(right\ axis)}")
    
    vline!(pl, [0.0];
          color = :black,
          linestyle = :dot,
          linewidth = 2.0,
          label = "",
      )

    ax2 = twinx(pl)
    plot!(ax2, t_full, eq; color = C_HC, 
    # linewidth = 3.5,
    ylabel = L"x\ \mathrm{(p.u.)}", legend = false, yguidefontcolor = C_HC,
    # ylims = (-0.0085, 0.0022),
    label = "",
    guidefontsize = 14, 
    tickfontsize = 12,
    ylims=(-0.027, 0.003),
    )
# plot!(ax2, 1e3 .* t_post, eq_post; color = C_HC, linewidth = 3.5, label = "")
    return pl
end

# ── driver ─────────────────────────────────────────────────────────────────────

function main_panels()
    mkpath(OUT); gr(); style!()
    p = P_F

    chk = assert_matches_seven_state(p)
    @printf("nine-state agrees with seven-state to %.2e (%d Newton iterations)\n\n",
            chk.agree, chk.iters)

    u_pre = pre_event_state_filtered(p)
    u_hc  = copy(u_pre)
    hc = solve_homotopy!(u_hc, p, ADDRESS_F; tol=1e-9, max_iter=100, Δλ=DLAM,
                         vd_idx=VD_F, vq_idx=VQ_F, always_new=true,
                         model! = converter_filtered!)
    @assert hc.converged "continuation failed"
    @printf("re-init: V_pcc = %.6f  |i_c| = %.6f  |w| = %.6f  (%d Newton iterations)\n",
            voltage_f(u_hc), current_f(u_hc), terminal_f(u_hc), hc.total_iters)

    step_trace = step_newton_trace(u_pre, post_event_f(p), H_DIVERGE)
    step_trace2 = step_newton_trace(u_pre, post_event_f(p), H_DIVERGE2)
    reinit_trace, nr_iters = reinit_newton_trace(u_pre, post_event_f(p))
    @printf("step Newton (h=10us) : %d iterations, final residual %.3e\n",
            length(step_trace), last(step_trace))
    @printf("direct re-init Newton: %d iterations, final residual %.3e\n",
            length(reinit_trace), last(reinit_trace))
    @assert last(step_trace)   > 1.0 "step Newton did not diverge"
    @assert last(reinit_trace) > 1.0 "re-init Newton did not diverge"

    mass = mass_matrix_f()
    sol_pre = MyDiffEq.Solve(
        MyDiffEq.ODEProblem(converter_filtered!, copy(u_pre), (0.0, T_PRE),
                            pre_event_f(p), mass), H_SIM;
        method = :Euler, adaptive = false, always_new = true)
    sol_post = MyDiffEq.Solve(
        MyDiffEq.ODEProblem(converter_filtered!, copy(u_hc), (0.0, T_POST),
                            post_event_f(p), mass), H_SIM;
        method = :Euler, adaptive = false, always_new = true)
    @assert sol_pre.retcode == :Success && sol_post.retcode == :Success
    @printf("pre-fault run  : %.0f ms, |i_c| = %.4f, e_q = %+.5f\n",
            1e3*T_PRE, current_f(sol_pre.u[end]), sol_pre.u[end][EQ_F])
    @printf("post-event run : %.0f ms, |i_c| = %.4f, e_q = %+.5f, V = %.4f\n",
            1e3*T_POST, current_f(sol_post.u[end]), sol_post.u[end][EQ_F],
            voltage_f(sol_post.u[end]))
    @printf("e_q continuity across the event: |Δe_q| = %.2e\n",
            abs(sol_post.u[1][EQ_F] - sol_pre.u[end][EQ_F]))

    pb = panel_current_plane(p, u_pre, u_hc, hc, nr_iters)
    pz = panel_zoom(p, u_hc)
    pc = panel_divergence(step_trace, step_trace2, reinit_trace)
    pd = panel_continuity(sol_pre.time, sol_pre.u, sol_post.time, sol_post.u, p)

    save(pb, "panel_b_current_plane"; w = 760, h = 660)
    save(pz, "panel_b_zoom";          w = 700, h = 520)
    save(pc, "panel_c_divergence";    w = 760, h = 560)
    save(pd, "panel_d_continuity";    w = 900, h = 560)

    pz = plot(pc, pb, pd,
         layout=(1,3),
         size=(1800,380),
         tickfontsize = 15,  # Bump up fonts so they are legible
        guidefontsize = 15,
        legendfontsize=13,
        # legend = :none,     # You likely won't have room for 3 legends
        bottom_margin = 11Plots.mm,
        top_margin = 3Plots.mm,
        left_margin = 10Plots.mm,
        right_margin=9Plots.mm,
    )

    folder = "/Users/aali27/Work/repos/emt_simulator_basic/figures/review1/"
    savefig(pz, folder*"3panels.pdf")

    println("\nWrote panels to $(normpath(OUT))")
end

main_panels()
