# Manifold figures for the GFM converter re-initialization case study.
#
#   julia --project=. scripts/plot_converter_manifold.jl [output_directory]
#
# Reproduces the presentation used for the static case study, but for the DAE
# case. The post-event algebraic system has six unknowns, so it cannot be drawn
# directly against a single variable. It is reduced exactly to one scalar
# equation as follows.
#
# Of the six algebraic rows of H(y,λ), four are linear in (i_ld, i_lq, i_cd, i_cq)
# for a given voltage and are eliminated in closed form. Of the two remaining
# mode-defining rows, row 6 is strictly monotone in v_d over the region of
# interest (d row6/d v_d runs from 1.0 at λ=0 to 2.0 at λ=1), so it defines a
# unique v_d = v_d(v_q; λ). Substituting into row 7 leaves
#
#     R(v_q; λ) = row7( v_d(v_q; λ), v_q; λ ),
#
# a scalar function whose zeros are exactly the solutions of the full
# six-variable system. No approximation is made: the reduced path reproduces the
# six-variable solver's continuation path to all printed digits.
#
# At λ=0, R(v_q;0) = v_q, so the pre-event solution sits at v_q = 0. At λ=1 the
# curve has two zeros — the physical branch and a spurious high-voltage branch
# that direct Newton can be attracted to.

ENV["GKSwstype"] = "100"

using LinearAlgebra, Printf
using Plots, LaTeXStrings
using Barq

include(joinpath(@__DIR__, "converter_convex_homotopy.jl"))

const OUT = isempty(ARGS) ? joinpath(@__DIR__, "..", "figures", "review1") : abspath(ARGS[1])
const LAMBDAS = collect(range(0.0, 1.0; length=6))

# Palette chosen to match the static case-study figure.
const C_PRE   = "#0072BD"   # pre-event manifold
const C_POST  = "#E8762C"   # post-event manifold
const C_MID   = "#17A398"   # intermediate manifolds
const C_PATH  = "#7F7F7F"   # continuation path
const C_ITER  = "#EF8A7A"   # continuation iterates
# const C_ITER = "DA291C"
const C_SOL   = "#E0219A"   # post-event solution
const C_INIT  = "#17A398"   # pre-event solution

# ── Exact scalar reduction ─────────────────────────────────────────────────────

"""Mode-defining rows 6 and 7 of H(y,λ) with the four linear rows eliminated."""
function mode_rows(v_d, v_q, λ, p)
    v_ref, v_slack, r_line, x_line = p[1], p[2], p[3], p[4]
    kp, i_L, g_f, i_max = p[8], p[9], p[10], p[11]
    e_q = 0.0                                     # frozen differential state
    z2 = r_line^2 + x_line^2
    V = hypot(v_d, v_q)
    i_ld = (r_line*(v_d - v_slack) + x_line*v_q)/z2
    i_lq = (-x_line*(v_d - v_slack) + r_line*v_q)/z2
    i_Ld, i_Lq = _load_currents(v_d, v_q, i_L)
    i_cd = i_ld + i_Ld + λ*g_f*v_d
    i_cq = i_lq + i_Lq + λ*g_f*v_q
    row6 = (1-λ)*(v_d - v_ref) + λ*(i_cq - e_q - kp*(v_ref - V))
    row7 = (1-λ)*v_q           + λ*(i_cd^2 + i_cq^2 - i_max^2)
    return row6, row7
end

"""Unique v_d solving row6(v_d; v_q, λ) = 0, by scalar Newton with warm start."""
function v_d_of(v_q, λ, p; v_d0=1.0, tol=1e-13, max_iter=80)
    v_d = v_d0
    for _ in 1:max_iter
        f = mode_rows(v_d, v_q, λ, p)[1]
        h = 1e-7
        df = (mode_rows(v_d+h, v_q, λ, p)[1] - mode_rows(v_d-h, v_q, λ, p)[1])/(2h)
        abs(df) < 1e-12 && return NaN
        Δ = clamp(-f/df, -0.5, 0.5)
        v_d += Δ
        abs(Δ) < tol && return v_d
    end
    return NaN
end

"""Reduced scalar residual R(v_q; λ)."""
function R_reduced(v_q, λ, p; v_d0=1.0)
    v_d = v_d_of(v_q, λ, p; v_d0=v_d0)
    return isnan(v_d) ? NaN : mode_rows(v_d, v_q, λ, p)[2], v_d
end

"""R(v_q; λ) sampled over `grid`, sweeping with a warm start to stay on one branch."""
function reduced_curve(grid, λ, p; v_d0=1.0)
    vals = similar(grid)
    v_d = v_d0
    for (i, v_q) in enumerate(grid)
        r, v_d_new = R_reduced(v_q, λ, p; v_d0=v_d)
        vals[i] = r
        isnan(v_d_new) || (v_d = v_d_new)
    end
    return vals
end

"""Zeros of R(·;λ) found by scanning `grid` for sign changes and refining."""
function reduced_roots(grid, λ, p)
    vals = reduced_curve(grid, λ, p)
    roots = Float64[]
    for i in 1:length(grid)-1
        (isnan(vals[i]) || isnan(vals[i+1])) && continue
        vals[i] == 0 && push!(roots, grid[i])
        if vals[i]*vals[i+1] < 0
            a, b = grid[i], grid[i+1]
            for _ in 1:80                         # bisection: robust near the fold
                m = (a+b)/2
                fm = R_reduced(m, λ, p; v_d0=1.0)[1]
                (R_reduced(a, λ, p; v_d0=1.0)[1])*fm <= 0 ? (b = m) : (a = m)
            end
            push!(roots, (a+b)/2)
        end
    end
    return unique(r -> round(r; digits=9), roots)
end

# ── Figures ────────────────────────────────────────────────────────────────────

function style!()
    default(fontfamily = "Computer Modern", linewidth = 3.5, markersize = 7,
            markerstrokewidth = 0, legendfontsize = 13, guidefontsize = 15,
            tickfontsize = 15, titlefontsize = 15, framestyle = :box,
            grid = false, dpi = 300)
end

"""
Panel (a): the event moves the manifold, so the inherited point is no longer a root.

Zoomed to the root region. The post-event curve is shallow between its two zeros
(|R| <= 0.09) while rising to O(1) just outside, so the wider window used in
panel (b) would hide both crossings.
"""
function panel_event(grid, p, v_init, v_sol, v_spurious)
    # xl, yl = (-0.26, 0.10), (-0.24, 0.26)
    # xl, yl = (-0.8, 0.35), (-1.5, 1.4)
    xl, yl = (-0.1, 0.05), (-0.5, 0.5)
    pl = plot(title = "(a)", xlabel = L"v_q\ \mathrm{(p.u.)}", ylabel = L"H(v_q,\lambda)",
            #   legend = :topleft,
            legend = (0.2, 0.92),
               xlims = xl, ylims = yl)
    plot!(pl, collect(xl), [0.0, 0.0]; color = :black, linewidth = 1.6,
          label = L"R(v_q;\lambda) = 0")
    plot!(pl, grid, reduced_curve(grid, 0.0, p); color = C_PRE,  label = "Pre-event manifold")
    plot!(pl, grid, reduced_curve(grid, 1.0, p); color = C_POST, label = "Post-event manifold")
    if !isnan(v_spurious)
        scatter!(pl, [v_spurious], [0.0]; marker = :circle, markersize = 10,
                 color = :white, markerstrokewidth = 3.0, markerstrokecolor = C_POST,
                 label = "Spurious solution")
        # annotate!(pl, v_spurious, 0.055, text(L"v_{\mathrm{spur}}", 15, C_POST))
    end
    scatter!(pl, [v_init], [0.0]; marker = :utriangle, markersize = 13, color = C_INIT,
             label = "Pre-event solution")
    scatter!(pl, [v_sol], [0.0]; marker = :star5, markersize = 15, color = C_SOL,
             label = "Post-event solution")
    # The two roots are only 0.011 p.u. apart, so both labels are lifted above the
    # axis and separated horizontally with short leader lines.
    # plot!(pl, [v_init+0.002, v_init+0.028], [0.012, 0.090]; color = C_INIT, linewidth = 1.2, label = "")
    # annotate!(pl, v_init+0.036, 0.108, text(L"v_{\mathrm{init}}", 15, C_INIT))
    # plot!(pl, [v_sol-0.003, v_sol-0.062], [0.012, 0.112]; color = C_SOL, linewidth = 1.2, label = "")
    # annotate!(pl, v_sol-0.078, 0.132, text(L"v_{\mathrm{sol}}", 15, C_SOL))
    return pl
end

"""Panel (b): the homotopy family and the continuation path that walks along it."""
function panel_homotopy(grid, p, path)
    xl, yl = (-0.15, 0.4), (-0.75, 0.75)
    pl = plot(title = "(b)", xlabel = L"v_q\ \mathrm{(p.u.)}", ylabel = L"H(v_q,\lambda)",
              legend = :topleft, xlims = xl, ylims = yl)
    plot!(pl, collect(xl), [0.0, 0.0]; color = :black, linewidth = 1.6,
        #   label = L"R(v_q;\lambda) = 0",
          label="",)
    plot!(pl, grid, reduced_curve(grid, 0.0, p); color = C_PRE, label = "")
    for (j, λ) in enumerate(LAMBDAS[2:end-1])
        plot!(pl, grid, reduced_curve(grid, λ, p); color = C_MID, linestyle = :dash,
              linewidth = 2.6, label = j == 1 ? "Intermediate manifolds" : "")
    end
    plot!(pl, grid, reduced_curve(grid, 1.0, p); color = C_POST,
        #   label = L"\mathrm{Post\!-\!event\ manifold}\ (\lambda = 1)",
          label="",)

    # Continuation path: advance lambda (vertical, the residual the step creates),
    # then correct back to the axis (diagonal, the Newton correction).
    xs, ys = Float64[], Float64[]
    for k in 1:length(LAMBDAS)-1
        v_here, λ_next = path[k], LAMBDAS[k+1]
        r_here = R_reduced(v_here, λ_next, p; v_d0=1.0)[1]
        append!(xs, [v_here, v_here, path[k+1]])
        append!(ys, [0.0, r_here, 0.0])
    end
    plot!(pl, xs, ys; color = C_PATH, linewidth = 2.0, label = "Continuation path")
    scatter!(pl, xs, ys; marker = :circle, markersize = 5.5, color = C_ITER, label = "")
    scatter!(pl, [path[1]], [0.0]; marker = :utriangle, markersize = 13, color = C_INIT, label = "")
    scatter!(pl, [path[end]], [0.0]; marker = :star5, markersize = 15, color = C_SOL, label = "")
    # plot!(pl, [path[1], path[1]+0.030], [0.05, 0.30]; color = C_INIT, linewidth = 1.2, label = "")
    # annotate!(pl, path[1]+0.042, 0.38, text(L"v_{\mathrm{init}}", 15, C_INIT))
    # plot!(pl, [path[end], path[end]-0.045], [-0.05, -0.32]; color = C_SOL, linewidth = 1.2, label = "")
    # annotate!(pl, path[end]-0.062, -0.41, text(L"v_{\mathrm{sol}}", 15, C_SOL))
    return pl
end

"""Panel (c): what the geometry costs numerically."""
function panel_convergence(nr, hc, stage_bounds)
    pl = plot(title = "(c)", xlabel = "Total NR iterations",
              ylabel = L"\|G(\mathbf{z})\|_\infty", yscale = :log10,
              legend = :bottomleft, ylims = (1e-10, 1e8))
    floor_ = 1e-10
    # Each continuation stage restarts its own residual, so the homotopy trace is
    # sawtoothed by construction; the dotted lines mark the stage boundaries.
    for b in stage_bounds
        vline!(pl, [b]; color = :gray70, linewidth = 1.0, linestyle = :dot, label = "")
    end
    # logabsx = log10.(abs.(nr.residuals))
    # plot!(pl, 1:length(nr.residuals), max.(nr.residuals, floor_); color = C_POST,
    #       linestyle = :dash, marker = :circle, markersize = 7,
    #       label = "Direct Newton (diverges)")
    # plot!(pl, 1:20, logabsx[1:20]; color = C_POST,
    #       linestyle = :dash, marker = :circle, markersize = 7,
    #       label = "Direct Newton (diverges)")
    plot!(pl, 1:25, nr.residuals[1:25]; color = C_POST,
          linestyle = :dash, marker = :circle, markersize = 7,
          label = "Direct Newton (diverges)")      
    plot!(pl, 1:length(hc.residuals), max.(hc.residuals, floor_); color = C_MID,
          marker = :square, markersize = 4, label = "Convex homotopy")
    return pl
end

function main_plot()
    mkpath(OUT)
    gr(); style!()
    p = P
    # grid = collect(range(-0.26, 0.30; length = 601))
    grid = collect(range(-0.8, 0.5; length=601))
    # grid = collect(range(-0.9, 0.3; length=601))

    # Solve once with the six-variable solvers; the reduced curves are only for display.
    u_pre = pre_event_state(p)
    u_nr = copy(u_pre)
    nr = solve_newton!(u_nr, post_event(p), ADDRESS; tol=1e-9, max_iter=100,
                       always_new=true, model! = converter_homotopy!)
    u_hc = copy(u_pre)
    hc = solve_homotopy!(u_hc, p, ADDRESS; tol=1e-9, max_iter=100, Δλ=0.2,
                         vd_idx=VD, vq_idx=VQ, always_new=true,
                         model! = converter_homotopy!)
    @assert hc.converged "continuation failed; nothing to plot"
    path = copy(hc.vq_hist)

    # Cross-check the reduction against the six-variable path before drawing anything.
    for (k, λ) in enumerate(LAMBDAS)
        r = R_reduced(path[k], λ, p; v_d0=1.0)[1]
        @assert abs(r) < 1e-7 "reduced residual $(r) at λ=$(λ) disagrees with the solver"
    end

    roots1 = reduced_roots(grid, 1.0, p)
    v_sol = path[end]
    others = filter(r -> abs(r - v_sol) > 1e-5, roots1)
    v_spurious = isempty(others) ? NaN : others[argmin(abs.(others .- v_sol))]

    @printf("pre-event solution   v_q = %+.6f\n", path[1])
    @printf("post-event solution  v_q = %+.6f  (V = %.6f)\n", v_sol, voltage(u_hc))
    if !isnan(v_spurious)
        v_d_sp = v_d_of(v_spurious, 1.0, p)
        @printf("spurious branch      v_q = %+.6f  (V = %.6f)\n", v_spurious, hypot(v_d_sp, v_spurious))
    end

    pa = panel_event(grid, p, path[1], v_sol, v_spurious)
    pb = panel_homotopy(grid, p, path)
    # Stage boundaries: solve_homotopy! concatenates residuals across lambda stages.
    bounds = Float64[]; acc = 0
    for k in 1:length(hc.stage_iters)-1
        acc += hc.stage_iters[k]; push!(bounds, acc + 0.5)
    end
    pc = panel_convergence(nr, hc, bounds)

    two = plot(pa, pb; layout = (1, 2), size = (1350, 580), left_margin = 10Plots.mm, bottom_margin = 9Plots.mm, top_margin = 4Plots.mm, right_margin = 4Plots.mm)
    savefig(two, joinpath(OUT, "converter_manifold_reinit.pdf"))
    savefig(two, joinpath(OUT, "converter_manifold_reinit.png"))

    three = plot(pa, pb, pc; layout = (1, 3), size = (1850, 600), left_margin = 11Plots.mm, bottom_margin = 10Plots.mm, top_margin = 4Plots.mm, right_margin = 4Plots.mm)
    savefig(three, joinpath(OUT, "converter_manifold_reinit_3panel.pdf"))
    savefig(three, joinpath(OUT, "converter_manifold_reinit_3panel.png"))

    println("\nWrote figures to $(normpath(OUT))")
end

main_plot()
