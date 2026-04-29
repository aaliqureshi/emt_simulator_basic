using Plots, Printf
using LaTeXStrings



const V_ref  = 1.0
const R_g    = 0.1
const P_load = 0.9
const I_max  = 3.0

res_all = []
v_all=[]

function gfm_residual(V)
    return (V_ref - V)/R_g - P_load/V
end

function gfm_jacobian(V)
    return -1.0/R_g + P_load/V^2
end

function gfl_residual(V)
    return I_max - P_load/V
end

function gfl_jacobian(V)
    return P_load/V^2
end

# Newton homotopy anchored at the inherited pre-switch state V0
function homotopy_residual(V, λ, V0)
    return gfl_residual(V) - (1 - λ)*gfl_residual(V0)
end

function homotopy_jacobian(V, λ, V0)
    return gfl_jacobian(V)
end


function solve_nr(residual_func, jacobian_func, initial_guess; max_iter=20, tol=1e-10)
    V = initial_guess
    for k in 1:max_iter
        if V <= 0 || !isfinite(V)
            return NaN
        end

        r = residual_func(V)
        J = jacobian_func(V)

        push!(v_all, V)
        push!(res_all, r)

        @printf("Iter %2d | V = %.8f | r = %.4e | J = %.4e\n", k, V, r, J)

        if abs(r) < tol
            # push!(v_all, V)
            # push!(res_all, r)
            return V
        end

        if abs(J) < 1e-12 || !isfinite(J)
            return NaN
        end

        ΔV = -r/J

        # optional damping
        α = 1.0
        Vtrial = V + α*ΔV
        # while Vtrial <= 0
        #     α *= 0.5
        #     if α < 1e-8
        #         return NaN
        #     end
        #     Vtrial = V + α*ΔV
        # end

        V = Vtrial
    end
    return NaN
end

function solve_homotopy(V0; nsteps=5)
    V = V0
    λvals = range(0.0, 1.0, length=nsteps+1)

    for λ in λvals
        f(Vx) = homotopy_residual(Vx, λ, V0)
        J(Vx) = homotopy_jacobian(Vx, λ, V0)

        Vnew = solve_nr(f, J, V; max_iter=20, tol=1e-10)

        if isnan(Vnew)
            @printf("Homotopy failed at λ = %.3f\n", λ)
            return NaN
        end

        @printf("λ = %.3f  ->  V = %.8f\n", λ, Vnew)
        V = Vnew
    end

    return V
end

V_hom = solve_homotopy(0.9, nsteps=4)


gr()
# ============================================================
# Colors
# ============================================================
c_gfm   = RGB(0.00, 0.32, 0.84)   # deep blue
c_gfl   = RGB(0.90, 0.52, 0.00)   # warm orange
c_homo  = RGB(0.00, 0.62, 0.62)   # teal
c_path  = RGB(0.45, 0.25, 0.15)   # brown
c_init  = RGB(0.00, 0.70, 0.75)   # cyan
c_guess   = RGB(0.00, 0.75, 0.75)   # cyan-ish
c_sol   = RGB(0.90, 0.0, 0.9)   # cyan-ish
c_iter  = RGB(0.86, 0.25, 0.25)   # red for Newton points
c_zero  = :black

size = (1200, 500)

default(
    fontfamily = "Computer Modern",
    linewidth = 3.5,
    markersize = 9,
    legendfontsize = 13,
    guidefontsize = 14,
    tickfontsize = 14,
    titlefontsize = 14,
    dpi = 300,
    # size = (1200, 500),
    size = size,
    # margin = 6Plots.mm
)

v = @. 0.1:0.001:1.0
# v = @. 0.1:0.01:3.0
V0 = 0.9
r1 = @. homotopy_residual(v, 0.0, V0)

r2 = @. homotopy_residual(v, 0.2, V0)
r3 = @. homotopy_residual(v, 0.4, V0)
r4 = @. homotopy_residual(v, 0.6, V0)
r5 = @. homotopy_residual(v, 0.8, V0)
r6 = @. homotopy_residual(v, 1.0, V0)

p1 = plot(
    v, gfm_residual.(v),
    color = c_gfm,
    # label = L"\mathrm{GFM\ residual}",
    # label = L"\mathrm{Pre-event\ manifold}",
    label="Pre-event manifold",
    # xlabel = L"V\ (\mathrm{p.u})",
    xlabel = "V (p.u)",
    ylabel = L"H(V,\lambda)",
    framestyle = :box,
    legend = :bottomright,
    grid = false,
    left_margin = 1Plots.mm,
    right_margin = 1Plots.mm,
    bottom_margin = 1Plots.mm,
    top_margin = 1Plots.mm,
)

plot!(
    p1, v, gfl_residual.(v),
    color = c_gfl,
    # label = L"\mathrm{GFL\ residual}\ "
    # label = L"\mathrm{Current-limiting modeGFL\ residual}\ "
    label = "Post-event manifold"
)

hline!(
    p1, [0.0],
    color = c_zero,
    linestyle = :solid,
    linewidth = 1.6,
    label = L"R(V)=0"
)

# plot!(p1, [0.9, 0.9], 
#           [0.0, gfl_residual(0.9)],
#           linestyle=:dash,
#           linewidth = 2.0,
#           label="",
#           color=:orangered,
#       )

scatter!(
    p1, [0.9], [0.0],
    marker = (:utriangle, 13, stroke(1.2, :black)),
    color = c_init,
    label = "Pre-event solution",
    # zorder = 5
)

# scatter!(
#     p1, [0.9], [gfl_residual(0.9)],
#     marker = (:diamond, 11, stroke(1.2, :black)),
#     color = c_guess,
#     label = "",
#     # zorder = 5
# )

scatter!(
    p1, [0.3], [0.0],
    marker = (:star5, 13, stroke(1.2, :black)),
    color = c_sol,
    label = "Post-event solution",
    # zorder = 5
)

# vline!(p1, [0.9])

x_init = 0.90    # pre-switch GFM solution; used as post-switch initial guess
x_sol  = 0.30    # post-switch consistent solution

annotate!(p1, x_init + 0.015, 0.73, text(L"V_{\mathrm{init}}", 21))
annotate!(p1, x_sol - 0.010, 0.73, text(L"V_{\mathrm{sol}}", 21))

plot!(
    p1,
    # title = "(a) Mode switch creates a difficult re-initialization problem",
    xlims = (0.1, 1.00),
    # ylims = (-3.5, 4.2)
    ylims = (-3.1, 4.2)
)

# ============================================================
# Panel (b): homotopy continuation
# ============================================================
p2 = plot(
    v, gfm_residual.(v),
    color = c_gfm,
    # label = L"\mathrm{GFM\ residual}",
    label="Pre-event manifold",
    # xlabel = L"V\ (\mathrm{p.u})",
    xlabel = "V (p.u)",
    ylabel = L"H(V,\lambda)",
    framestyle = :box,
    legend = :bottomright,
    grid = false,
    left_margin = 1Plots.mm,
    right_margin = 1Plots.mm,
    bottom_margin = 2Plots.mm,
    top_margin = 2Plots.mm,
    
)

# λvals = [0.0, 0.20, 0.4, 0.60, 0.8]
λvals = [0.0, 0.25, 0.5, 0.75]

for (k, λ) in enumerate(λvals)
    lbl = k == 1 ? "Intermediate manifolds" : ""
    plot!(
        p2, v, homotopy_residual.(v, λ, 0.9),
        color = c_homo,
        linestyle = :dash,
        linewidth = 2.5,
        alpha = 0.75,
        label = lbl
    )
end

plot!(
    p2, v, gfl_residual.(v),
    color = c_gfl,
    # label = L"\mathrm{GFL\ residual}\ (\lambda=1.0\))"
    label="Post-event manifold (λ = 1)"
)

hline!(
    p2, [0.0],
    color = c_zero,
    linestyle = :solid,
    linewidth = 1.6,
    label = L"R(V\; ; λ)=0"
)

# Draw path between tracked roots
plot!(
    p2,
    # vcat(0.3, tracked_x, x_sol),
    v_all,
    res_all,
    # zeros(length(vcat(x_init, tracked_x, x_sol))),
    # color = c_path,
    # color = :sienna,
    color=:grey,
    linewidth = 2.4,
    marker = (:none),
    label = "Continuation path",
    # zorder = 4
)

scatter!(p2, v_all, res_all, 
    #   color=c_iter,
    color=:tomato,
      marker = (:circle, 5, stroke(0.8, :black)),
      alpha = 0.75,
      label = "")

# Continuation path
scatter!(
    p2, [0.9], [0.0],
    marker = (:utriangle, 13, stroke(1.2, :black)),
    color = c_init,
    label = "",
    # zorder = 6
)

scatter!(
    p2, [0.3], [0.0],
    marker = (:star5, 13, stroke(1.2, :black)),
    color = c_sol,
    label = "",
    # zorder = 6
)





# plot!(
#     p1,
#     # title = "(a) Mode switch creates a difficult re-initialization problem",
#     xlims = (0.1, 1.00),
#     ylims = (-3.5, 4.2)
# )
plot!(
    p2,
    title = "",
    xlims = (0.10, 1.00),
    ylims = (-3.1, 2.2)
)

annotate!(p2, x_init + 0.011, -0.4, text(L"V_{\mathrm{init}}", 21))
annotate!(p2, x_sol-0.009, 0.4, text(L"V_{\mathrm{sol}}", 21))
# ============================================================
# Combine panels
# ============================================================
plt = plot(
    p1, p2,
    layout = (1, 2),
    # size = (1300, 510)
    size = size,

    )




# savefig(plt, "homotopy_reinitialization_publication_ready.png")
# savefig(plt, "homotopy_reinitialization_converter_font5.pdf") 

display(plt)

####### Newton Divergnec
using Plots
using LaTeXStrings

gr()

# ------------------------------------------------------------
# Data
# ------------------------------------------------------------
iters = 1:5
xvals = [0.9, -0.9, -4.5, -76.5, -19660.5]
resid = [2.00, 4.00, 3.20, 3.01, 3.00]

# Useful transformed quantity for showing divergence clearly
logabsx = log10.(abs.(xvals))

size = (520, 380)

# tosave = false
tosave = true

# ------------------------------------------------------------
# Publication-style defaults
# ------------------------------------------------------------
default(
    fontfamily = "Computer Modern",
    linewidth = 3.0,
    markersize = 6,
    markerstrokewidth = 0.,
    legendfontsize = 13,
    guidefontsize = 14,
    tickfontsize = 14,
    titlefontsize = 13,
    dpi = 300,
    # size = (1000, 420),
    size = size,
    # margin = 5Plots.mm
)

# ------------------------------------------------------------
# Colors
# ------------------------------------------------------------
c_iter = RGB(0.78, 0.16, 0.16)   # deep red
c_res  = RGB(0.00, 0.45, 0.74)   # blue
c_ref  = :black

# ------------------------------------------------------------
# Subplot 1: iterate magnitude growth
# Plot log10(|x_n|) so early and late iterations are both visible
# ------------------------------------------------------------
p1 = plot(
    iters, logabsx,
    # iters, xvals,
    # yscale=:log10,
    color = c_iter,
    marker = :circle,
    xlabel = "NR Iteration (m)",
    ylabel = L"\log_{10}(|V^m|)",
    # title = "(a) Newton iterate magnitude",
    # label = L"\log_{10}(|V_k|)",
    label="",
    framestyle = :box,
    legend = :topleft,
    grid = true,
    gridalpha = 0.1,
    # xticks = collect(iters)
    left_margin = 1Plots.mm,
    right_margin = 1Plots.mm,
    bottom_margin = 1Plots.mm,
    top_margin = 3Plots.mm,
)

# Optional annotations for exact x_n values
for (n, y, xv) in zip(iters, logabsx, xvals)
    annotate!(p1, n, y + 0.3, text("$(xv)", 13))
end

display(p1)

tosave && savefig(p1, "converter_newton_divergence_iterates_single.pdf")
# ------------------------------------------------------------
# Subplot 2: residual history
# ------------------------------------------------------------
p2 = plot(
    iters, resid,
    color = c_res,
    marker = :circle,
    xlabel = "NR Iteration (m)",
    ylabel = L"|R(V^m)|",
    # title = "(b) Residual history",
    # label = L"|f(x_n)|",
    label = "",
    framestyle = :box,
    legend = :topright,
    grid = true,
    gridalpha = 0.1,
    xticks = collect(iters),
    left_margin = 1Plots.mm,
    right_margin = 1Plots.mm,
    bottom_margin = 1Plots.mm,
    top_margin = 2Plots.mm,
    size = size,
)

# Reference line at zero residual
# hline!(
#     p2, [0.0],
#     color = c_ref,
#     linestyle = :dash,
#     linewidth = 1.2,
#     label = "Zero residual"
# )

# ------------------------------------------------------------
# Combine subplots
# ------------------------------------------------------------
# plt = plot(
#     p1, p2,
#     layout = (1, 2),
#     size = (1000, 380)
# )

# savefig(plt, "newton_divergence_two_subplots.png")
# savefig(plt, "converter_newton_divergence_font4.pdf")

display(p2)

tosave && savefig(p2, "converter_newton_divergence_residuals_single.pdf")

