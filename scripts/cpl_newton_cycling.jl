using LinearAlgebra, ForwardDiff, Plots, LaTeXStrings

# ══════════════════════════════════════════════════════════════════════════════
# Series RL circuit with switching shunt resistance and constant power load
#
#   V_s ──[ R ]──[ L ]──(node v)──┬──[ R_shunt ]
#                                  └──[ P_load = const ]
#
# DAE:
#   Differential:  L · di/dt = V_s - R·i - v
#   Algebraic:     0 = i - v/R_shunt - P₀/v     (KCL at node)
#
# At steady state:  i = (V_s - v) / R
# Substituting into KCL:
#   (V_s - v)/R = v/R_shunt + P₀/v
#   → v² (1/R + 1/R_shunt) - V_s/R · v + P₀ = 0
#
# This quadratic has two solutions (high-V and low-V) when the discriminant > 0.
# A switching event (R_shunt changes) shifts both roots. If the pre-event voltage
# lies outside the basin of convergence of the post-event high-V root, Newton
# cycles or diverges.
# ══════════════════════════════════════════════════════════════════════════════

# ── Circuit parameters ──
V_s   = 100.0       # source voltage (V)
R     = 1.0         # series resistance (Ω)
L     = 0.01        # inductance (H)
P0    = 50.0       # constant power load (W)

R_shunt_pre  = 50.0    # shunt resistance before event (Ω)
R_shunt_post = 0.300   # shunt resistance after event (Ω) — load shed

# ── Algebraic equation: g(v, i, R_sh) = i - v/R_sh - P₀/v ──
# At steady state i = (V_s - v)/R, so everything is in terms of v:
function g_ss(v, R_sh)
    i_ss = (V_s - v) / R
    return i_ss - v / R_sh - P0 / v
end

# Full DAE residual for the 2-state system [i, v]
function residual!(F, x, R_sh)
    i, v = x
    F[1] = V_s - R * i - v          # L·di/dt = ... (= 0 at steady state)
    F[2] = i - v / R_sh - P0 / v    # KCL (algebraic)
end

# ── Analytical roots of the steady-state quadratic ──
# v²·(1/R + 1/R_sh) - (V_s/R)·v + P₀ = 0
function analytical_roots(R_sh)
    a = 1/R + 1/R_sh
    b = -V_s / R
    c = P0
    disc = b^2 - 4*a*c
    if disc < 0
        return nothing, nothing  # no real solution
    end
    v_high = (-b + sqrt(disc)) / (2a)
    v_low  = (-b - sqrt(disc)) / (2a)
    return v_high, v_low
end

# ── Newton-Raphson solver with full iteration history ──
function newton_solve(x0, R_sh; tol=1e-10, max_iter=50)
    x = copy(x0)
    hist_x = [copy(x)]
    hist_res = Float64[]

    for k in 1:max_iter
        F = zeros(2)
        residual!(F, x, R_sh)
        push!(hist_res, norm(F))

        J = ForwardDiff.jacobian(z -> begin
            FF = similar(z)
            residual!(FF, z, R_sh)
            FF
        end, x)

        dx = J \ (-F)
        x .+= dx
        push!(hist_x, copy(x))

        norm(dx) < tol * (1 + norm(x)) && break
    end

    F_final = zeros(2)
    residual!(F_final, x, R_sh)
    push!(hist_res, norm(F_final))

    return (x=x, hist_x=hist_x, hist_res=hist_res,
            converged=hist_res[end] < tol)
end

# ══════════════════════════════════════════════════════════════════════════════
# 1. Compute pre-event and post-event equilibria
# ══════════════════════════════════════════════════════════════════════════════

v_hi_pre, v_lo_pre = analytical_roots(R_shunt_pre)
v_hi_post, v_lo_post = analytical_roots(R_shunt_post)

i_hi_pre = (V_s - v_hi_pre) / R
i_hi_post = (V_s - v_hi_post) / R

println("Pre-event:  v_high = $(round(v_hi_pre, digits=4)),  v_low = $(round(v_lo_pre, digits=4))")
println("Post-event: v_high = $(round(v_hi_post, digits=4)),  v_low = $(round(v_lo_post, digits=4))")
println("Jump in v: $(round(v_hi_pre, digits=4)) → $(round(v_hi_post, digits=4))")

# ══════════════════════════════════════════════════════════════════════════════
# 2. Newton from pre-event state into post-event system (NO reinitialization)
# ══════════════════════════════════════════════════════════════════════════════

x_pre = [i_hi_pre, v_hi_pre]   # pre-event equilibrium
r_no_reinit = newton_solve(x_pre, R_shunt_post; max_iter=50)

println("\n── Newton WITHOUT reinitialization ──")
println("Starting from pre-event: v = $(round(v_hi_pre, digits=4))")
println("Converged: $(r_no_reinit.converged)")
println("Final v: $(round(r_no_reinit.x[2], digits=4))")
println("Iterations: $(length(r_no_reinit.hist_res)-1)")

# ══════════════════════════════════════════════════════════════════════════════
# 3. Newton from a good initial guess (WITH reinitialization / near post-event)
# ══════════════════════════════════════════════════════════════════════════════

x_good = [i_hi_post, v_hi_post] .* 0.95  # slightly perturbed post-event solution
r_reinit = newton_solve(x_good, R_shunt_post; max_iter=50)

println("\n── Newton WITH reinitialization ──")
println("Starting from near post-event: v = $(round(x_good[2], digits=4))")
println("Converged: $(r_reinit.converged)")
println("Final v: $(round(r_reinit.x[2], digits=4))")
println("Iterations: $(length(r_reinit.hist_res)-1)")

# ══════════════════════════════════════════════════════════════════════════════
# 4. Figures
# ══════════════════════════════════════════════════════════════════════════════

# ── Figure 1: PV curves and Newton iteration path ──
v_range = range(5.0, 99.0, length=500)

# Power delivered to node as function of v:  P_node(v) = v · i = v·(V_s - v)/R
# Power consumed:  P_consumed(v) = v²/R_sh + P₀
P_available(v) = v * (V_s - v) / R
P_consumed_pre(v) = v^2 / R_shunt_pre + P0
P_consumed_post(v) = v^2 / R_shunt_post + P0

# Alternatively, plot g_ss(v) = 0 curve directly
g_pre  = [g_ss(v, R_shunt_pre)  for v in v_range]
g_post = [g_ss(v, R_shunt_post) for v in v_range]

lw = 1.5
ms = 5

p_pv = plot(v_range, g_pre, lw=lw, color=:gray, linestyle=:dash,
     label="Pre-event: g(v)", ylabel=L"g(v)", xlabel=L"v \; \mathrm{(V)}")
plot!(p_pv, v_range, g_post, lw=lw, color=:black, label="Post-event: g(v)")
hline!(p_pv, [0.0], color=:black, linestyle=:dot, lw=0.8, label="")

# Mark equilibria
scatter!(p_pv, [v_hi_pre], [0.0], color=:green, markersize=7, label="Pre-event equil.", markerstrokewidth=0)
scatter!(p_pv, [v_hi_post, v_lo_post], [0.0, 0.0], color=:blue,
         marker=[:circle :diamond], markersize=7, label="Post-event equil.", markerstrokewidth=0)

# Newton iteration path (voltage component only)
v_newton = [h[2] for h in r_no_reinit.hist_x]
n_show = min(length(v_newton), 20)
for k in 1:n_show-1
    vk = v_newton[k]
    gk = g_ss(vk, R_shunt_post)
    # vertical line from (vk, g(vk)) to (vk, 0), then horizontal to next iterate
    plot!(p_pv, [vk, vk], [0, gk], color=:red, lw=0.8, alpha=0.6, label=(k==1 ? "Newton path" : ""))
    plot!(p_pv, [vk, v_newton[k+1]], [0, 0], color=:red, lw=0.8, alpha=0.6, arrow=true, label="")
end

plot!(p_pv, framestyle=:box, grid=true, gridalpha=0.2,
     legend=:topright, legendfontsize=6, tickfontsize=8, guidefontsize=9)

# ── Figure 2: Convergence comparison ──
n1 = length(r_no_reinit.hist_res)
n2 = length(r_reinit.hist_res)

p_conv = plot(1:n1, r_no_reinit.hist_res,
     yscale=:log10, lw=lw, color=:red, linestyle=:dash,
     marker=:diamond, markersize=4, markerstrokewidth=0,
     label="Without reinitialization",
     xlabel="Newton iteration", ylabel=L"\| F(\mathbf{x}) \|_2")

plot!(p_conv, 1:n2, r_reinit.hist_res,
     lw=lw, color=:blue,
     marker=:circle, markersize=4, markerstrokewidth=0,
     label="With reinitialization")

hline!(p_conv, [1e-10], color=:black, linestyle=:dot, lw=0.8, label="")
plot!(p_conv, framestyle=:box, grid=true, gridalpha=0.2,
     legend=:right, legendfontsize=6, tickfontsize=8, guidefontsize=9)

# ── Combined figure ──
fig = plot(p_pv, p_conv, layout=(1, 2), size=(750, 320), dpi=300,
          margin=4Plots.mm, bottom_margin=5Plots.mm,
          plot_title="")

display(fig)
# savefig(fig, "cpl_newton_cycling.pdf")


i_cov(e) = 5.0 * tanh(0.5*e + 0.5)

e = -7.0:0.01:7.0

i = @. i_cov(e)

plot(e, i)

f(x) = 5.0 * tanh(0.5*x + 0.5)

df(x) = 2.5 * sech(0.5*x + 0.5)^2

x_k(x) = x - f(x)/df(x)

j = 0
x0 = 2.0
while j<10
    x1 = x_k(x0)
    x0 = x1
    j+=1
    @show x0
end 


i_cov(e) = 5.0 * tanh(0.5*e - 2.5)

e = -7.0:0.01:7.0

i = @. i_cov(e)
plot(e,i)