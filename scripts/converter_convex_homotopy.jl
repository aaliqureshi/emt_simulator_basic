# Convex homotopy re-initialization for the GFM converter tutorial case.
#
#   julia --project=. scripts/converter_convex_homotopy.jl
#
# The event switches the converter from voltage-forming (GFM) operation to the
# current-limited mode and applies a fault conductance g_f. That changes the
# algebraic constraint manifold, so the inherited pre-event algebraic state is no
# longer consistent and the DAE cannot be resumed from it.
#
# Re-initialization deforms the pre-event manifold into the post-event one,
#
#     H(y, λ) = (1-λ)·G_GFM(y) + λ·G_GFL(y),
#
# and tracks the solution branch from λ=0 to λ=1. Unlike the natural-parameter
# ramp of g_f, λ=0 here is the pre-event manifold *exactly*, so the inherited
# state is already a root and no auxiliary solve is needed to start the path.
#
# The differential state e_q is held fixed across the event (it is continuous),
# which is what the solvers in algebraic_solvers.jl do: they freeze the leading
# n_mech entries of u and solve only the algebraic block.

using LinearAlgebra, Printf
using Barq
using MyDiffEq

# u = (e_q, i_ld, i_lq, v_d, v_q, i_cd, i_cq); only e_q is differential.
const EQ, ILD, ILQ, VD, VQ, ICD, ICQ = 1, 2, 3, 4, 5, 6, 7

# p_base = (v_ref, v_slack, r_line, x_line, ω, tq, ki, kp, i_L, g_f, i_max)
# const P = (1.0, 1.0, 0.08, 0.35, 2π*60, 0.2, 6.5, 4.7, 0.32, 1.119, 1.43)

# const P = (1.0, 1.0, 0.08, 0.35, 2π*60, 0.2, 6.5, 4.7, 0.32, 0.6, 1.43)
# const P = (1.0, 1.0, 0.08, 0.35, 2π*60, 0.2, 6.5, 4.7, 0.32, 0.6, 0.9)
# const P = (1.0, 1.0, 0.08, 0.35, 2π*60, 0.2, 6.5, 8.0, 0.2, 1.70, 1.60)
const P = (1.0, 1.0, 0.08, 0.35, 2π*60, 0.2, 6.5, 4.9, 0.32, 1.144, 1.43)

# The solvers derive the algebraic block from address["delta"] and address["omega"].
# One differential state (e_q) => alg_idx = 2:7.
const ADDRESS = Dict("delta" => [1], "omega" => Int[])

"""
    _load_currents(v_d, v_q, i_L)

Constant-current load with conversion to constant impedance below 0.7 p.u.
"""
function _load_currents(v_d, v_q, i_L)
    V = sqrt(v_d^2 + v_q^2)
    m = max(V, 0.7)
    return i_L*v_d/m, i_L*v_q/m
end

"""
    converter_homotopy!(du, u, p, t)

Residual of the convex homotopy H(y,λ) = (1-λ)·G_GFM + λ·G_GFL.

`p` is the base parameter tuple with λ appended as the final entry, matching the
convention used by `solve_homotopy!` and `solve_adaptive_homotopy!`
(`p = (p_base..., λ)`). λ=0 recovers the pre-event GFM model and λ=1 the
post-event current-limited model, so this single function also serves as the
plain GFM and GFL residual.

Rows 2-3 (line dynamics) are common to both modes and are therefore unaffected
by λ. Rows 4-5 blend in the fault current, and rows 6-7 blend the voltage-forming
constraints into the reactive-current command and the circular current limit.
"""
function converter_homotopy!(du, u, p, t)
    e_q, i_ld, i_lq, v_d, v_q, i_cd, i_cq = u
    v_ref, v_slack, r_line, x_line, ω, tq, ki, kp, i_L, g_f, i_max = p[1:11]
    λ = p[12]

    V = sqrt(v_d^2 + v_q^2)
    i_Ld, i_Lq = _load_currents(v_d, v_q, i_L)

    # Differential state: reactive-voltage control lag. Frozen during re-init.
    du[EQ] = (ki*(v_ref - V) - e_q)/tq

    # Line equations: identical in both modes.
    du[ILD] = (v_d - v_slack - r_line*i_ld + x_line*i_lq)*ω/x_line
    du[ILQ] = (v_q - r_line*i_lq - x_line*i_ld)*ω/x_line

    # KCL at the PCC. The fault branch exists only in the post-event model, so
    # the convex combination scales it by λ.
    du[VD] = i_cd - i_ld - i_Ld - λ*g_f*v_d
    du[VQ] = i_cq - i_lq - i_Lq - λ*g_f*v_q

    # Mode-defining rows. λ=0: v_d = v_ref, v_q = 0 (voltage forming).
    #                    λ=1: reactive-current command + circular current limit.
    du[ICD] = (1-λ)*(v_d - v_ref) + λ*(i_cq - (e_q + kp*(v_ref - V)))
    du[ICQ] = (1-λ)*v_q          + λ*(i_cd^2 + i_cq^2 - i_max^2)

    return nothing
end

"""Pre-event (λ=0) and post-event (λ=1) parameter tuples."""
pre_event(p)  = (p..., 0.0)
post_event(p) = (p..., 1.0)

"""
    pre_event_state(p)

Closed-form GFM equilibrium. With v_d = v_ref, v_q = 0 the line equations give
i_ld, i_lq directly and KCL gives the converter currents.
"""
function pre_event_state(p)
    v_ref, v_slack, r_line, x_line = p[1], p[2], p[3], p[4]
    i_L = p[9]
    z2 = r_line^2 + x_line^2
    i_ld = r_line*(v_ref - v_slack)/z2
    i_lq = -x_line*(v_ref - v_slack)/z2
    i_Ld, i_Lq = _load_currents(v_ref, 0.0, i_L)
    return [0.0, i_ld, i_lq, v_ref, 0.0, i_ld + i_Ld, i_lq + i_Lq]
end

"""Infinity norm of the algebraic residual (rows 2:7) at parameters p."""
function algebraic_residual(u, p)
    du = zeros(eltype(u), 7)
    converter_homotopy!(du, u, p, 0.0)
    return norm(du[2:end], Inf)
end

voltage(u) = hypot(u[VD], u[VQ])
current(u) = hypot(u[ICD], u[ICQ])

"""
    converter_current(v_d, v_q, p, λ)

Converter current magnitude implied by (v_d, v_q) at homotopy level λ. The line
equations and PCC KCL are linear in (i_ld, i_lq, i_cd, i_cq) given the voltage,
so they can be solved in closed form; this holds exactly at every converged
continuation stage and is used only for reporting.
"""
function converter_current(v_d, v_q, p, λ)
    v_slack, r_line, x_line, i_L, g_f = p[2], p[3], p[4], p[9], p[10]
    z2 = r_line^2 + x_line^2
    i_ld = (r_line*(v_d - v_slack) + x_line*v_q)/z2
    i_lq = (-x_line*(v_d - v_slack) + r_line*v_q)/z2
    i_Ld, i_Lq = _load_currents(v_d, v_q, i_L)
    return hypot(i_ld + i_Ld + λ*g_f*v_d, i_lq + i_Lq + λ*g_f*v_q)
end

function main()
    p_pre, p_post = pre_event(P), post_event(P)

    u_pre = pre_event_state(P)
    @printf("Pre-event GFM equilibrium:  V = %.6f p.u., |i_c| = %.6f p.u.\n",
            voltage(u_pre), current(u_pre))
    @printf("  |G_GFM(u_pre)|_inf = %.3e   (inherited state is on the pre-event manifold)\n",
            algebraic_residual(u_pre, p_pre))
    @printf("  |G_GFL(u_pre)|_inf = %.3e   (and off the post-event manifold)\n\n",
            algebraic_residual(u_pre, p_post))

    # ── Direct Newton re-initialization from the inherited state ──
    u_nr = copy(u_pre)
    nr = solve_newton!(u_nr, p_post, ADDRESS; tol=1e-9, max_iter=100,
                       always_new=true, model! = converter_homotopy!)
    @printf("Direct Newton re-init : converged = %-5s  iters = %3d  final |G|_inf = %.3e\n",
            nr.converged, nr.iters, last(nr.residuals))

    # ── Convex homotopy re-initialization ──
    # λ=0 is the pre-event manifold, so the inherited state is already a root and
    # the first continuation stage costs a single residual evaluation.
    u_hc = copy(u_pre)
    hc = solve_homotopy!(u_hc, P, ADDRESS; tol=1e-9, max_iter=100, Δλ=0.2,
                         vd_idx=VD, vq_idx=VQ, always_new=true,
                         model! = converter_homotopy!)
    @printf("Convex homotopy re-init: converged = %-5s  iters = %3d  λ_failed = %s\n",
            hc.converged, hc.total_iters, string(hc.λ_failed))

    if !hc.converged
        @printf("\nHomotopy path terminated at λ = %s: the post-event manifold is\n", string(hc.λ_failed))
        println("unreachable along this branch. Reduce the event severity or use arclength tracking.")
        return
    end

    # The λ=0 stage costs a single residual check because the inherited state is
    # already a root there; the natural-parameter ramp of g_f needs a real solve.
    println("\nContinuation path:")
    println("    λ        V (p.u.)   |i_c| (p.u.)")
    for (i, λ) in enumerate(0.0:0.2:1.0)
        v_d, v_q = hc.vd_hist[i], hc.vq_hist[i]
        @printf("  %.1f     %9.6f   %9.6f\n",
                λ, hypot(v_d, v_q), converter_current(v_d, v_q, P, λ))
    end

    # ── Consistency of the re-initialized point ──
    @printf("\nRe-initialized state: V = %.6f p.u., |i_c| = %.6f p.u. (i_max = %.2f)\n",
            voltage(u_hc), current(u_hc), P[11])
    @printf("  |G_GFL(u_hc)|_inf = %.3e\n", algebraic_residual(u_hc, p_post))
    @printf("  e_q preserved across the event: %.3e\n", abs(u_hc[EQ] - u_pre[EQ]))
    @assert u_hc[EQ] == u_pre[EQ] "differential state must be continuous across the event"
    @assert algebraic_residual(u_hc, p_post) < 1e-8 "re-initialized point is not consistent"

    # ── Resume the DAE integration from the consistent point ──
    mass = zeros(7, 7)
    mass[EQ, EQ] = 1.0
    for h in (1e-4, 1e-5, 1e-6)
        prob = MyDiffEq.ODEProblem(converter_homotopy!, copy(u_hc), (0.0, 0.01), p_post, mass)
        sol = MyDiffEq.Solve(prob, h; method=:Euler, adaptive=false, always_new=true)
        ok = sol.retcode == :Success
        @printf("Resumed integration h = %4d μs: retcode = %-8s", round(Int, h*1e6), string(sol.retcode))
        ok ? @printf("  V(10 ms) = %.6f  |i_c| = %.6f\n", voltage(sol.u[end]), current(sol.u[end])) : println()
    end

    # ── Same integration attempted from the inconsistent inherited state ──
    println()
    for h in (1e-4, 1e-5, 1e-6)
        prob = MyDiffEq.ODEProblem(converter_homotopy!, copy(u_pre), (0.0, 0.01), p_post, mass)
        sol = MyDiffEq.Solve(prob, h; method=:Euler, adaptive=false, always_new=true)
        @printf("No re-init,    h = %4d μs: retcode = %s\n", round(Int, h*1e6), string(sol.retcode))
    end
end

# Only run when executed as a script, so the model above can be `include`d by the
# plotting code without triggering the solve.
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
