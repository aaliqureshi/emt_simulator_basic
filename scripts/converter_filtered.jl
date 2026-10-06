# Nine-state GFM converter with the coupling filter represented in phasor form.
#
# This is a SEPARATE model from the seven-state `converter_homotopy!` in
# converter_convex_homotopy.jl, which is left untouched.
#
#   slack --(r_l + j x_l)-- PCC --(r_c + j x_c)-- converter terminal
#                            |
#                     load i_L + fault g_f
#
# The voltage-forming control regulates the PCC, not the converter terminal, so
# the terminal voltage w floats to whatever the filter requires. That matters:
# regulating the terminal instead would starve the fault current and the current
# limit would never bind (see the note at the bottom of this file).
#
# Because w appears ONLY in the two filter equations, the system is block
# triangular: the six equations in (i_ld, i_lq, v_d, v_q, i_cd, i_cq) are exactly
# those of the seven-state model, and the filter rows then give w explicitly.
# The filter is therefore physically meaningful but numerically inert for the
# re-initialization problem — `assert_matches_seven_state` checks this on demand.

using LinearAlgebra, Printf
using Barq

include(joinpath(@__DIR__, "converter_convex_homotopy.jl"))

# u = (e_q, i_ld, i_lq, v_d, v_q, i_cd, i_cq, w_d, w_q); only e_q is differential.
const EQ_F, ILD_F, ILQ_F, VD_F, VQ_F, ICD_F, ICQ_F, WD_F, WQ_F = 1, 2, 3, 4, 5, 6, 7, 8, 9

# p_base = (v_ref, v_slack, r_line, x_line, ω, tq, ki, kp, i_L, g_f, i_max, r_c, x_c)
const P_F = (P..., 0.015, 0.15)

const ADDRESS_F = Dict("delta" => [1], "omega" => Int[])

"""
    converter_filtered!(du, u, p, t)

Convex homotopy residual H(y,λ) = (1-λ)·G_GFM + λ·G_GFL for the nine-state model.
`p` is the 13-entry base tuple with λ appended, matching the convention used by
`solve_homotopy!` (`p = (p_base..., λ)`).
"""
function converter_filtered!(du, u, p, t)
    e_q, i_ld, i_lq, v_d, v_q, i_cd, i_cq, w_d, w_q = u
    v_ref, v_slack, r_line, x_line, ω, tq, ki, kp, i_L, g_f, i_max, r_c, x_c = p[1:13]
    λ = p[14]

    V = sqrt(v_d^2 + v_q^2)
    i_Ld, i_Lq = _load_currents(v_d, v_q, i_L)

    du[EQ_F]  = (ki*(v_ref - V) - e_q)/tq

    du[ILD_F] = (v_d - v_slack - r_line*i_ld + x_line*i_lq)*ω/x_line
    du[ILQ_F] = (v_q           - r_line*i_lq - x_line*i_ld)*ω/x_line

    du[VD_F]  = i_cd - i_ld - i_Ld - λ*g_f*v_d
    du[VQ_F]  = i_cq - i_lq - i_Lq - λ*g_f*v_q

    # Coupling filter, converter terminal -> PCC. Not scaled by ω/x_c: the scaling
    # is arbitrary row weighting and a small x_c would wreck the residual norms.
    du[ICD_F] = w_d - v_d - r_c*i_cd + x_c*i_cq
    du[ICQ_F] = w_q - v_q - r_c*i_cq - x_c*i_cd

    # Mode rows act on the PCC, exactly as in the seven-state model.
    du[WD_F]  = (1-λ)*(v_d - v_ref) + λ*(i_cq - (e_q + kp*(v_ref - V)))
    du[WQ_F]  = (1-λ)*v_q           + λ*(i_cd^2 + i_cq^2 - i_max^2)

    return nothing
end

pre_event_f(p)  = (p..., 0.0)
post_event_f(p) = (p..., 1.0)

"""Pre-event equilibrium: the seven-state solution, with w recovered from the filter."""
function pre_event_state_filtered(p)
    u7 = pre_event_state(p[1:11])
    r_c, x_c = p[12], p[13]
    w_d = u7[VD] + r_c*u7[ICD] - x_c*u7[ICQ]
    w_q = u7[VQ] + r_c*u7[ICQ] + x_c*u7[ICD]
    return vcat(u7, w_d, w_q)
end

voltage_f(u)  = hypot(u[VD_F], u[VQ_F])
current_f(u)  = hypot(u[ICD_F], u[ICQ_F])
terminal_f(u) = hypot(u[WD_F], u[WQ_F])

"""Infinity norm of the algebraic residual (rows 2:9)."""
function algebraic_residual_f(u, p)
    du = zeros(eltype(u), 9)
    converter_filtered!(du, u, p, 0.0)
    return norm(du[2:end], Inf)
end

function mass_matrix_f()
    m = zeros(9, 9)
    m[EQ_F, EQ_F] = 1.0
    return m
end

"""
    assert_matches_seven_state(p_f)

Check that adding the filter changes nothing about the re-initialization problem:
the continuation path and the re-initialized state must agree with the seven-state
model to machine precision in states 1..7.
"""
function assert_matches_seven_state(p_f = P_F; tol = 1e-12)
    p7 = p_f[1:11]

    u9 = copy(pre_event_state_filtered(p_f))
    hc9 = solve_homotopy!(u9, p_f, ADDRESS_F; tol=1e-9, max_iter=100, Δλ=0.2,
                          vd_idx=VD_F, vq_idx=VQ_F, always_new=true,
                          model! = converter_filtered!)
    u7 = copy(pre_event_state(p7))
    hc7 = solve_homotopy!(u7, p7, ADDRESS; tol=1e-9, max_iter=100, Δλ=0.2,
                          vd_idx=VD, vq_idx=VQ, always_new=true,
                          model! = converter_homotopy!)

    @assert hc9.converged && hc7.converged "continuation failed"
    @assert norm(u9[1:7] - u7) < tol "nine-state solution differs from seven-state"
    @assert norm(hc9.vd_hist - hc7.vd_hist) < tol "continuation paths differ (v_d)"
    @assert norm(hc9.vq_hist - hc7.vq_hist) < tol "continuation paths differ (v_q)"
    @assert hc9.total_iters == hc7.total_iters "iteration counts differ"
    return (; agree = norm(u9[1:7] - u7), iters = hc9.total_iters)
end

if abspath(PROGRAM_FILE) == @__FILE__
    r = assert_matches_seven_state()
    @printf("nine-state matches seven-state to %.2e in states 1..7 (%d Newton iterations)\n",
            r.agree, r.iters)
    u = copy(pre_event_state_filtered(P_F))
    @printf("pre-event : V_pcc = %.6f  |i_c| = %.6f  |w| = %.6f\n",
            voltage_f(u), current_f(u), terminal_f(u))
    solve_homotopy!(u, P_F, ADDRESS_F; tol=1e-9, max_iter=100, Δλ=0.2,
                    vd_idx=VD_F, vq_idx=VQ_F, always_new=true, model! = converter_filtered!)
    @printf("post-event: V_pcc = %.6f  |i_c| = %.6f  |w| = %.6f\n",
            voltage_f(u), current_f(u), terminal_f(u))
end
