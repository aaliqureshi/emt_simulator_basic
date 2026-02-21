module Homotopy

export solve_algebraic!, pseudo_arclength!

include("dynamic_sim.jl"); using .DynamicSim

using ForwardDiff
using LinearAlgebra

function solve_algebraic!(u, p, mass_matrix, address; tol=1e-6, max_iter=1000)
    n = length(u)
    # Freeze only mechanical states (delta, omega); solve everything else
    n_mech = length(address["delta"]) + length(address["omega"])
    alg_idx = n_mech+1:n

    for k in 1:max_iter
        # Full residual evaluation (only algebraic part matters)
        du = zeros(n)
        solve_dynamic_sim!(du, u, p, 0.0)
        g = du[alg_idx]

        # Jacobian of algebraic equations w.r.t. algebraic variables only
        jac_alg = ForwardDiff.jacobian(y -> begin
            T = eltype(y)
            u_tmp = Vector{T}(u)
            u_tmp[alg_idx] .= y
            du_tmp = similar(u_tmp)
            solve_dynamic_sim!(du_tmp, u_tmp, p, 0.0)
            du_tmp[alg_idx]
        end, u[alg_idx])

        # Newton correction
        delta_y = jac_alg \ (-g)
        u[alg_idx] .+= delta_y

        if norm(delta_y, 2) < tol * (1 + norm(u[alg_idx], 2))
            return true, k
        end
    end
    return false, max_iter
end

"""
    pseudo_arclength!(u₁, λ₁, u₂, λ₂, p_base, mass_matrix, address; kwargs...)

Pseudo-arclength continuation for tracing the solution curve g(u, λ) = 0
past turning points (saddle-node bifurcations).

Uses true tangent computation (SVD of bordered Jacobian [dg/du | dg/dλ])
instead of secant approximation for robustness near turning points.

Takes two pre-converged points (u₁, λ₁) and (u₂, λ₂).

Returns:
    (path_λ, path_u, success)
"""
function pseudo_arclength!(u₁, λ₁, u₂, λ₂, p_base, mass_matrix, address;
                           λ_target=1.0, Δs=nothing,
                           tol=1e-6, max_newton=50, max_steps=500,
                           Δs_max=0.05, Δs_min=1e-12)

    n = length(u₁)
    n_mech = length(address["delta"]) + length(address["omega"])
    alg_idx = n_mech+1:n
    n_alg = length(alg_idx)

    make_p(λ) = (p_base..., λ)

    # Evaluate algebraic residual
    function eval_g(u_eval, λ_eval)
        du = zeros(n)
        solve_dynamic_sim!(du, u_eval, make_p(λ_eval), 0.0)
        return du[alg_idx]
    end

    # Full Jacobian [dg/du_alg | dg/dλ] via ForwardDiff (AD for both u and λ)
    function eval_full_jac(u_eval, λ_eval)
        z0 = vcat(u_eval[alg_idx], λ_eval)
        ForwardDiff.jacobian(z -> begin
            T = eltype(z)
            u_tmp = Vector{T}(u_eval)
            u_tmp[alg_idx] .= z[1:n_alg]
            λ_val = z[end]
            du_tmp = similar(u_tmp)
            solve_dynamic_sim!(du_tmp, u_tmp, make_p(λ_val), zero(T))
            du_tmp[alg_idx]
        end, z0)
    end

    # Tangent to solution curve via SVD of [dg/du_alg | dg/dλ]
    function compute_tangent(u_eval, λ_eval, direction_ref=nothing)
        G = eval_full_jac(u_eval, λ_eval)   # n_alg × (n_alg+1)
        F = svd(G)
        tangent = F.V[:, end]   # null-space vector (right singular vector for smallest σ)
        tangent ./= norm(tangent)
        println("  σ_min = $(round(F.S[end], sigdigits=4))")

        # Orient direction
        if !isnothing(direction_ref)
            dot(tangent, direction_ref) < 0 && (tangent .*= -1)
        else
            tangent[end] < 0 && (tangent .*= -1)   # default: λ increasing
        end
        return tangent
    end

    # ── Initialization ──
    u_curr = copy(u₂)
    λ_curr = λ₂

    println("Computing initial tangent at λ=$λ_curr...")
    tangent = compute_tangent(u_curr, λ_curr)
    println("Initial tangent: λ̇ = $(round(tangent[end], digits=6))")

    if isnothing(Δs)
        sec_len = norm(vcat(u₂[alg_idx] - u₁[alg_idx], λ₂ - λ₁))
        Δs = clamp(sec_len, Δs_min, Δs_max)
        println("Auto Δs = $Δs")
    end

    path_λ = [λ₁, λ₂]
    path_u = [copy(u₁), copy(u₂)]

    for step in 1:max_steps
        ẋ = tangent[1:n_alg]
        λ̇ = tangent[end]

        # ── Predictor ──
        u_c = copy(u_curr)
        u_c[alg_idx] .+= Δs * ẋ
        λ_c = λ_curr + Δs * λ̇

        # ── Corrector (Newton on augmented system) ──
        converged = false
        iters = 0
        for k in 1:max_newton
            g = eval_g(u_c, λ_c)
            N_val = dot(ẋ, u_c[alg_idx] - u_curr[alg_idx]) + λ̇ * (λ_c - λ_curr) - Δs

            Jfull = eval_full_jac(u_c, λ_c)    # n_alg × (n_alg+1)

            # Augmented Jacobian: append arclength constraint row
            J_aug = zeros(n_alg + 1, n_alg + 1)
            J_aug[1:n_alg, :] = Jfull
            J_aug[end, 1:n_alg] = ẋ
            J_aug[end, end] = λ̇

            rhs = vcat(g, N_val)
            Δ_corr = J_aug \ (-rhs)
            u_c[alg_idx] .+= Δ_corr[1:n_alg]
            λ_c += Δ_corr[end]
            iters = k

            if norm(Δ_corr) < tol * (1 + norm(u_c[alg_idx]))
                converged = true
                break
            end
        end

        # ── Step quality checks ──
        if !converged
            Δs = max(Δs * 0.5, Δs_min)
            println("  Step $step: corrector failed ($(iters) iters), Δs → $(round(Δs, sigdigits=3))")
            Δs ≤ 2*Δs_min && (println("Δs too small, stopping."); break)
            continue
        end

        # Reject if corrector drifted too far from predictor (branch jumping)
        # Absolute floor (1e-2) prevents death spiral when Δs is tiny
        pred_alg = u_curr[alg_idx] .+ Δs * ẋ
        pred_λ = λ_curr + Δs * λ̇
        drift = norm(vcat(u_c[alg_idx] - pred_alg, λ_c - pred_λ))
        drift_tol = max(3 * Δs, 1e-2)
        if drift > drift_tol
            Δs = max(Δs * 0.5, Δs_min)
            println("  Step $step: corrector drift $(round(drift, sigdigits=3)) > $(round(drift_tol, sigdigits=3)), Δs → $(round(Δs, sigdigits=3))")
            Δs ≤ 2*Δs_min && (println("Δs too small, stopping."); break)
            continue
        end

        # ── Accept step ──
        push!(path_λ, λ_c)
        push!(path_u, copy(u_c))
        println("Step $step: λ=$(round(λ_c, digits=6)), λ̇=$(round(λ̇, digits=4)), iters=$iters, Δs=$(round(Δs, sigdigits=3))")

        # True tangent at new point, oriented consistently with previous
        tangent = compute_tangent(u_c, λ_c, tangent)

        # Shift
        u_curr = copy(u_c)
        λ_curr = λ_c

        # Adaptive Δs
        if iters <= 3
            Δs = min(Δs * 1.5, Δs_max)
        elseif iters >= 10
            Δs = max(Δs * 0.5, Δs_min)
        end

        if λ_c >= λ_target
            println("Reached λ=$(round(λ_c, digits=6)) ≥ target $λ_target")
            return path_λ, path_u, true
        end
    end

    return path_λ, path_u, false
end

end # module
