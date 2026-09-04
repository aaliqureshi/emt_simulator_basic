using ForwardDiff, LinearAlgebra

# ── Shared helpers ──

function _setup_alg(u, address)
    n = length(u)
    n_mech = length(address["delta"]) + length(address["omega"])
    # n_mech = length(address["delta"]) + length(address["omega"]) + length(address["line_id"]) + length(address["line_iq"])
    alg_idx = n_mech+1:n
    return n, alg_idx
end

function _eval_g(u, p, n, alg_idx)
    du = zeros(n)
    solve_dynamic_sim!(du, u, p, 30.0)
    return du[alg_idx]
end

function _eval_g_jac(u, p, n, alg_idx)
    du = zeros(n)
    solve_dynamic_sim!(du, u, p, 30.0)
    g = du[alg_idx]

    jac = ForwardDiff.jacobian(y -> begin
        T = eltype(y)
        u_tmp = Vector{T}(u)
        u_tmp[alg_idx] .= y
        du_tmp = similar(u_tmp)
        solve_dynamic_sim!(du_tmp, u_tmp, p, 30.0)
        du_tmp[alg_idx]
    end, u[alg_idx])

    return g, jac
end

# ── 1. Standard Newton-Raphson ──

function solve_newton!(u, p, address; tol=1e-6, max_iter=100, always_new=false)
    n, alg_idx = _setup_alg(u, address)
    hist = Float64[]
    re_eval_count = 0
    jac_updates = 0
    # if !always_new
    #     _, J = _eval_g_jac(u, p, n, alg_idx)
    # end
    # jac_updates += 1
    correction_norm = Float64[]
    local J
    for k in 1:max_iter
        if always_new || k == 1
            g, J = _eval_g_jac(u, p, n, alg_idx)
            jac_updates += 1
        else
            g, _ = _eval_g_jac(u, p, n, alg_idx)
        end
        # g, J = _eval_g_jac(u, p, n, alg_idx)
        push!(hist, norm(g, Inf))
        any(!isfinite, g) && return (converged=false, iters=k, residuals=hist, jac_updates=jac_updates)

        # if k == 30 || k == 60
        #     _, J = _eval_g_jac(u, p, n, alg_idx)
        #     jac_updates += 1
        # end

        Δy = J \ (-g)
        push!(correction_norm, norm(Δy))
        any(!isfinite, Δy) && return (converged=false, iters=k, residuals=hist, jac_updates=jac_updates, correction_norm=correction_norm)
        u[alg_idx] .+= Δy

        norm(Δy) < tol * (1 + norm(u[alg_idx])) &&
            return (converged=true, iters=k, residuals=hist, jac_updates=jac_updates, correction_norm=correction_norm)
    end
    push!(hist, norm(_eval_g(u, p, n, alg_idx), Inf))
    return (converged=false, iters=max_iter, residuals=hist, jac_updates=jac_updates, correction_norm=correction_norm)
end

# ── 2. Damped Newton (fixed α) ──

function solve_damped_newton!(u, p, address; α=0.5, tol=1e-6, max_iter=100)
    n, alg_idx = _setup_alg(u, address)
    hist = Float64[]

    for k in 1:max_iter
        g, J = _eval_g_jac(u, p, n, alg_idx)
        push!(hist, norm(g, Inf))
        any(!isfinite, g) && return (converged=false, iters=k, residuals=hist)

        Δy = J \ (-g)
        any(!isfinite, Δy) && return (converged=false, iters=k, residuals=hist)
        u[alg_idx] .+= α * Δy

        norm(α * Δy) < tol * (1 + norm(u[alg_idx])) &&
            return (converged=true, iters=k, residuals=hist)
    end
    push!(hist, norm(_eval_g(u, p, n, alg_idx), Inf))
    return (converged=false, iters=max_iter, residuals=hist)
end

# ── 3. Newton with Armijo backtracking line search ──

function solve_backtracking_newton!(u, p, address; tol=1e-6, max_iter=100, c₁=1e-4, ρ=0.5)
    n, alg_idx = _setup_alg(u, address)
    hist = Float64[]

    for k in 1:max_iter
        g, J = _eval_g_jac(u, p, n, alg_idx)
        push!(hist, norm(g, Inf))
        any(!isfinite, g) && return (converged=false, iters=k, residuals=hist)

        Δy = J \ (-g)
        any(!isfinite, Δy) && return (converged=false, iters=k, residuals=hist)

        # Merit function: φ(α) = ½||g(u + α·Δy)||²
        φ₀ = 0.5 * norm(g)^2
        dφ₀ = -norm(g)^2   # directional derivative: gᵀ·J·Δy = gᵀ·(-g)

        α = 1.0
        for _ in 1:30
            u_trial = copy(u)
            u_trial[alg_idx] .+= α * Δy
            g_trial = _eval_g(u_trial, p, n, alg_idx)
            any(!isfinite, g_trial) && (α *= ρ; continue)

            0.5 * norm(g_trial)^2 ≤ φ₀ + c₁ * α * dφ₀ && break
            α *= ρ
            α < 1e-14 && break
        end

        u[alg_idx] .+= α * Δy

        norm(α * Δy) < tol * (1 + norm(u[alg_idx])) &&
            return (converged=true, iters=k, residuals=hist)
    end
    push!(hist, norm(_eval_g(u, p, n, alg_idx), Inf))
    return (converged=false, iters=max_iter, residuals=hist)
end

# ── 4. Levenberg-Marquardt (adaptive damping) ──

function solve_levenberg_marquardt!(u, p, address; tol=1e-6, max_iter=100, μ₀=1e-3)
    n, alg_idx = _setup_alg(u, address)
    hist = Float64[]
    μ = μ₀

    for k in 1:max_iter
        g, J = _eval_g_jac(u, p, n, alg_idx)
        push!(hist, norm(g, Inf))
        any(!isfinite, g) && return (converged=false, iters=k, residuals=hist)

        # (JᵀJ + μI)Δy = -Jᵀg
        Δy = (J' * J + μ * I) \ (-J' * g)
        any(!isfinite, Δy) && return (converged=false, iters=k, residuals=hist)

        # Evaluate trial point
        u_trial = copy(u)
        u_trial[alg_idx] .+= Δy
        g_trial = _eval_g(u_trial, p, n, alg_idx)

        # Gain ratio: actual vs predicted reduction
        actual    = norm(g)^2 - norm(g_trial)^2
        predicted = norm(g)^2 - norm(g + J * Δy)^2
        ρ_gain = predicted > 0 ? actual / predicted : -1.0

        if ρ_gain > 0.25       # good step → accept, decrease μ
            u[alg_idx] .+= Δy
            μ = max(μ / 3, 1e-15)
        elseif ρ_gain > 0      # acceptable → accept, keep μ
            u[alg_idx] .+= Δy
        else                   # bad step → reject, increase μ
            μ = min(μ * 2, 1e10)
        end

        if ρ_gain > 0 && norm(Δy) < tol * (1 + norm(u[alg_idx]))
            return (converged=true, iters=k, residuals=hist)
        end
    end
    push!(hist, norm(_eval_g(u, p, n, alg_idx), Inf))
    return (converged=false, iters=max_iter, residuals=hist)
end

# ── 5. Natural parameter homotopy (proposed method) ──

function solve_homotopy!(u, p_base, address; 
                         tol=1e-6, 
                         max_iter=100, 
                         Δλ=0.01, 
                         λ_target=1.0,
                         always_new=false)
    n, alg_idx = _setup_alg(u, address)
    hist = Float64[]
    total_iters = 0
    jac_updates = 0
    correction_norm = Float64[]
    # p = (p_base..., 0.0)

    for λ in 0.0:Δλ:λ_target
        p = (p_base..., λ)
        local J
        for k in 1:max_iter
            # local J
            if always_new || k == 1
                g, J = _eval_g_jac(u, p, n, alg_idx)
                jac_updates += 1
            else
                g, _ = _eval_g_jac(u, p, n, alg_idx)
            end
            # g, _ = _eval_g_jac(u, p, n, alg_idx)
            # g, J = _eval_g_jac(u, p, n, alg_idx)
            push!(hist, norm(g, Inf))

            # if k == 30 || k == 60
            #     _, J = _eval_g_jac(u, p, n, alg_idx)
            #     jac_updates += 1
            # end

            Δy = J \ (-g)
            u[alg_idx] .+= Δy
            total_iters += 1

            push!(correction_norm, norm(Δy, 2))

            if norm(Δy) < tol * (1 + norm(u[alg_idx]))
                break
            end

            if k == max_iter
                return (converged=false, iters=total_iters, residuals=hist, correction_norm=correction_norm, λ_failed=λ, jac_updates=jac_updates)
            end
        end
        println("converged in $(total_iters) iterations")
    end

    return (converged=true, iters=total_iters, residuals=hist, correction_norm=correction_norm, λ_failed=nothing, jac_updates=jac_updates)
end

# function solve_homotopy2!(u, p_base, address; tol=1e-6, max_iter=1000, Δλ=0.001, λ_target=1.0)
#     n, alg_idx = _setup_alg(u, address)
#     hist = Float64[]
#     total_iters = 0

#     for λ in 0.0:Δλ:λ_target
#         p = (p_base..., λ)

#         for k in 1:max_iter
#             r = solve_levenberg_marquardt!(u, p_base, address)

#         end
#     end

#     return (converged=true, iters=total_iters, residuals=hist, λ_failed=nothing)
# end

# ── 6. Adaptive LM-Homotopy hybrid (proposed method) ──

function solve_homotopy_lm!(u, p_base, address;
                            tol=1e-6, max_iter=100, μ₀=1e-3,
                            Δλ_init=0.01, Δλ_min=1e-6, Δλ_max=0.2, λ_target=1.0)
    n, alg_idx = _setup_alg(u, address)
    hist = Float64[]
    λ_hist = Float64[]
    μ = μ₀
    total_iters = 0
    λ = 0.0
    Δλ = Δλ_init

    p = (p_base..., 0.0)
    # _, J = _eval_g_jac(u, p, n, alg_idx)

    while λ < λ_target
        Δλ = min(Δλ, λ_target - λ)   # don't overshoot target
        λ_next = λ + Δλ
        p = (p_base..., λ_next)

        # Save state for rollback on failure
        u_saved = copy(u)
        μ_saved = μ

        step_iters = 0
        step_converged = false

        for k in 1:max_iter
            g, J = _eval_g_jac(u, p, n, alg_idx)
            # g, _ = _eval_g_jac(u, p, n, alg_idx)
            total_iters += 1
            step_iters = k
            push!(hist, norm(g, Inf))
            any(!isfinite, g) && break

            # LM step: (JᵀJ + μI)Δy = -Jᵀg
            Δy = (J' * J + μ * I) \ (-J' * g)
            any(!isfinite, Δy) && break

            # Trial evaluation
            u_trial = copy(u)
            u_trial[alg_idx] .+= Δy
            g_trial = _eval_g(u_trial, p, n, alg_idx)

            # Gain ratio
            actual    = norm(g)^2 - norm(g_trial)^2
            predicted = norm(g)^2 - norm(g + J * Δy)^2
            ρ_gain = predicted > 0 ? actual / predicted : -1.0

            if ρ_gain > 0.25
                u[alg_idx] .+= Δy
                μ = max(μ / 3, 1e-15)
            elseif ρ_gain > 0
                u[alg_idx] .+= Δy
            else
                μ = min(μ * 2, 1e10)
            end

            if ρ_gain > 0 && norm(Δy) < tol * (1 + norm(u[alg_idx]))
                step_converged = true
                break
            end
        end

        if step_converged
            # Accept step, advance λ
            λ = λ_next
            push!(λ_hist, λ)

            # Adaptive Δλ: grow if easy, keep if moderate
            if step_iters <= 3
                Δλ = min(Δλ * 2.0, Δλ_max)
            elseif step_iters <= 6
                Δλ = min(Δλ * 1.2, Δλ_max)
            end
            # keep μ warm-started for next step
        else
            # Reject step: rollback state and shrink Δλ
            u[alg_idx] .= u_saved[alg_idx]
            μ = μ_saved
            Δλ = Δλ * 0.5

            if Δλ < Δλ_min
                println("Δλ too small at λ=$(round(λ, digits=6)), stopping.")
                return (converged=false, iters=total_iters, residuals=hist,
                        λ_hist=λ_hist, λ_failed=λ + Δλ)
            end
        end
    end

    return (converged=true, iters=total_iters, residuals=hist,
            λ_hist=λ_hist, λ_failed=nothing)
end


# ── 7. Predictor-Corrector Homotopy with PI Step Size Control ──
#
# Predictor:  secant extrapolation from two previous solutions
# Corrector:  Newton-Raphson using predicted value as initial guess
# Step size:  PI controller on log(Δλ), driven by Newton difficulty
#
# Newton difficulty metric:
#   d_k = n_iters / n_target
# where n_target is the desired number of corrector iterations.
#
# PI controller (Söderlind digital form):
#   log(Δλ_{k+1}/Δλ_k) = k_p * (e_k - e_{k-1}) + k_i * e_k
# where e_k = -log(d_k) (positive when step is easy, negative when hard).

function _newton_corrector!(u, p, n, alg_idx; tol=1e-6, max_iter=50, always_new=false)
    converged = false
    n_iters = 0
    first_residual = NaN
    contraction = NaN
    prev_corr_norm = NaN
    jac_updates = 0

    local J
    for k in 1:max_iter
        if always_new || k == 1
            g, J = _eval_g_jac(u, p, n, alg_idx)
            jac_updates += 1
        else
            g, _ = _eval_g_jac(u, p, n, alg_idx)
        end
        # g, J = _eval_g_jac(u, p, n, alg_idx)
        res_norm = norm(g, Inf)

        if k == 1
            first_residual = res_norm
        end
        any(!isfinite, g) && break

        Δy = J \ (-g)
        any(!isfinite, Δy) && break

        corr_norm = norm(Δy)
        if k >= 2 && isfinite(prev_corr_norm) && prev_corr_norm > 0
            contraction = corr_norm / prev_corr_norm
        end
        prev_corr_norm = corr_norm

        u[alg_idx] .+= Δy
        n_iters = k

        if corr_norm < tol * (1 + norm(u[alg_idx]))
            converged = true
            break
        end
    end

    return (converged=converged, iters=n_iters,
            first_residual=first_residual, contraction=contraction)
end

function solve_adaptive_homotopy!(u, p_base, address;
                                   tol=1e-6, max_iter=50, max_steps=500,
                                   Δλ_init=0.05, Δλ_min=1e-10, Δλ_max=0.5,
                                   λ_target=1.0, n_target=4,
                                   k_p=0.7, k_i=0.4,
                                   always_new=false)

    n, alg_idx = _setup_alg(u, address)
    n_alg = length(alg_idx)

    # History arrays
    λ_hist    = Float64[]
    iter_hist = Int[]
    Δλ_hist   = Float64[]
    total_newton_iters = 0
    total_jac_evals    = 0

    # ── Bootstrap point 1: solve at λ = 0 ──
    λ_prev = 0.0
    p = (p_base..., λ_prev)
    r = _newton_corrector!(u, p, n, alg_idx; tol=tol, max_iter=max_iter, always_new=always_new)
    total_newton_iters += r.iters
    total_jac_evals    += r.iters
    if !r.converged
        println("Homotopy bootstrap failed at λ=0")
        return (converged=false, λ_hist=λ_hist, iter_hist=iter_hist,
                Δλ_hist=Δλ_hist, total_newton_iters=total_newton_iters,
                total_jac_evals=total_jac_evals, λ_failed=0.0)
    end
    if λ_target == 0.0
        push!(λ_hist, λ_prev)
        return (converged=true, λ_hist=λ_hist, iter_hist=iter_hist,
        Δλ_hist=Δλ_hist, total_newton_iters=total_newton_iters,
        total_jac_evals=total_jac_evals, λ_failed=nothing)
    end
    y_prev = copy(u[alg_idx])
    push!(λ_hist, λ_prev)

    # ── Bootstrap point 2: first small step ──
    Δλ = Δλ_init
    λ_curr = Δλ
    p = (p_base..., λ_curr)
    r = _newton_corrector!(u, p, n, alg_idx; tol=tol, max_iter=max_iter, always_new=always_new)
    total_newton_iters += r.iters
    total_jac_evals    += r.iters
    if !r.converged
        println("Homotopy bootstrap failed at λ=$λ_curr")
        return (converged=false, λ_hist=λ_hist, iter_hist=iter_hist,
                Δλ_hist=Δλ_hist, total_newton_iters=total_newton_iters,
                total_jac_evals=total_jac_evals, λ_failed=λ_curr)
    end
    y_curr = copy(u[alg_idx])
    push!(λ_hist, λ_curr)
    push!(iter_hist, r.iters)
    push!(Δλ_hist, Δλ)

    # PI controller state
    e_prev = 0.0

    # ── Main continuation loop ──
    for step in 1:max_steps
        # Don't overshoot target
        Δλ_step = min(Δλ, λ_target - λ_curr)
        Δλ_step <= 0 && break
        λ_next = λ_curr + Δλ_step

        # ── Predictor: secant extrapolation ──
        dλ = λ_curr - λ_prev
        if dλ > 0
            slope = (y_curr - y_prev) / dλ
            y_pred = y_curr + slope * Δλ_step
        else
            y_pred = copy(y_curr)
        end

        # Save state for rollback
        u_saved = copy(u[alg_idx])
        u[alg_idx] .= y_pred

        # ── Corrector: Newton-Raphson from predicted initial guess ──
        p = (p_base..., λ_next)
        r = _newton_corrector!(u, p, n, alg_idx; tol=tol, max_iter=max_iter, always_new=always_new)
        total_newton_iters += r.iters
        total_jac_evals    += r.iters

        if !r.converged
            # Reject step: rollback and shrink Δλ aggressively
            u[alg_idx] .= u_saved
            Δλ = Δλ * 0.25
            if Δλ < Δλ_min
                println("Homotopy failed: Δλ below minimum at λ=$(round(λ_curr, digits=8))")
                return (converged=false, λ_hist=λ_hist, iter_hist=iter_hist,
                        Δλ_hist=Δλ_hist, total_newton_iters=total_newton_iters,
                        total_jac_evals=total_jac_evals, λ_failed=λ_next)
            end
            println("  Step $step: corrector failed at λ=$(round(λ_next, digits=6)), shrinking Δλ → $(round(Δλ, sigdigits=3))")
            continue
        end

        # ── Accept step ──
        y_prev = copy(y_curr)
        λ_prev = λ_curr
        y_curr = copy(u[alg_idx])
        λ_curr = λ_next

        push!(λ_hist, λ_curr)
        push!(iter_hist, r.iters)
        push!(Δλ_hist, Δλ_step)

        # ── PI step size control ──
        # Newton difficulty: ratio of iterations taken to target iterations
        difficulty = r.iters / n_target
        # Clamp difficulty away from zero to avoid log(0)
        difficulty = max(difficulty, 0.1)

        # Error signal in log-space (positive = easy, negative = hard)
        e_k = -log(difficulty)

        # PI update on log(Δλ)
        Δ_log = k_p * (e_k - e_prev) + k_i * e_k
        Δλ = Δλ * exp(Δ_log)
        Δλ = clamp(Δλ, Δλ_min, Δλ_max)

        e_prev = e_k

        println("Step $step: λ=$(round(λ_curr, digits=6)), iters=$(r.iters), " *
                "difficulty=$(round(difficulty, digits=2)), " *
                "contraction=$(isnan(r.contraction) ? "N/A" : string(round(r.contraction, digits=3))), " *
                "Δλ_next=$(round(Δλ, sigdigits=3))")

        if λ_curr >= λ_target
            println("Homotopy converged to λ=$λ_target in $step steps, " *
                    "$(total_newton_iters) total Newton iterations")
            return (converged=true, λ_hist=λ_hist, iter_hist=iter_hist,
                    Δλ_hist=Δλ_hist, total_newton_iters=total_newton_iters,
                    total_jac_evals=total_jac_evals, λ_failed=nothing)
        end
    end

    println("Homotopy: max steps ($max_steps) reached at λ=$(round(λ_curr, digits=6))")
    return (converged=false, λ_hist=λ_hist, iter_hist=iter_hist,
            Δλ_hist=Δλ_hist, total_newton_iters=total_newton_iters,
            total_jac_evals=total_jac_evals, λ_failed=λ_curr)
end