using LinearAlgebra: norm, SingularException, ZeroPivotException

# Analytic derivative of S = v .* conj(Y*v), restricted to the same
# [Q_PQ; P_non_slack] rows and [V_PQ; theta_non_slack] columns as powerflow!.
# Sparse throughout: a dense 3012wp Jacobian has over 32 million entries.
function _powerflow_jacobian(u, cache::PowerFlowCache)
    V = copy(cache.models.bus.v)
    theta = copy(cache.models.bus.theta)
    pq, ns = cache.v_update_buses, cache.non_slack_buses
    V[pq] = u[cache.balance_reactive_power]
    theta[ns] = u[cache.balance_real_power]
    v = V .* cis.(theta)
    current = cache.Y * v
    Dv = spdiagm(0 => v)
    Di = spdiagm(0 => current)
    Dunit = spdiagm(0 => cis.(theta))
    dS_dV = Dv * conj.(cache.Y * Dunit) + conj.(Di) * Dunit
    dS_dtheta = 1im * Dv * conj.(Di - cache.Y * Dv)
    return [imag.(dS_dV[pq, pq]) imag.(dS_dtheta[pq, ns]);
            real.(dS_dV[ns, pq]) real.(dS_dtheta[ns, ns])]
end

_positive_pq(u, cache) = all(>(0), @view u[cache.balance_reactive_power])
_singular_pf_error(e) = e isa SingularException || e isa ZeroPivotException

function _correct_power_flow(u0, cache, offset; abstol, maxiters)
    u = copy(u0)
    r, trial_r = similar(u), similar(u)
    for iteration in 0:maxiters
        powerflow!(r, u, cache)
        r .-= offset
        residual_norm = norm(r, Inf)
        if isfinite(residual_norm) && residual_norm <= abstol && _positive_pq(u, cache)
            return (; u, converged=true, iterations=iteration, residual_norm)
        end
        if !isfinite(residual_norm) || iteration == maxiters
            return (; u, converged=false, iterations=iteration, residual_norm)
        end
        direction = try
            -(_powerflow_jacobian(u, cache) \ r)
        catch e
            _singular_pf_error(e) || rethrow()
            return (; u, converged=false, iterations=iteration, residual_norm)
        end
        all(isfinite, direction) ||
            return (; u, converged=false, iterations=iteration, residual_norm)

        # Armijo backtracking on the infinity norm, rejecting negative PQ voltages.
        alpha = 1.0
        accepted = false
        for _ in 0:15
            trial = u + alpha * direction
            if _positive_pq(trial, cache)
                powerflow!(trial_r, trial, cache)
                trial_r .-= offset
                if norm(trial_r, Inf) <= (1 - 1e-4 * alpha) * residual_norm
                    u = trial
                    accepted = true
                    break
                end
            end
            alpha /= 2
        end
        accepted || return (; u, converged=false, iterations=iteration, residual_norm)
    end
end

"""
    solve_power_flow_continuation!(sys; initial_step=0.05, min_step=1e-6,
        max_step=0.2, abstol=1e-9, maxiters=20, maxsteps=200, verbose=true)

Solve the existing power-flow equations from a flat start using
`H(u, lambda) = F(u) - (1-lambda)*F(u_flat)`. PQ magnitudes start at 1 p.u.
and non-slack angles at zero; prescribed PV/slack magnitudes and the slack
angle stay fixed. No PQ-voltage or non-slack-angle guess from the file is used.
The flat state solves H at lambda=0 exactly, including taps and shunts.

Trace to lambda=1 using a tangent predictor, sparse damped Newton corrector,
and adaptive steps. Failed corrections are discarded and the step is halved.
`maxiters` limits each corrector; `maxsteps` limits all attempts, including
rejections. This is natural-parameter continuation, not a fold-tracing solver.

Returns a named tuple with `u`, `converged`, `retcode` (a Symbol), `lambda`,
`residual_norm` (the ORIGINAL F), `iterations` (all corrector iterations),
and `history` (accepted and rejected attempts). Only a successful lambda=1
solution updates bus voltages, angles, and d-q components. Device limits and
bus classifications are unchanged from the existing power-flow formulation.
"""
function solve_power_flow_continuation!(sys; initial_step=0.05, min_step=1e-6,
        max_step=0.2, abstol=1e-9, maxiters=20, maxsteps=200, verbose=true)
    all(isfinite, (min_step, initial_step, max_step)) &&
        0 < min_step <= initial_step <= max_step <= 1 ||
        throw(ArgumentError("require 0 < min_step <= initial_step <= max_step <= 1"))
    isfinite(abstol) && abstol > 0 || throw(ArgumentError("abstol must be finite and positive"))
    maxiters isa Integer && maxiters >= 0 || throw(ArgumentError("maxiters must be a nonnegative integer"))
    maxsteps isa Integer && maxsteps > 0 || throw(ArgumentError("maxsteps must be a positive integer"))

    u = vcat(ones(length(sys.v_update_buses)), zeros(length(sys.non_slack_buses)))
    cache = PowerFlowCache(sys, u)
    f_flat = similar(u)
    powerflow!(f_flat, u, cache)
    all(isfinite, f_flat) || throw(ArgumentError("nonfinite residual at flat start"))
    lambda = 0.0
    step = Float64(initial_step)
    history = NamedTuple[]
    iterations = 0
    retcode = :MaxSteps

    for _ in 1:maxsteps
        target = min(1.0, lambda + step)
        actual_step = target - lambda
        # H_u * du/dlambda + H_lambda = 0, with H_lambda = F(u_flat).
        tangent = try
            _powerflow_jacobian(u, cache) \ (-f_flat)
        catch e
            _singular_pf_error(e) || rethrow()
            retcode = :SingularJacobian
            break
        end
        if !all(isfinite, tangent)
            retcode = :SingularJacobian
            break
        end
        predicted = u + actual_step * tangent
        if !all(isfinite, predicted) || !_positive_pq(predicted, cache)
            predicted = copy(u)
        end
        correction = _correct_power_flow(predicted, cache, (1 - target) * f_flat;
                                         abstol, maxiters)
        iterations += correction.iterations
        push!(history, (; lambda=target, step=actual_step, accepted=correction.converged,
                         iterations=correction.iterations,
                         residual_norm=correction.residual_norm))
        if verbose
            println("PF continuation: lambda=$target, accepted=$(correction.converged), " *
                    "iterations=$(correction.iterations), |H|inf=$(correction.residual_norm)")
        end
        if correction.converged
            u = correction.u
            lambda = target
            if lambda == 1.0
                retcode = :Success
                break
            end
            if correction.iterations <= 4
                step = min(max_step, step * 1.5)
            elseif correction.iterations >= 8
                step = max(min_step, step / 2)
            end
        else
            step = actual_step / 2
            if step < min_step
                retcode = :StepSizeTooSmall
                break
            end
        end
    end

    residual = similar(u)
    powerflow!(residual, u, cache)
    residual_norm = norm(residual, Inf)
    converged = lambda == 1.0 && isfinite(residual_norm) && residual_norm <= abstol
    if converged
        @views sys.models.bus.v[sys.v_update_buses] .= u[cache.balance_reactive_power]
        @views sys.models.bus.theta[sys.non_slack_buses] .= u[cache.balance_real_power]
        phasor2DP!(sys.models.bus)
    elseif retcode == :Success
        retcode = :ResidualFailure
    end
    return (; u, converged, retcode, lambda, residual_norm, iterations, history)
end
