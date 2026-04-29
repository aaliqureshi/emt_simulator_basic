using ForwardDiff, LinearAlgebra


"""Residual g = du[alg_idx] from fx!(du, u, p)."""
function _eval_g(u, p, n, fx!, alg_idx)
    du = zeros(n)
    fx!(du, u, p)
    return du[alg_idx]
end

"""Residual g = du[alg_idx] and its Jacobian w.r.t. u[alg_idx], using fx!(du, u, p)."""
function _eval_g_jac(u, p, n, fx!, alg_idx)
    du = zeros(n)
    fx!(du, u, p)
    g = du[alg_idx]
    y = u[alg_idx]

    if y isa Number
        jac = ForwardDiff.derivative(y_val -> begin
            T = eltype(y_val)
            u_tmp = Vector{T}(u)
            u_tmp[alg_idx] = y_val
            du_tmp = similar(u_tmp)
            fx!(du_tmp, u_tmp, p)
            du_tmp[alg_idx]
        end, y)
    else
        jac = ForwardDiff.jacobian(y_vec -> begin
            T = eltype(y_vec)
            u_tmp = Vector{T}(u)
            u_tmp[alg_idx] .= y_vec
            du_tmp = similar(u_tmp)
            fx!(du_tmp, u_tmp, p)
            du_tmp[alg_idx]
        end, y)
    end

    return g, jac
end

# ── 1. Standard Newton-Raphson ──

function solve_newton!(u, p, fx!, n, alg_idx; tol=1e-6, max_iter=100)
    # n, alg_idx = _setup_alg(u, address)
    hist = Float64[]
    re_eval_count = 0
    jac_updates = 0
    # _, J = _eval_g_jac(u, p, n, alg_idx)
    jac_updates += 1
    correction_norm = Float64[]
    for k in 1:max_iter
        # g, _ = _eval_g_jac(u, p, n, alg_idx)
        g, J = _eval_g_jac(u, p, n, fx!, alg_idx)
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
    push!(hist, norm(_eval_g(u, p, n, fx!, alg_idx), Inf))
    return (converged=false, iters=max_iter, residuals=hist, jac_updates=jac_updates, correction_norm=correction_norm)
end



# ── 5. Natural parameter homotopy (proposed method) ──

function solve_homotopy!(u, p_base, fx!, n, alg_idx; tol=1e-6, max_iter=100, Δλ=0.01, λ_target=1.0)
    hist = Float64[]
    total_iters = 0
    jac_updates = 0
    correction_norm = Float64[]

    p = (p_base..., 0.0)

    # _, J = _eval_g_jac(u, p, n, alg_idx)
    # jac_updates += 1

    for λ in 0.0:Δλ:λ_target
    # for λ in [0.0, 0.6, 0.69, 0.7]
        p = (p_base..., λ)
        # _, J = _eval_g_jac(u, p, n, alg_idx)
        # jac_updates += 1

        # if λ == λ_target
        #     tol = 1e-6
        # else
        #     tol = 1e-3
        # end

        for k in 1:max_iter
            # g, _ = _eval_g_jac(u, p, n, alg_idx)
            g, J = _eval_g_jac(u, p, n, fx!, alg_idx)
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


function fx!(du,u,p)
    x, y = u
    r = p[1]
    du[1] = -x + cos(y)
    du[2] = x^2 + y^2 - r
end

r = 10

u = Float64[4, -4]
p = (r,)
n = 2
alg_idx = 2
solve_newton!(u, p, fx!, n, alg_idx; tol=1e-6, max_iter=100)

using Plots

q(t) = 1e10^(1-t) * 1e-2^(t)
p(t) = 10^(10 - 12*t)

t = collect(0:0.1:1.0)
f1 = @. p(t)
f2 = @. q(t)

plot(t[3:end], f1[3:end])
plot!(t[3:end], f2[3:end])

plot(t, log10.(f1))
plot!(t, log10.(f2))