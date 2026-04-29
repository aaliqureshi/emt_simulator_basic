using Pkg; Pkg.activate(".")
using LinearAlgebra, ForwardDiff


function _eval_jac(u, p, f!)
    jac = ForwardDiff.jacobian(y -> begin
        du = similar(y)
        f!(du, y, p)
        du
    end, u)

    return jac
end


function solve_newton!(u, p, f!; tol=1e-6, max_iter=1000)
    hist = Float64[]
    g = similar(u)

    for k in 1:max_iter
        g = similar(u)
        f!(g, u, p)
        J = _eval_jac(u, p, f!)
        push!(hist, norm(g, Inf))
        any(!isfinite, g) && return (converged=false, iters=k, residuals=hist)

        Δy = J \ (-g)
        any(!isfinite, Δy) && return (converged=false, iters=k, residuals=hist)
        u .+= Δy

        norm(Δy) < tol * (1 + norm(u)) &&
            return (converged=true, iters=k, residuals=hist)
    end
    push!(hist, norm(f!(g, u, p), Inf))
    return (converged=false, iters=max_iter, residuals=hist)
end

function solve_homotopy!(u, p_base, f!; tol=1e-6, max_iter=1000, Δλ=0.001, λ_target=1.0)
    hist = Float64[]
    total_iters = 0

    for λ in 0.0:Δλ:λ_target
        p = (p_base..., λ)

        for k in 1:max_iter
            g = similar(u)
            f!(g, u, (p_base..., λ))
            push!(hist, norm(g, Inf))
            J = _eval_jac(u, (p_base..., λ), f!)

            Δy = J \ (-g)
            u .+= Δy
            total_iters += 1

            if norm(Δy) < tol * (1 + norm(u))
                break
            end

            if k == max_iter
                return (converged=false, iters=total_iters, residuals=hist, λ_failed=λ)
            end
        end
    end

    return (converged=true, iters=total_iters, residuals=hist, λ_failed=nothing)
end


function simple_cct!(du, u, p)
    i = u[1]

    Vs, R = p

    du[1] = Vs - R * i
end

function r_cct_sat!(du, u, p)
    i = u[1]

    Vs, R = p

    if i > 10.0
        R = 1.2
    end

    du[1] = Vs - R * i
end

function r_cct_sat_homotopy!(du, u, p)
    i = u[1]

    Vs, R, λ = p

    if i > 10.0
        R = 1.2
    end

    @show Vs

    du[1] = Vs - R * i
end


u = [0.0]
Vs = 10.0
R = 0.5
p = (Vs, R)
# solve_newton!(u, p, simple_cct!)
r = solve_newton!(u, p, r_cct_sat!)

λ_target = 1.0
λ_step = 0.01
λ_history = [0.0]
r_homotopy = solve_homotopy!(u, p, r_cct_sat_homotopy!, λ_target=λ_target, Δλ=λ_step)


using Plots
plot(r.residuals)