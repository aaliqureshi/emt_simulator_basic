using ForwardDiff
using LinearAlgebra

# =====================================================================
# Constraint-stabilization reset for inconsistent DAE initial conditions
# (White et al., "Stabilized Neural Differential Equations ...", 2023).
#
# After an event the network algebraic equations change, so the carried-
# over state (x^+, y^-) is OFF the new manifold:  g(x^+, y^-) != 0.
# Physically the differential states x are continuous, so we FREEZE x = x^+
# and drive only the algebraic variables y back onto M^+ = {g = 0} in a
# fictitious time tau via the gradient-flow / continuous-Newton dynamics
#
#       dy/dtau = -gamma * g_y(x^+, y)^T * H * g(x^+, y)
#
#   mode = :gradient -> H = I                 (steepest descent on 1/2 g'Hg)
#   mode = :newton   -> dy/dtau = -gamma * g_y \ g   (continuous Newton)
#
# State layout (from build_dynamic_address): the differential states
# (delta, omega, line_id, line_iq; nonzero mass) come first, the algebraic
# states (balance_d/q, fault_id/iq, gen_id/iq; zero mass) come last, so the
# algebraic block is u[nx+1:end]. `alg_mask` makes this robust to reordering.
# =====================================================================

"Boolean mask of algebraic states (zero diagonal in the mass matrix)."
algebraic_mask(mass_matrix) = [mass_matrix[i, i] == 0.0 for i in 1:size(mass_matrix, 1)]

"""
    stabilize_algebraic!(u, p, t, mass_matrix; ...)

Project the algebraic part of state `u` onto the manifold g = 0, holding the
differential states fixed. Mutates `u` in place and returns a NamedTuple with
the residual-norm history. The full residual is `du = solve_dynamic_sim!`,
whose algebraic rows are the constraint g.
"""
function stabilize_algebraic!(u0, p, t, mass_matrix;
                              gamma = 1.0, dtau = 1e-3, tau_max = 50.0,
                              tol = 1e-5, iter_max = 50000, mode = :gradient)

    alg   = findall(algebraic_mask(mass_matrix))   # indices of algebraic vars/eqs
    u = copy(u0)
    n     = length(u)
    du    = zeros(n)
    mu = 1/(gamma * dtau)
    eye = Matrix{Float64}(I, length(alg), length(alg))

    # g(y): set the algebraic entries of u to y, return the algebraic residual.
    # Differentiating this gives g_y directly (ng x ny), with no wasted g_x columns.
    residual_g = function (y)
        uu = eltype(y).(u)        # promote for ForwardDiff duals
        uu[alg] = y
        d = similar(uu)
        solve_dynamic_sim!(d, uu, p, t)
        d[alg]
    end

    g_hist = Float64[]; tau = 0.0; iter = 0
    g = residual_g(u[alg])
    push!(g_hist, norm(g))

    gy = ForwardDiff.jacobian(residual_g, u[alg])
    D = inv(diagm(diag(gy)))
    


    while norm(g) > tol && tau < tau_max && iter < iter_max
        gy = ForwardDiff.jacobian(residual_g, u[alg])      # g_y, square (ng x ny)
        if mode === :newton
            dy = -(gy \ g)                                  # continuous Newton step
        elseif mode === :lm
            dy = -((transpose(gy) * gy + mu*eye)\(transpose(gy) * g))
        elseif mode === :new
            # G = eye*g
            # G = diagm(diag(gy))
            # dy = -(transpose(gy) * inv(G * transpose(G)) * g)
            # D = diagm(diag(gy))
            dy = -g
        else
            dy = -transpose(gy) * g                       # H = I gradient step
        end
        u[alg] .+= dtau .* gamma .* dy
        # u[alg] .+= gamma .* dy
        g = residual_g(u[alg])
        tau += dtau 
        iter += 1
        push!(g_hist, norm(g))
        @show norm(g, 1)
    end

    return (u = u, g_norm = norm(g), iters = iter, tau = tau, history = g_hist,
            converged = norm(g) <= tol)
end

res = stabilize_algebraic!(u1, p, t, mass_matrix; mode = :newton, tol = 1e-10)
@show res.converged, res.iters, res.g_norm