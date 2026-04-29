using Printf, Plots

# --- System Parameters ---
const V_ref  = 1.0  # GFM internal voltage reference (pu)
const R_g    = 0.1  # Output virtual/physical impedance (pu)
const P_load = 0.9  # Constant power demand (pu)
const I_max  = 3.0  # Saturated fault current limit (GFL mode)

# --- Manifold 1: Grid Forming Mode ---
function gfm_residual(V)
    return (V_ref - V) / R_g - (P_load / V)
end

function gfm_jacobian(V)
    return -1.0 / R_g + (P_load / (V^2))
end

function gfl_residual(V)
   return I_max - (P_load / V)
end

function gfl_jacobian(V)
    return (P_load / V^2)
end

function homotopy_residual_mode_switch(V, lambda)
    return (1-lambda) * gfm_residual(V) + lambda * gfl_residual(V)
end

function homotopy_jacobian_mode_switch(V, lambda)
    return (1-lambda)*gfm_jacobian(V) + lambda*gfl_jacobian(V)
end

# --- Newton-Raphson Solver ---
function solve_nr(residual_func, jacobian_func, initial_guess; max_iter=10)
    V = initial_guess
    @printf("Initial Guess: V = %.4f pu\n", V)
    
    for i in 1:max_iter
        res = residual_func(V)
        jac = jacobian_func(V)
        
        @printf("Iter %d: V = %8.4f | Residual = %8.4f | Jacobian = %8.4f\n", i, V, res, jac)
        
        if abs(res) < 1e-6
            println("--> Converged!\n")
            return V
        end
        
        # The Convergence Killer
        if abs(jac) < 1e-12
            println("--> FATAL ERROR: Jacobian is singular (divide by zero)!\n")
            return NaN
        end
        
        # Newton Update Step
        V = V - (res / jac)
    end
    println("--> Failed to converge after max iterations.\n")
    return V
end

# --- Robust NR Solver for Homotopy Steps ---
function solve_nr_homotopy(V_guess, lambda; max_iter=10)
    V = V_guess
    for i in 1:max_iter
        res = homotopy_residual_mode_switch(V, lambda)
        jac = homotopy_jacobian_mode_switch(V, lambda)
        
        if abs(res) < 1e-8
            println("--> $(lambda) ---> Converged!\n")
            return V
        end
        V = V - (res / jac)
    end
    return NaN # Failure
end




####
lw = 2.5
ms = 7.0
V = 0.1:0.01:1.0
# V = -0.0:0.01:3.0
vv = @. gfm_residual(V)
gfm_curve= vv
kk = @. 3.0 - P_load/V
pa = plot(V, vv, lw=lw, label= "GFM Mode")




using Printf

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

println("=== Pre-fault GFM solve ===")
V_prefault = solve_nr(gfm_residual, gfm_jacobian, 1.0)
@printf("Pre-fault solution: %.8f\n\n", V_prefault)

println("=== Direct post-switch Newton on GFL ===")
V_direct = solve_nr(gfl_residual, gfl_jacobian, V_prefault)
@printf("Direct GFL Newton result: %.8f\n\n", V_direct)

println("=== Homotopy continuation to GFL manifold ===")
V_hom = solve_homotopy(V_prefault, nsteps=5)
@printf("Homotopy final solution: %.8f\n", V_hom)

println("\nExact GFL solution should be:")
@printf("V_exact = %.8f\n", P_load / I_max)

vals = [0.9, 0.818, 0.75, 0.69, 0.64, 0.6, 0.56, 0.53, 0.5, 0.47, 0.45, 0.43, 0.41, 0.39, 0.375, 0.36, 0.346, 0.333, 0.321, 0.310, 0.3]

scatter!([vals], [zeros(21)])

l0 = 0.0
v = @. 0.1:0.01:1.0
# v = @. 0.1:0.01:3.0
V0 = V_prefault
r1 = @. homotopy_residual(v, 0.0, V0)

r2 = @. homotopy_residual(v, 0.2, V0)
r3 = @. homotopy_residual(v, 0.4, V0)
r4 = @. homotopy_residual(v, 0.6, V0)
r5 = @. homotopy_residual(v, 0.8, V0)
r6 = @. homotopy_residual(v, 1.0, V0)


using Plots
begin
lw= 2.5
pa = plot(V, vv, lw=lw, label= "GFM Mode", color=:blue)
pb = plot!(v ,r1, label="", lw=lw, ls=:dot, color=:teal, alpha=0.7)
plot!(v ,r2, label="", ls=:dot, lw=lw, color=:teal, alpha=0.7)
plot!(v ,r3, label="", ls=:dot, lw=lw, color=:teal, alpha=0.7)
plot!(v ,r4, label="", ls=:dot, lw=lw, color=:teal, alpha=0.7)
plot!(v ,r5, label="intermediate manifolds", ls=:dot, lw=lw, color=:teal, alpha=0.7)
plot!(v ,r6, label="GFL Mode (λ = 1.0)",lw=lw, color=:orange)
hline!([0.0], lw = 1.5, color=:black, label="")
# scatter!([0.6428, 0.5, 0.409, 0.3462], zeros(4), label="intermediate solutions")
scatter!(v_all[2:end], res_all[2:end], label= "Newton Continuation points", markeralpha = 0.5, color=:red)
plot!(v_all, res_all, lw = 1.5, label="continuation path", color=:sienna)
scatter!([0.9], zeros(1), marker=:utriangle, label="initial guess", markersize=6.0, markercolor=:cyan)
scatter!([0.3], zeros(1), marker=:star5, label="solution", markersize=8.0, markercolor=:cyan)
ylims!(pb, -4.0, 4.1)
# annotate!([
#     (0.9, r1[end]+0.1,   text("λ=0.0", 10, :left, :black)),
#     # (0.9, r1[end],   text("λ=0.0", 10, :left, :black)),
#     # (0.9, r5[end],   text("λ=0.8", 10, :left, :black)),   # optional, keeps it clean
# ])
end