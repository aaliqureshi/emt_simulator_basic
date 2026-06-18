
const V_ref  = 1.0
const R_g    = 0.1
const P_load = 0.9
const I_max  = 3.0

const V_init = 0.9


g1(V)  = I_max - P_load / V          
dg1(V) = P_load / V^2                

const V_star = P_load / I_max        

# ---- Stabilization flow:  dV/dτ = -γ * g_y * H * g ----
#   mode = :gradient  -> H = 1            (Baumgarte / penalty, Boesen §3.4)
#   mode = :newton    -> H = 1/g_y^2      (Gauss–Newton / continuous-Newton: dg/dτ = -γ g)
function rhs(V; γ = 1.0, mode = :gradient)
    gy = dg1(V)
    g  = g1(V)
    H  = mode === :newton ? 1.0 / gy^2 : 1.0
    return -γ * gy * H * g
end

# ---- Integrate to steady state (explicit Euler) ----
function stabilize(V0; γ = 5.0, dτ = 1e-3, τmax = 20.0, tol = 1e-10, mode = :gradient)
    V = V0
    τs = [0.0]; Vs = [V]; gs = [g1(V)]
    τ = 0.0
    while τ < τmax && abs(g1(V)) > tol
        V += dτ * rhs(V; γ = γ, mode = mode)
        τ += dτ
        push!(τs, τ); push!(Vs, V); push!(gs, g1(V))
    end
    return τs, Vs, gs
end

# ---- Baseline: direct Newton on g1 (this is what fails) ----
function newton(V0; maxit = 10)
    V = V0; hist = [V]
    for _ in 1:maxit
        V -= g1(V) / dg1(V)
        push!(hist, V)
        (V <= 0 || !isfinite(V)) && break    # left the physical domain V>0
    end
    return hist
end

# ---- Run ----
τg, Vg, gg = stabilize(V_init; mode = :gradient)
τn, Vn, gn = stabilize(V_init; mode = :newton)
nh = newton(V_init)

println("Target V*               = ", V_star)
println("Gradient flow (H=I):  V = ", Vg[end], "  |g| = ", abs(gg[end]), "  steps = ", length(Vg))
println("Newton  flow  (H=1/gy²): V = ", Vn[end], "  |g| = ", abs(gn[end]), "  steps = ", length(Vn))
println("Direct Newton iterates  = ", nh, "  -> leaves V>0 domain (fails)")

# ---- Export trajectories to CSV (for plotting) ----
function writecsv(path, τs, Vs, gs)
    open(path, "w") do io
        println(io, "tau,V,absg")
        for i in eachindex(τs)
            println(io, τs[i], ",", Vs[i], ",", abs(gs[i]))
        end
    end
end
writecsv("case1_gradient.csv", τg, Vg, gg)
writecsv("case1_newton.csv",   τn, Vn, gn)
println("Wrote case1_gradient.csv and case1_newton.csv")
