# Fault-clearing re-init of the PowerSAS 39-bus case (d_039_mod, fault at bus 1)
# solved with Barq's solvers. Data exported from PowerSAS's own NR setup.
using Barq, JSON, LinearAlgebra, Printf

vecf(x) = Float64.(x isa AbstractVector ? x : [x])
matf(x) = permutedims(reduce(hcat, [Float64.(r) for r in x]))   # JSON rows -> matrix

function load_case(dir)
    ex = JSON.parsefile(joinpath(dir, "nr_call1.json")); he = JSON.parsefile(joinpath(dir, "he_clear.json"))
    Y = matf(ex["Yr"]) + im*matf(ex["Yi"])
    S = vecf(ex["Sr"]) + im*vecf(ex["Si"]); I = vecf(ex["Ir"]) + im*vecf(ex["Ii"])
    gb = Int.(vecf(ex["syn_bus"])); GV = matf(ex["MatGV"]); GR = matf(ex["MatGRhs"])
    V0 = vecf(ex["V0r"]) + im*vecf(ex["V0i"]); Vnr = vecf(ex["Vr"]) + im*vecf(ex["Vi"])
    Vhe = vecf(he["Vr"]) + im*vecf(he["Vi"])
    sw = Int.(vecf(get(ex, "sw_bus", Float64[])))
    return (; Y, S, I, gb, GV, GR, V0, Vnr, Vhe, sw, loop=ex["loop"], flag=ex["flag"])
end

# generator current injection, linear in V at the frozen state (PowerSAS MatGV/MatGRhs)
function gen_current(c, Vr, Vi)
    T = promote_type(eltype(Vr), eltype(Vi)); n = length(Vr)
    Igr = zeros(T, n); Igi = zeros(T, n)
    for (k, b) in enumerate(c.gb)
        Igr[b] += c.GR[k, 1] - (c.GV[k, 1]*Vr[b] + c.GV[k, 2]*Vi[b])
        Igi[b] += c.GR[k, 2] - (c.GV[k, 3]*Vr[b] + c.GV[k, 4]*Vi[b])
    end
    return Igr + im*Igi
end

# form = :power (PowerSAS residual) or :current (KCL). Y(λ) = Yf + λ (Y - Yf).
function make_model(c, free, form, Yf)
    nb = length(c.V0); nf = length(free)
    (du, u, p, t) -> begin
        λ = p[end]; T = eltype(u)
        Vr = Vector{T}(real.(c.V0)); Vi = Vector{T}(imag.(c.V0))
        Vr[free] = u[1:nf]; Vi[free] = u[nf+1:2nf]
        V = Vr + im*Vi
        Yl = Yf + λ*(c.Y - Yf)
        Ig = gen_current(c, Vr, Vi)
        d = conj.(c.S) .+ c.I .* abs.(V) .+ Ig .* conj.(V) .- (Yl*V) .* conj.(V)
        r = form == :power ? d : d ./ conj.(V)
        du[1:nf] = real.(r[free]); du[nf+1:2nf] = imag.(r[free])
    end
end

address = Dict("delta"=>1:0, "omega"=>1:0, "line_id"=>1:0, "line_iq"=>1:0, "balance_d"=>1:0, "balance_q"=>1:0)
kcl_mismatch(c, V) = abs.((conj.(c.S) .+ c.I .* abs.(V) .+ gen_current(c, real.(V), imag.(V)) .* conj.(V) .- (c.Y*V) .* conj.(V)) ./ conj.(V))

for zf in parse.(Float64, ARGS)
    dir = joinpath(@__DIR__, @sprintf("export_zf%.4f", zf)); isdir(dir) || continue
    c = load_case(dir); nb = length(c.V0); free = setdiff(1:nb, c.sw)
    Yf = copy(c.Y); Yf[1, 1] += 1/(im*zf)            # fault-on admittance: bolted-reactance shunt at bus 1
    u0 = [real.(c.V0[free]); imag.(c.V0[free])]
    @printf "\n=== |Zf| = %.4f  (%d buses, %d free, slack %s) ===\n" zf nb length(free) string(c.sw)
    @printf "PowerSAS NR at clearing: loop=%d flag=%d  |V1|=%.4f  min|V|=%.4f   HE: |V1|=%.4f min|V|=%.4f\n" c.loop c.flag abs(c.Vnr[1]) minimum(abs.(c.Vnr)) abs(c.Vhe[1]) minimum(abs.(c.Vhe))
    m = kcl_mismatch(c, c.Vnr); @printf "  KCL mismatch at PowerSAS NR point: bus 1 = %.3g pu, max elsewhere = %.3g   (at HE point: %.3g)\n" m[1] maximum(m[2:end]) maximum(kcl_mismatch(c, c.Vhe))
    for form in (:power, :current)
        f = make_model(c, free, form, Yf)
        du = zeros(2length(free)); f(du, u0, (address, 0.0), 0.0)
        u = copy(u0); r = solve_newton!(u, (address, 1.0), address; max_iter=50, always_new=true, model! = f)
        V = copy(c.V0); V[free] = u[1:length(free)] + im*u[length(free)+1:end]
        @printf "  [%-7s] fault-on residual at start %.1e | Newton: conv=%s it=%2d |V1|=%.4f min|V|=%.4f dist to HE=%.2e\n" form norm(du, Inf) r.converged r.iters abs(V[1]) minimum(abs.(V)) norm(V - c.Vhe, Inf)
        uh = copy(u0); rh = solve_homotopy!(uh, (address,), address; Δλ=0.02, always_new=true, model! = f)
        Vh = copy(c.V0); Vh[free] = uh[1:length(free)] + im*uh[length(free)+1:end]
        @printf "            Y-homotopy (50 steps): conv=%s total_it=%d |V1|=%.4f dist to HE=%.2e\n" rh.converged rh.total_iters abs(Vh[1]) norm(Vh - c.Vhe, Inf)
    end
end
