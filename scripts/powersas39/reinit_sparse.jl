# Re-init at a PowerSAS switching event (sparse export), solved with Barq's solvers.
# usage: julia --project=. reinit_sparse.jl <case dir> <app|clr> <zf>
using Barq, JSON, LinearAlgebra, SparseArrays, Printf

vecf(x) = Float64.(x isa AbstractVector ? x : [x])
matf(x) = permutedims(reduce(hcat, [Float64.(r) for r in x]))
cvec(d, a, b) = vecf(d[a]) + im*vecf(d[b])

casedir, ev, zf = ARGS[1], ARGS[2], parse(Float64, ARGS[3])
run = ev == "app" ? "appNR" : "clrNR"
ex = JSON.parsefile(joinpath(casedir, run, "nr_call1.json"))
st_he = JSON.parsefile(joinpath(casedir, "he", "states.json"))
nb = Int(ex["nbus"])
Y = sparse(Int.(vecf(ex["Yi_idx"])), Int.(vecf(ex["Yj_idx"])), vecf(ex["Yv_r"]) + im*vecf(ex["Yv_i"]), nb, nb)
S = cvec(ex, "Sr", "Si"); I = cvec(ex, "Ir", "Ii")
gb = Int.(vecf(ex["syn_bus"])); GV = matf(ex["MatGV"]); GR = matf(ex["MatGRhs"])
V0 = cvec(ex, "V0r", "V0i"); Vnr = cvec(ex, "Vr", "Vi")
fb = Int(st_he["fb"])
Vref = ev == "app" ? cvec(st_he, "Va_post_r", "Va_post_i") : cvec(st_he, "Vc_post_r", "Vc_post_i")
sw = Int.(vecf(get(ex, "sw_bus", Float64[])))
free = setdiff(1:nb, sw); nf = length(free)
shunt = spzeros(ComplexF64, nb, nb); shunt[fb, fb] = 1/(im*zf)
Yprev = ev == "app" ? Y - shunt : Y + shunt

function gen_current(Vr, Vi)
    T = promote_type(eltype(Vr), eltype(Vi)); Igr = zeros(T, nb); Igi = zeros(T, nb)
    for (k, b) in enumerate(gb)
        Igr[b] += GR[k, 1] - (GV[k, 1]*Vr[b] + GV[k, 2]*Vi[b])
        Igi[b] += GR[k, 2] - (GV[k, 3]*Vr[b] + GV[k, 4]*Vi[b])
    end
    Igr + im*Igi
end
function residual(V, Yl, form)
    d = conj.(S) .+ I .* abs.(V) .+ gen_current(real.(V), imag.(V)) .* conj.(V) .- (Yl*V) .* conj.(V)
    form == :power ? d : d ./ conj.(V)
end
function make_model(form)
    (du, u, p, t) -> begin
        T = eltype(u); Vr = Vector{T}(real.(V0)); Vi = Vector{T}(imag.(V0))
        Vr[free] = u[1:nf]; Vi[free] = u[nf+1:2nf]
        r = residual(Vr + im*Vi, Yprev + p[end]*(Y - Yprev), form)
        du[1:nf] = real.(r[free]); du[nf+1:2nf] = imag.(r[free])
    end
end
toV(u) = (V = copy(V0); V[free] = u[1:nf] + im*u[nf+1:end]; V)
address = Dict("delta"=>1:0, "omega"=>1:0, "line_id"=>1:0, "line_iq"=>1:0, "balance_d"=>1:0, "balance_q"=>1:0)
kcl(V) = abs.(residual(V, Y, :current))

@printf "%s event, %d buses, fault bus %d, |Zf| = %.3f\n" ev nb fb zf
@printf "start: |V_fb| = %.4f; residual at start with previous Y (sanity, should be ~HE tol): %.1e\n" abs(V0[fb]) norm(residual(V0, Yprev, :current), Inf)
@printf "PowerSAS NR: loop=%d flag=%d resid=%.2e |V_fb|=%.4f min|V|=%.4f | max KCL mismatch %.3g at bus %d | dist to HE %.2e\n" ex["loop"] ex["flag"] ex["resid"] abs(Vnr[fb]) minimum(abs.(Vnr)) maximum(kcl(Vnr)) argmax(kcl(Vnr)) norm(Vnr - Vref, Inf)
@printf "HE reference: |V_fb|=%.4f min|V|=%.4f max KCL mismatch %.2e | displacement |V_HE - V_start|_inf = %.3f\n" abs(Vref[fb]) minimum(abs.(Vref)) maximum(kcl(Vref)) norm(Vref - V0, Inf)
for form in (:power, :current)
    f = make_model(form); u0 = [real.(V0[free]); imag.(V0[free])]
    u = copy(u0); t0 = time()
    r = solve_newton!(u, (address, 1.0), address; max_iter=30, always_new=true, model! = f)
    V = toV(u)
    @printf "[%-7s] Newton: conv=%s it=%2d (%.0fs) |V_fb|=%.4f min|V|=%.4f maxKCL=%.2e dist to HE=%.2e\n" form r.converged r.iters time()-t0 abs(V[fb]) minimum(abs.(V)) maximum(kcl(V)) norm(V - Vref, Inf)
    @printf "          residual history: %s\n" join([@sprintf("%.1e", x) for x in r.residuals], " ")
end

# natural-parameter homotopy Y_prev -> Y on the same equations, from the same inherited point
for form in (:power, :current)
    f = make_model(form); u = [real.(V0[free]); imag.(V0[free])]; t0 = time()
    rh = solve_homotopy!(u, (address,), address; Δλ=0.02, always_new=true, model! = f)
    V = toV(u)
    @printf "[%-7s] Y-homotopy (50 steps): conv=%s total_it=%d (%.0fs) |V_fb|=%.4f min|V|=%.4f maxKCL=%.2e dist to HE=%.2e\n" form rh.converged rh.total_iters time()-t0 abs(V[fb]) minimum(abs.(V)) maximum(kcl(V)) norm(V - Vref, Inf)
end
