# PowerSAS clearing re-init with PQBRAK added to the constant-power load part.
# The inherited pre-clearing state is PowerSAS's (fault-on trajectory without PQBRAK);
# only the post-clearing algebraic equations get the low-voltage factor
# S_i -> S_i * kP(|V_i|) (OpenIPSL characteristic 1, PQBRAK = 0.7 pu).
#   julia --project=. scripts/powersas39/pqbrak_reinit.jl <export dir> <zf> [pqbrak] [hc=true]
#   optional: zip=kz,ki,kp (re-split the constant-power loads), hc=false
# Works with the dense 39-bus exports (nr_call1.json + he_clear.json) and the sparse
# Polish export (clrNR/nr_call1.json + he/states.json).
using Barq, JSON, LinearAlgebra, SparseArrays, Printf
const BM = Barq.Models.BusModel
vecf(x) = Float64.(x isa AbstractVector ? x : [x]); matf(x) = permutedims(reduce(hcat, [Float64.(r) for r in x]))
cvec(d, a, b) = vecf(d[a]) + im*vecf(d[b])
dir, zf = ARGS[1], parse(Float64, ARGS[2]); pqbrak = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.7
do_hc = !("hc=false" in ARGS)
# zip=kz,ki,kp re-splits each exported constant-power load S0 into ZIP parts on the
# base |V_HE| (post-clearing HE voltages), so the HE root stays an exact root.
zarg = filter(a -> startswith(a, "zip="), ARGS)
kz, ki, kp = isempty(zarg) ? (0.0, 0.0, 1.0) : Tuple(parse.(Float64, split(zarg[1][5:end], ",")))
sparse_fmt = isfile(joinpath(dir, "clrNR", "nr_call1.json"))
ex = JSON.parsefile(sparse_fmt ? joinpath(dir, "clrNR", "nr_call1.json") : joinpath(dir, "nr_call1.json"))
if sparse_fmt
    nb = Int(ex["nbus"]); Y = sparse(Int.(vecf(ex["Yi_idx"])), Int.(vecf(ex["Yj_idx"])), vecf(ex["Yv_r"]) + im*vecf(ex["Yv_i"]), nb, nb)
    st = JSON.parsefile(joinpath(dir, "he", "states.json")); fb = Int(st["fb"]); Vhe = cvec(st, "Vc_post_r", "Vc_post_i")
else
    Y = sparse(matf(ex["Yr"]) + im*matf(ex["Yi"])); nb = size(Y, 1); fb = 1
    he = JSON.parsefile(joinpath(dir, "he_clear.json")); Vhe = cvec(he, "Vr", "Vi")
end
S = cvec(ex, "Sr", "Si"); I = cvec(ex, "Ir", "Ii"); gb = Int.(vecf(ex["syn_bus"])); GV = matf(ex["MatGV"]); GR = matf(ex["MatGRhs"])
V0 = cvec(ex, "V0r", "V0i"); sw = Int.(vecf(get(ex, "sw_bus", Float64[]))); free = setdiff(1:nb, sw); nf = length(free)
shunt = spzeros(ComplexF64, nb, nb); shunt[fb, fb] = 1/(im*zf); Yprev = Y + shunt      # clearing: previous = fault on
function gen_current(Vr, Vi)
    T = promote_type(eltype(Vr), eltype(Vi)); Igr = zeros(T, nb); Igi = zeros(T, nb)
    for (k, b) in enumerate(gb)
        Igr[b] += GR[k, 1] - (GV[k, 1]*Vr[b] + GV[k, 2]*Vi[b]); Igi[b] += GR[k, 2] - (GV[k, 3]*Vr[b] + GV[k, 4]*Vi[b])
    end
    Igr + im*Igi
end
kP(v) = pqbrak > 0 ? BM._openipsl_load_factors(v, pqbrak, 1)[1] : one(v)
Vb = abs.(Vhe)
function residual(V, Yl, form)
    vm = sqrt.(abs2.(V) .+ 1e-24)
    Sload = conj.(S) .* (kz .* (vm ./ Vb).^2 .+ ki .* (vm ./ Vb) .+ kp .* kP.(vm))
    d = Sload .+ I .* vm .+ gen_current(real.(V), imag.(V)) .* conj.(V) .- (Yl*V) .* conj.(V)
    form == :power ? d : d ./ conj.(V)
end
kcl(V) = abs.(residual(V, Y, :current))
make_model(form) = (du, u, p, t) -> begin
    T = eltype(u); Vr = Vector{T}(real.(V0)); Vi = Vector{T}(imag.(V0)); Vr[free] = u[1:nf]; Vi[free] = u[nf+1:2nf]
    r = residual(Vr + im*Vi, Yprev + p[end]*(Y - Yprev), form); du[1:nf] = real.(r[free]); du[nf+1:2nf] = imag.(r[free])
end
address = Dict("delta"=>1:0, "omega"=>1:0, "line_id"=>1:0, "line_iq"=>1:0, "balance_d"=>1:0, "balance_q"=>1:0)
toV(u) = (V = copy(V0); V[free] = u[1:nf] + im*u[nf+1:end]; V)
report(tag, conv, it, V) = @printf "  %-26s conv=%-5s it=%-3d |V_fb|=%.4f  min|V|=%.4f  max KCL=%.2e  dist to HE=%.2e\n" tag conv it abs(V[fb]) minimum(abs.(V)) maximum(kcl(V)) norm(V - Vhe, Inf)
@printf "%s  |Zf|=%.3f  ZIP=(%.2f, %.2f, %.2f)  PQBRAK=%s  (fault bus %d, |V_fb| before clearing %.3f)\n" dir zf kz ki kp pqbrak fb abs(V0[fb])
@printf "  HE root (PowerSAS, no PQBRAK) with PQBRAK equations: max KCL = %.2e, min|V| = %.3f\n" maximum(kcl(Vhe)) minimum(abs.(Vhe))
u0 = [real.(V0[free]); imag.(V0[free])]
for form in (:power, :current)
    u = copy(u0); r = solve_newton!(u, (address, 1.0), address; max_iter=50, always_new=true, model! = make_model(form))
    report("Newton $(form) (inherited)", r.converged, r.iters, toV(u))
end
flat = [ones(nf); zeros(nf)]
u = copy(flat); r = solve_newton!(u, (address, 1.0), address; max_iter=50, always_new=true, model! = make_model(:power))
report("Newton power (flat 1∠0)", r.converged, r.iters, toV(u))
if do_hc
    u = copy(u0); rh = redirect_stdout(devnull) do; solve_homotopy!(u, (address,), address; Δλ=0.02, always_new=true, model! = make_model(:power)) end
    report("HC power (50 steps)", rh.converged, rh.total_iters, toV(u))
end
