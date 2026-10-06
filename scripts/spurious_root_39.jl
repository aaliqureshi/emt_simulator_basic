# Fault clearing on the Barq 39-bus in RMS form (algebraic lines and voltages),
# loads without PQBRAK, fault at a zero-injection bus. Newton from the inherited
# point vs. natural homotopy on the same (power-mismatch) equations, then
# integration from each root.
using Barq, MyDiffEq, LinearAlgebra, SparseArrays, Printf

fbus = parse(Int, get(ARGS, 1, "1"))
loadname = get(ARGS, 2, "zip")
load_kw = loadname == "zip" ? (; zip=(0.7, 0.1, 0.2), low_voltage=false) : (; zip=(0.0, 0.0, 1.0), low_voltage=false)
t_fault = 0.25; t_end = 1.5; h = parse(Float64, get(ARGS, 3, "0.01"))

models = load_data("cases/Fault_Cases/ieee39_fault.xlsx"); sys = build_system(models)
models.fault.bus = [fbus]
solve_power_flow!(sys); run_static_init!(sys)
address = build_dynamic_address(sys); n = length(build_initial_conditions(sys, address))
# RMS mass matrix: only delta and omega are differential
M = spzeros(n, n)
for i in address["delta"]; M[i, i] = 1.0; end
for (j, i) in enumerate(address["omega"]); M[i, i] = models.generator.M[j]; end
# solver address: everything after delta, omega is algebraic
alg_addr = Dict("delta"=>address["delta"], "omega"=>address["omega"], "line_id"=>1:0, "line_iq"=>1:0, "balance_d"=>1:0, "balance_q"=>1:0)
alg = (length(address["delta"]) + length(address["omega"]) + 1):n
nsb = sys.non_slack_buses; bd, bq = address["balance_d"], address["balance_q"]
Vmag(u) = (v = copy(models.bus.v); v[nsb] = hypot.(u[bd], u[bq]); v)
spread(u) = (d = u[address["delta"]]; (maximum(d) - minimum(d))*180/pi)
pb = (address, sys.models, sys.incidence_matrix, sys.C_eq, nsb)

# Power-mismatch bus balance (the formulation under test), kept here so the script
# does not depend on edits to src/models/bus.jl. ZIP load without low-voltage switch.
function balance_power!(du, u, p; zip)
    T = eltype(u); k_z, k_i, k_p = zip
    address, models, incidence_matrix, C_eq, non_slack_buses, _ = p
    bus, generator, fault, load = models.bus, models.generator, models.fault, models.load
    vd = Vector{T}(bus.vd); vq = Vector{T}(bus.vq)
    vd[non_slack_buses] = u[address["balance_d"]]; vq[non_slack_buses] = u[address["balance_q"]]
    delta = u[address["delta"]]; gid = u[address["gen_id"]]; giq = u[address["gen_iq"]]
    id = zeros(T, length(bus.idx)); iq = zeros(T, length(bus.idx))
    id[generator.bus] += @. gid*sin(delta) + giq*cos(delta)
    iq[generator.bus] += @. giq*sin(delta) - gid*cos(delta)
    id[fault.bus] -= u[address["fault_id"]]; iq[fault.bus] -= u[address["fault_iq"]]
    id .+= incidence_matrix * u[address["line_id"]]; iq .+= incidence_matrix * u[address["line_iq"]]
    w = 2*pi*60
    id .+= @. w*C_eq*vq; iq .-= @. w*C_eq*vd
    ph = @. id*vd + iq*vq; qh = @. id*vq - iq*vd
    for (k, b) in enumerate(load.bus)
        v = hypot(vd[b], vq[b]); v0 = bus.v[b]
        lf = k_z*(v/v0)^2 + k_i*(v/v0) + k_p
        ph[b] -= load.p[k]*lf; qh[b] -= load.q[k]*lf
    end
    du[address["balance_d"]] = ph[non_slack_buses]; du[address["balance_q"]] = qh[non_slack_buses]
end
model(du, u, p, t) = (solve_generator!(du, u, p); solve_line!(du, u, p, t); solve_fault!(du, u, p, t); balance_power!(du, u, p; zip=load_kw.zip))
reverse_model(du, u, p, t) = model(du, u, (p[1:end-1]..., 1 - p[end]), t)   # homotopy 1 -> 0 (clearing)
function sim(u, lam, T)
    try
        s = MyDiffEq.Solve(MyDiffEq.ODEProblem(model, u, (0.0, T), (pb..., lam), M), h, method=:Trap, adaptive=false, tstops=[], always_new=true)
        return s, s.retcode
    catch e
        return nothing, :error
    end
end
resid(u, lam) = (du = zeros(n); model(du, u, (pb..., lam), 0.0); norm(du[alg], Inf))

for xf in (length(ARGS) >= 4 ? parse.(Float64, split(ARGS[4], ",")) : (0.005, 0.01, 0.015, 0.02, 0.03))
    models.fault.x_fault[1] = xf
    u0 = build_initial_conditions(sys, address)
    ua = copy(u0); ra = solve_newton!(ua, (pb..., 1.0), alg_addr; max_iter=50, always_new=true, model! = model)
    line = @sprintf "bus %d %s xf=%.3f | prefault resid %.1e | app Newton conv=%s it=%d V%d=%.3f |" fbus loadname xf resid(u0, 0.0) ra.converged ra.iters fbus Vmag(ua)[fbus]
    if !ra.converged; println(line * " (application failed)"); continue; end
    sf, rf = sim(ua, 1.0, t_fault)
    if rf != :Success; println(line * " fault-on $rf"); continue; end
    uc = sf.u[end]
    un = copy(uc); rn = solve_newton!(un, (pb..., 0.0), alg_addr; max_iter=50, always_new=true, model! = model)
    uh = copy(uc); rh = solve_homotopy!(uh, pb, alg_addr; Δλ=0.02, always_new=true, model! = reverse_model)
    line *= @sprintf " clear: V%d pre %.3f | Newton conv=%s it=%d V%d=%.4f minV=%.3f | homotopy conv=%s V%d=%.4f minV=%.3f |" fbus Vmag(uc)[fbus] rn.converged rn.iters fbus Vmag(un)[fbus] minimum(Vmag(un)) rh.converged fbus Vmag(uh)[fbus] minimum(Vmag(uh))
    for (tag, u) in (("N", un), ("H", uh))
        s, r = sim(u, 0.0, t_end - t_fault)
        line *= s === nothing ? " $tag: $r |" : @sprintf(" %s: %s t=%.2f V%d_end=%.3f spread %.0f->%.0f deg |", tag, r, t_fault + s.time[end], fbus, Vmag(s.u[end])[fbus], spread(u), spread(s.u[end]))
    end
    println(line)
end
