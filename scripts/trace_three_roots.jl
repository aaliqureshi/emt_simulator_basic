# Integrate the post-clearing model from Newton's root, the HC root and the flat-start
# root; report V at the faulted bus, KCL mismatch and angle spread; save a CSV.
#   julia --project=. scripts/trace_three_roots.jl mode=rms form=power pqbrak=0.7 tfault=0.25 lf=0.97 "match=fault bus 8 xf=0.001 cleared by trip 8-9" T=2.0
push!(ARGS, "families=fault,faultline")
include("reinit_search.jl")
using Printf
T = parse(Float64, get(args, "T", "2.0")); h = parse(Float64, get(args, "h", "0.01"))
tag = LF_SPEC[1]; c = build_case(tag == "base" ? 1.0 : parse(Float64, tag)*lf_max_search())
ev = filter(e -> e.fault_mode == :clear, make_events(c))[1]
c.models.fault.bus[1] = ev.fault_bus; c.models.fault.x_fault[1] = ev.xf
f = event_model(c, ev); fk = event_model(c, ev; form=:kcl)
fa = event_model(c, Event(family=:fault, label="", fault_bus=ev.fault_bus, xf=ev.xf, fault_mode=:apply))
_, ua, _, _ = homotopy_reinit(c, fa, c.u0); _, _, _, sol = integrate(c, fa, ua, 1.0, H_FAULTON, T_FAULT); u = copy(sol.u[end])
ib = findfirst(==(ev.fault_bus), c.nsb)
flat = copy(u); flat[c.address["balance_d"]] .= 1.0; flat[c.address["balance_q"]] .= 0.0
roots = Dict("newton" => newton_reinit(c, f, u)[2], "hc" => homotopy_reinit(c, f, u)[2], "flat" => newton_reinit(c, f, flat)[2])
spread(x) = (d = x[c.address["delta"]]; (maximum(d) - minimum(d))*180/pi)
@printf "%s, LF = %.4f, V%d before clearing = %.3f\n" ev.label c.lf ev.fault_bus vmag(c, u)[ib]
outcsv = joinpath("outputs", "clearing_lowv", "trace_" * replace(ev.label, r"[^A-Za-z0-9]+" => "_") * "_lf$(tag).csv")
open(outcsv, "w") do io
    println(io, "root,t,V_faultbus,minV,kcl,angle_spread_deg")
    for name in ("newton", "hc", "flat")
        ok, rc, tf, s = integrate(c, f, roots[name], 1.0, h, T)
        @printf "\n[%s] start: V%d = %.3f, KCL = %.1e -> %s, t_end = %.2f\n" name ev.fault_bus vmag(c, roots[name])[ib] kcl_mismatch(c, fk, roots[name]) rc (s === nothing ? NaN : s.time[end])
        s === nothing && continue
        for (k, x) in enumerate(s.u)
            println(io, join([name, s.time[k], vmag(c, x)[ib], minimum(vmag(c, x)), kcl_mismatch(c, fk, x), spread(x)], ","))
        end
        for k in unique(round.(Int, range(1, length(s.u), length=9)))
            x = s.u[k]; @printf "   t=%.2f  V%d=%.3f  min V=%.3f  KCL=%.1e  spread=%.0f deg\n" s.time[k] ev.fault_bus vmag(c, x)[ib] minimum(vmag(c, x)) kcl_mismatch(c, fk, x) spread(x)
        end
    end
end
println("\nsaved ", outcsv)
