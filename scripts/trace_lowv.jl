# Follow the direct post-clearing run long enough to see whether it stays on the
# spurious low-voltage branch.
#   julia --project=. scripts/trace_lowv.jl mode=emt form=power lf=0.97 "match=fault bus 3 xf=0.001 clear" h=5e-5 T=0.02
push!(ARGS, "families=fault,faultline")
include("reinit_search.jl")
using Printf
h = parse(Float64, get(args, "h", "5e-5")); T = parse(Float64, get(args, "T", "0.02"))
tag = LF_SPEC[1]; lfmax = tag == "base" ? NaN : lf_max_search()
c = build_case(tag == "base" ? 1.0 : parse(Float64, tag)*lfmax)
ev = filter(e -> e.fault_mode == :clear, make_events(c))[1]; println("event: ", ev.label, "  LF = ", round(c.lf, digits=4), "  h = $h")
c.models.fault.bus[1] = ev.fault_bus; c.models.fault.x_fault[1] = ev.xf
f = event_model(c, ev); fk = event_model(c, ev; form=:kcl)
fa = event_model(c, Event(family=:fault, label="", fault_bus=ev.fault_bus, xf=ev.xf, fault_mode=:apply))
_, ua, _, _ = homotopy_reinit(c, fa, c.u0); _, _, _, sol = integrate(c, fa, ua, 1.0, H_FAULTON, T_FAULT); u = copy(sol.u[end])
ib = findfirst(==(ev.fault_bus), c.nsb)
_, uH, _, _ = homotopy_reinit(c, f, u)
for (name, u_start) in (("direct from inherited state", u), ("direct from HC root", uH))
    ok, rc, tf, s = integrate(c, f, u_start, 1.0, h, T)
    @printf "\n%s: %s, t_end = %.5f\n" name rc (s === nothing ? NaN : s.time[end])
    s === nothing && continue
    for k in unique(round.(Int, range(1, length(s.u), length=10)))
        x = s.u[k]; @printf "  t=%.5f  V%d=%.4f  min V=%.4f  KCL mismatch=%.1e\n" s.time[k] ev.fault_bus vmag(c, x)[ib] minimum(vmag(c, x)) kcl_mismatch(c, fk, x)
    end
end
