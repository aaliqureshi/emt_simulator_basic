# Trace one search case after the event: which bus goes to V = 0, when, and how
# the physical trajectory (small h) behaves.
#   julia --project=. scripts/trace_case.jl mode=rms form=power lf=0.99 "match=fault bus 6 xf=0.001 clear@0.1"
include("reinit_search.jl")
using Printf
tag = LF_SPEC[1]; lfmax = tag == "base" ? NaN : lf_max_search()
c = build_case(tag == "base" ? 1.0 : parse(Float64, tag)*lfmax)
ev = make_events(c)[1]; println("event: ", ev.label, "  LF = ", c.lf)
c.models.fault.bus[1] = ev.fault_bus > 0 ? ev.fault_bus : c.nsb[1]; c.models.fault.x_fault[1] = isnan(ev.xf) ? 0.01 : ev.xf
f = event_model(c, ev); fk = event_model(c, ev; form=:kcl)
u = copy(c.u0)
if ev.fault_mode == :clear
    fa = event_model(c, Event(family=:fault, label="", fault_bus=ev.fault_bus, xf=ev.xf, fault_mode=:apply))
    _, ua, _, _ = homotopy_reinit(c, fa, c.u0); _, _, _, sol = integrate(c, fa, ua, 1.0, H_FAULTON, T_FAULT); u = copy(sol.u[end])
end
spread(x) = (d = x[c.address["delta"]]; (maximum(d) - minimum(d))*180/pi)
bus_of(i) = c.nsb[i]
for (h, T) in ((1e-2, 0.3), (1e-3, 0.3), (1e-4, 0.3))
    ok, rc, tf, s = integrate(c, f, u, 1.0, h, T)
    println(@sprintf("\nh = %g: %s, t_end = %.4f", h, rc, s === nothing ? NaN : s.time[end]))
    s === nothing && continue
    for k in unique(round.(Int, range(1, length(s.u), length=12)))
        x = s.u[k]; v = vmag(c, x); i = argmin(v)
        @printf "  t=%.3f  kcl=%.1e  min V=%.4f at bus %d  V%d=%.3f  angle spread=%.0f deg\n" s.time[k] kcl_mismatch(c, fk, x) v[i] bus_of(i) ev.fault_bus v[findfirst(==(ev.fault_bus), c.nsb)] spread(x)
    end
end
