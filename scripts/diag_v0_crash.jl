# Why does integration from power-form Newton's V = 0 root throw?
push!(ARGS, "families=faultline")
include("reinit_search.jl")
using ForwardDiff, Printf
tag = LF_SPEC[1]; lfmax = lf_max_search()
c = build_case(parse(Float64, tag)*lfmax)
ev = filter(e -> e.fault_mode == :clear, make_events(c))[1]; println("event: ", ev.label)
c.models.fault.bus[1] = ev.fault_bus; c.models.fault.x_fault[1] = ev.xf
f = event_model(c, ev)
fa = event_model(c, Event(family=:fault, label="", fault_bus=ev.fault_bus, xf=ev.xf, fault_mode=:apply))
_, ua, _, _ = homotopy_reinit(c, fa, c.u0); _, _, _, sol = integrate(c, fa, ua, 1.0, H_FAULTON, T_FAULT); u = copy(sol.u[end])
_, uN, it = newton_reinit(c, f, u)
v = vmag(c, uN); i = argmin(v); b = c.nsb[i]
@printf "Newton root: %d it, min |V| = %.3e at bus %d (vd = %.3e, vq = %.3e)\n" it v[i] b uN[c.address["balance_d"][i]] uN[c.address["balance_q"][i]]
p = (c.pb..., 1.0)
J = ForwardDiff.jacobian(x -> (du = zeros(eltype(x), length(x)); f(du, x, p, 0.0); du), uN)
bad = findall(!isfinite, J)
@printf "Jacobian: %d non-finite entries\n" length(bad)
if !isempty(bad)
    rows = unique(first.(Tuple.(bad))); cols = unique(last.(Tuple.(bad)))
    name(k) = first([key for (key, r) in c.address if k in r])
    println("  rows: ", unique(name.(rows)), "  cols: ", unique(name.(cols)))
end
fk = event_model(c, ev; form=:kcl)
for T in (0.05, 0.3)
    try
        s = MyDiffEq.Solve(MyDiffEq.ODEProblem(f, uN, (0.0, T), p, c.M), 1e-2, method=:Euler, adaptive=false, tstops=[], always_new=true)
        @printf "T = %.2f: retcode %s, t_end %.3f\n" T s.retcode s.time[end]
        for k in eachindex(s.u)
            x = s.u[k]; vv = vmag(c, x)
            @printf "   t=%.2f  V%d=%.3e  min V=%.3e  KCL=%.1e  finite=%s\n" s.time[k] b vv[i] minimum(vv) kcl_mismatch(c, fk, x) all(isfinite, x)
        end
    catch e
        bt = catch_backtrace()
        println("T = $T: integrator exception: ", first(sprint(showerror, e), 300))
        for fr in stacktrace(bt)[1:min(8, end)]; println("     ", fr); end
    end
end
