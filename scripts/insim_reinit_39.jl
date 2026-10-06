# In-simulation re-initialization: ramp the fault in over a window T_c while the
# DAE integrates (line currents free to move), then continue with a normal step.
using Barq, MyDiffEq, LinearAlgebra, Printf
models = load_data("cases/Fault_Cases/ieee39_fault.xlsx"); sys = build_system(models)
models.fault.bus = [24]; x_f = 0.015; models.fault.x_fault[1] = x_f
solve_power_flow!(sys); run_static_init!(sys)
address = build_dynamic_address(sys); M = build_mass_matrix(sys, address)
u0 = build_initial_conditions(sys, address)
nsb = sys.non_slack_buses; k = findfirst(==(24), nsb)
bd, bq = address["balance_d"], address["balance_q"]; fid, fiq = address["fault_id"][1], address["fault_iq"][1]
p = (address, sys.models, sys.incidence_matrix, sys.C_eq, nsb, 1.0)
load_kw = (; zip=(0.7, 0.1, 0.2), low_voltage=false)
V24(u) = hypot(u[bd[k]], u[bq[k]])

# fault susceptance as a function of time within the window
b_of(t, Tc, shape) = begin
    s = Tc > 0 ? clamp(t/Tc, 0.0, 1.0) : 1.0
    shape == :geom ? 1/((1e10)^(1-s) * x_f^s) : s/x_f + 1e-10
end
function mk(bfun)
    (du, u, p, t) -> begin
        solve_generator!(du, u, p); solve_line!(du, u, p, t); balance!(du, u, p; load_kw...)
        b = bfun(t)
        du[fid] = b*u[bd[k]] + u[fiq]
        du[fiq] = b*u[bq[k]] - u[fid]
    end
end
function integrate(f, u, T, dt, method)
    try
        s = MyDiffEq.Solve(MyDiffEq.ODEProblem(f, u, (0.0, T), p, M), dt, method=method, adaptive=false, tstops=[], always_new=true)
        return s, s.retcode
    catch e
        return nothing, Symbol("error: " * first(sprint(showerror, e), 60))
    end
end

for shape in (:lin, :geom), Tc in (1e-5, 1e-4, 5e-4, 1e-3, 2e-3, 5e-3)
    sw, rw = integrate(mk(t -> b_of(t, Tc, shape)), u0, Tc, Tc/200, :Euler)
    if sw === nothing || rw != :Success
        tf = sw === nothing ? NaN : sw.time[end]
        @printf "%-4s Tc=%-7g window FAILED (%s) at t=%.3g (s=%.3f)\n" shape Tc rw tf tf/Tc
        continue
    end
    uw = sw.u[end]; minw = minimum(V24.(sw.u))
    line = @sprintf "%-4s Tc=%-7g window ok: V24 min %.3f, end %.3f |" shape Tc minw V24(uw)
    for (dt, meth) in ((5e-4, :Euler), (5e-5, :Euler), (5e-4, :Trap))
        sc, rc = integrate(mk(t -> 1/x_f), uw, 0.1, dt, meth)
        line *= sc === nothing || rc != :Success ? @sprintf(" %s h=%g: %s @%.3g |", meth, dt, rc, sc === nothing ? NaN : sc.time[end]) :
                @sprintf(" %s h=%g: ok V24(0.1)=%.3f min %.3f |", meth, dt, V24(sc.u[end]), minimum(V24.(sc.u)))
    end
    println(line)
end
