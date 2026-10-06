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

# consistent point with moved currents: end of a successful 1 ms linear window
sw, _ = integrate(mk(t -> b_of(t, 1e-3, :lin)), u0, 1e-3, 1e-3/200, :Euler)
uw = sw.u[end]
diffidx = vcat(collect(address["delta"]), collect(address["omega"]), collect(address["line_id"]), collect(address["line_iq"]))
ubrown = copy(uw); ubrown[diffidx] = u0[diffidx]      # y* with the original differential states
g0 = zeros(length(u0)); mk(t -> 1/x_f)(g0, ubrown, p, 0.0)
alg = setdiff(1:length(u0), diffidx)
@printf "window end: V24=%.3f | line currents moved by max %.3f pu | alg residual at (x0, y*) = %.2e
" V24(uw) maximum(abs.(uw[diffidx] - u0[diffidx])) norm(g0[alg], Inf)
for meth in (:Euler, :Trap), dt in (2e-3, 1e-3, 5e-4, 1e-4, 5e-5, 5e-6)
    sc, rc = integrate(mk(t -> 1/x_f), ubrown, 20*dt, dt, meth)
    @printf "  restart from (x0, y*) %-5s h=%-7g -> %-8s %s
" meth dt rc (sc === nothing || rc != :Success ? "" : @sprintf("V24 after 1st step %.3f, end %.3f", V24(sc.u[2]), V24(sc.u[end])))
end
