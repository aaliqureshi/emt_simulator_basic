using Barq, MyDiffEq, LinearAlgebra, SparseArrays, Printf
models = load_data("cases/Fault_Cases/ieee39_fault.xlsx"); sys = build_system(models)
fb = parse(Int, get(ARGS, 1, "20")); models.fault.bus = [fb]; models.fault.x_fault[1] = parse(Float64, get(ARGS, 2, "0.009"))
solve_power_flow!(sys); run_static_init!(sys)
address = build_dynamic_address(sys); M0 = build_mass_matrix(sys, address); u0 = build_initial_conditions(sys, address)
nsb = sys.non_slack_buses; k = findfirst(==(fb), nsb); bd, bq = address["balance_d"], address["balance_q"]
pb = (address, sys.models, sys.incidence_matrix, sys.C_eq, nsb)
Vf(u) = hypot(u[bd[k]], u[bq[k]])
mk(kw) = (du, u, p, t) -> (solve_generator!(du, u, p); solve_line!(du, u, p, t); solve_fault!(du, u, p, t); balance!(du, u, p; kw...))
function run(f, u, lam, dt, T, M)
    try
        s = MyDiffEq.Solve(MyDiffEq.ODEProblem(f, u, (0.0, T), (pb..., lam), M), dt, method=:Euler, adaptive=false, tstops=[], always_new=true)
        return s, s.retcode
    catch e
        return nothing, Symbol(first(sprint(showerror, e), 50))
    end
end
function case(label, kw; C_fb=nothing)
    M = copy(M0)
    if C_fb !== nothing; M[bd[k], bd[k]] = C_fb; M[bq[k], bq[k]] = C_fb; end
    f = mk(kw)
    ur = copy(u0); r = solve_newton!(ur, (pb..., 1.0), address; max_iter=50, always_new=true, model! = f)
    line = @sprintf "%-34s C%d=%.0e | reinit conv=%s it=%d |" label fb M[bd[k], bd[k]] r.converged r.iters
    for dt in (5e-4, 5e-5, 5e-6)
        s, rc = run(f, ur, 1.0, dt, 200*dt, M)
        if s === nothing; line *= " h=$dt: $rc |"; continue; end
        v = Vf.(s.u)
        line *= rc == :Success ? @sprintf(" h=%g ok V%d %.3f->%.3f |", dt, fb, v[1], v[end]) :
                @sprintf(" h=%g %s at step %d (V%d %.3f->%.4f) |", dt, rc, length(s.u)-1, fb, v[1], v[end])
    end
    println(line)
end
default_kw = (;)                                            # what solve_dynamic_sim! uses
case("default balance! kwargs", default_kw)
case("pure Z", (; zip=(1.0, 0.0, 0.0), low_voltage=false))
case("ZIP (.7,.1,.2), no low-voltage", (; zip=(0.7, 0.1, 0.2), low_voltage=false))
case("ZIP + PQBRAK 0.7 char 1", (; zip=(0.7, 0.1, 0.2), low_voltage=true, pqbrak=0.7, characteristic=1))
for C in (1e-4, 1e-3, 1e-2)
    case("default kwargs", default_kw; C_fb=C)
end

# in-simulation re-init: ramp the fault admittance linearly over Tc, then continue
x_f = models.fault.x_fault[1]; fid, fiq = address["fault_id"][1], address["fault_iq"][1]
function mk_ramp(kw, bfun)
    (du, u, p, t) -> begin
        solve_generator!(du, u, p); solve_line!(du, u, p, t); balance!(du, u, p; kw...)
        b = bfun(t); du[fid] = b*u[bd[k]] + u[fiq]; du[fiq] = b*u[bq[k]] - u[fid]
    end
end
for Tc in (1e-4, 5e-4, 1e-3, 2e-3, 5e-3)
    sw, rw = run(mk_ramp((;), t -> clamp(t/Tc, 0, 1)/x_f + 1e-10), u0, 1.0, Tc/200, Tc, M0)
    if rw != :Success
        @printf "window Tc=%-6g default kwargs: FAILED (%s) at s=%.3f, V%d=%.3f\n" Tc rw (sw === nothing ? NaN : sw.time[end]/Tc) fb (sw === nothing ? NaN : Vf(sw.u[end])); continue
    end
    line = @sprintf "window Tc=%-6g default kwargs: ok, V%d=%.3f |" Tc fb Vf(sw.u[end])
    for dt in (5e-4, 5e-5)
        s, rc = run(mk_ramp((;), t -> 1/x_f), sw.u[end], 1.0, dt, 0.05, M0)
        line *= rc == :Success ? @sprintf(" then h=%g to 50 ms: ok V%d=%.3f min %.3f |", dt, fb, Vf(s.u[end]), minimum(Vf.(s.u))) : " then h=$dt: $rc |"
    end
    println(line)
end
