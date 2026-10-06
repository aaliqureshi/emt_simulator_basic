# Trace a fault-clearing event where Newton re-initialization from the inherited
# (fault-on) state does not converge, direct integration of the post-clearing model
# fails at every step size, and the homotopy re-initialization reaches the
# high-voltage root. Records per-iteration histories for the letter figure.
#
#   julia --project=. scripts/trace_divergence.jl mode=rms form=power pqbrak=1e-6 tfault=0.1 lf=base \
#         "match=fault bus 16 xf=0.001 cleared by trip 16-21" T=2.0
#
# Outputs (outputs/clearing_lowv/divergence_<label>_lf<tag>/):
#   newton_reinit.csv   k, |g|_inf, |dy|, V_faultbus, minV, maxV        Newton on the algebraic equations from the inherited state
#   direct_h<h>.csv     k, |G|_2, |dx|, V_faultbus, minV, maxV          Newton of the first implicit-Euler step (no re-init)
#   hc_path.csv         lambda, V_faultbus, minV, stage_iters           natural-parameter homotopy path
#   hv_run.csv          t, V_faultbus, minV, kcl, angle_spread_deg      integration from the homotopy root
#   summary.txt
push!(ARGS, "families=fault,faultline")
include("reinit_search.jl")
using Printf, ForwardDiff

T = parse(Float64, get(args, "T", "2.0")); h_run = parse(Float64, get(args, "h", "0.01"))
NIT = parse(Int, get(args, "nit", "50"))
tag = LF_SPEC[1]; c = build_case(tag == "base" ? 1.0 : parse(Float64, tag)*lf_max_search())
ev = filter(e -> e.fault_mode == :clear, make_events(c))[1]
c.models.fault.bus[1] = ev.fault_bus; c.models.fault.x_fault[1] = ev.xf
f = event_model(c, ev); fk = event_model(c, ev; form=:kcl)
fa = event_model(c, Event(family=:fault, label="", fault_bus=ev.fault_bus, xf=ev.xf, fault_mode=:apply))
ok, ua, _, _ = homotopy_reinit(c, fa, c.u0); ok || error("fault application re-init failed")
ok, _, _, sol = integrate(c, fa, ua, 1.0, H_FAULTON, T_FAULT); ok || error("fault-on integration failed")
u = copy(sol.u[end])
ib = findfirst(==(ev.fault_bus), c.nsb)
sol_faulton = sol
n, alg_idx = Barq.DynamicSim._setup_alg(u, c.alg_addr)
p1 = (c.pb..., 1.0)
spread(x) = (d = x[c.address["delta"]]; (maximum(d) - minimum(d))*180/pi)
vfb(x) = vmag(c, x)[ib]
outdir = joinpath("outputs", "clearing_lowv", "divergence_" * replace(ev.label, r"[^A-Za-z0-9]+" => "_") * "_lf$(tag)_$(FORM)" * (haskey(args, "tag") ? "_" * args["tag"] : ""))
mkpath(outdir)
io = open(joinpath(outdir, "summary.txt"), "w")
say(args...) = (s = sprint(print, args...); println(s); println(io, s))

say("event: ", ev.label, "  LF = ", @sprintf("%.4f", c.lf), "  loads = ", LOAD_KW, "  form = ", FORM, "  mode = ", MODE)
say(@sprintf("inherited state: V%d = %.4f, min V = %.4f, angle spread = %.1f deg, post-clearing algebraic residual |g|_inf = %.3e",
             ev.fault_bus, vfb(u), minimum(vmag(c, u)), spread(u), alg_resid(c, f, u, 1.0)))

# time axis: t = 0 at clearing; pre-fault segment of 50 ms, then the fault-on trajectory
open(joinpath(outdir, "faulton_run.csv"), "w") do o
    println(o, "t,V_faultbus,minV")
    v0 = vfb(c.u0); m0 = minimum(vmag(c, c.u0))
    println(o, join([-T_FAULT - 0.05, v0, m0], ",")); println(o, join([-T_FAULT, v0, m0], ","))
    for (k, xk) in enumerate(sol_faulton.u); println(o, join([sol_faulton.time[k] - T_FAULT, vfb(xk), minimum(vmag(c, xk))], ",")); end
end

# ---------------------------------------------------------------- Newton re-init, full iterate history
say("\n[Newton re-init from inherited state]")
x = copy(u)
open(joinpath(outdir, "newton_reinit.csv"), "w") do o
    println(o, "k,g_inf,g_2,dy_2,V_faultbus,minV,maxV,condJ")
    for k in 1:NIT
        g, J = Barq.DynamicSim._eval_g_jac(x, p1, n, alg_idx, f)
        if any(!isfinite, g) || any(!isfinite, J)
            say(@sprintf("  k=%2d non-finite residual/Jacobian (V%d = %.3e)", k, ev.fault_bus, vfb(x))); break
        end
        dy = -(J \ g)
        println(o, join([k, norm(g, Inf), norm(g), norm(dy), vfb(x), minimum(vmag(c, x)), maximum(vmag(c, x)), cond(J)], ","))
        (k <= 12 || k % 5 == 0) && say(@sprintf("  k=%2d  |g|=%.3e  |dy|=%.3e  V%d=%.4f  minV=%.4f  maxV=%.3f", k, norm(g, Inf), norm(dy), ev.fault_bus, vfb(x), minimum(vmag(c, x)), maximum(vmag(c, x))))
        any(!isfinite, dy) && (say("  non-finite Newton step"); break)
        x[alg_idx] .+= dy
        norm(dy) < 1e-12*(1 + norm(x[alg_idx])) && (say("  converged at k=$k"); break)
    end
end
uN = copy(u); rN = solve_newton!(uN, p1, c.alg_addr; max_iter=NIT, always_new=true, model! = f)
say(@sprintf("  solve_newton!: converged=%s iters=%d last residuals: %s", rN.converged, rN.iters, join([@sprintf("%.2e", r) for r in rN.residuals[max(1,end-4):end]], " ")))

# ---------------------------------------------------------------- direct integration (no re-init): first implicit-Euler step
say("\n[Direct integration, first step from inherited state]")
M = Matrix(c.M)
for h in H_DIRECT
    okd, rc, tf, sd = try
        s = MyDiffEq.Solve(MyDiffEq.ODEProblem(f, u, (0.0, N_DIRECT*h), p1, c.M), h, method=METHOD, adaptive=false, tstops=[], always_new=true)
        (s.retcode == :Success, s.retcode, s.time[end], s)
    catch e
        (false, Symbol("exception: " * sprint(showerror, e)[1:min(end, 160)]), NaN, nothing)
    end
    say(@sprintf("  h=%g: %s (t_end=%.4g)", h, rc, tf))
    # replicate the step's Newton iteration as MyDiffEq forms it: differential rows
    # M (x - u) - h f(x), algebraic rows g(x) (no h scaling)
    G(x) = (du = zeros(eltype(x), n); f(du, x, p1, h); r = M*(x - u) - h*du; r[alg_idx] .= du[alg_idx]; r)
    y = copy(u)
    open(joinpath(outdir, @sprintf("direct_h%g.csv", h)), "w") do o
        println(o, "k,G_2,G_inf,dx_2,V_faultbus,minV,maxV,condJ")
        for k in 1:NIT
            g = G(y); J = ForwardDiff.jacobian(G, y)
            if any(!isfinite, g) || any(!isfinite, J)
                say(@sprintf("    k=%2d non-finite residual/Jacobian", k)); break
            end
            dx = try -(J \ g) catch e; say("    k=$k linear solve failed: $(sprint(showerror, e))"); break end
            println(o, join([k, norm(g), norm(g, Inf), norm(dx), vfb(y), minimum(vmag(c, y)), maximum(vmag(c, y)), cond(J)], ","))
            (k <= 8 || k % 10 == 0) && say(@sprintf("    k=%2d  |G|=%.3e  |dx|=%.3e  V%d=%.4f  minV=%.4f  maxV=%.3f", k, norm(g), norm(dx), ev.fault_bus, vfb(y), minimum(vmag(c, y)), maximum(vmag(c, y))))
            y .+= dx
            norm(dx) < 1e-7*(1 + norm(y)) && (say("    converged at k=$k, V$(ev.fault_bus) = $(round(vfb(y), digits=4)), KCL = $(kcl_mismatch(c, fk, y))"); break)
        end
    end
end

# ---------------------------------------------------------------- homotopy re-init
say("\n[Homotopy re-init (clearing parameter 0 -> 1)]")
uh = copy(u)
rH = solve_homotopy!(uh, c.pb, c.alg_addr; Δλ=0.01, always_new=true, model! = f,
                     vd_idx=c.address["balance_d"][ib], vq_idx=c.address["balance_q"][ib])
say(@sprintf("  converged=%s total Newton iterations=%d lambda_failed=%s", rH.converged, rH.total_iters, string(rH.λ_failed)))
say(@sprintf("  root: V%d = %.4f, min V = %.4f, KCL = %.2e, |g|_inf = %.2e", ev.fault_bus, vfb(uh), minimum(vmag(c, uh)), kcl_mismatch(c, fk, uh), alg_resid(c, f, uh, 1.0)))
open(joinpath(outdir, "hc_path.csv"), "w") do o
    println(o, "lambda,V_faultbus,stage_iters")
    for (l, vd, vq, it) in zip(rH.λ_hist, rH.vd_hist, rH.vq_hist, rH.stage_iters)
        println(o, join([l, hypot(vd, vq), it], ","))
    end
end

# flat start for comparison
flat = copy(u); flat[c.address["balance_d"]] .= 1.0; flat[c.address["balance_q"]] .= 0.0
uF = copy(flat); rF = solve_newton!(uF, p1, c.alg_addr; max_iter=NIT, always_new=true, model! = f); okF = rF.converged; itF = rF.iters
say(@sprintf("  flat start Newton: converged=%s iters=%d, distance to homotopy root (alg states, inf) = %.2e", okF, itF, okF ? norm(uF[alg_idx] - uh[alg_idx], Inf) : NaN))
open(joinpath(outdir, "flat_reinit.csv"), "w") do o
    println(o, "k,g_inf,dy_2")
    for k in eachindex(rF.correction_norm); println(o, join([k, rF.residuals[k], rF.correction_norm[k]], ",")); end
end
# Newton residual history along the homotopy, flattened (stage boundaries in hc_path.csv)
open(joinpath(outdir, "hc_residuals.csv"), "w") do o
    println(o, "k,g_inf"); for (k, r) in enumerate(rH.residuals); println(o, join([k, r], ",")); end
end

# ---------------------------------------------------------------- run from the homotopy root
say("\n[Integration from the homotopy root]")
okR, rcR, tfR, sR = integrate(c, f, uh, 1.0, h_run, T)
say(@sprintf("  h=%g T=%g: %s, t_end=%.3f", h_run, T, rcR, tfR))
if sR !== nothing
    open(joinpath(outdir, "hv_run.csv"), "w") do o
        println(o, "t,V_faultbus,minV,kcl,angle_spread_deg")
        for (k, xk) in enumerate(sR.u)
            println(o, join([sR.time[k], vfb(xk), minimum(vmag(c, xk)), kcl_mismatch(c, fk, xk), spread(xk)], ","))
        end
    end
    say(@sprintf("  max angle spread = %.1f deg, min V over run = %.3f, max KCL = %.1e, V%d(end) = %.3f",
                 maximum(spread(xk) for xk in sR.u), minimum(minimum(vmag(c, xk)) for xk in sR.u),
                 maximum(kcl_mismatch(c, fk, xk) for xk in sR.u), ev.fault_bus, vfb(sR.u[end])))
end
# ---------------------------------------------------------------- direct post-clearing run (no re-init) at h_run for T
say("\n[Direct integration, no re-init, h=$(h_run) for T=$(T)]")
okD, rcD, tfD, sD = integrate(c, f, u, 1.0, h_run, T)
say(@sprintf("  %s, t_end=%.3f", rcD, tfD))
if sD !== nothing && length(sD.u) > 1
    open(joinpath(outdir, "direct_run.csv"), "w") do o
        println(o, "t,V_faultbus,minV,kcl,angle_spread_deg")
        for (k, xk) in enumerate(sD.u); println(o, join([sD.time[k], vfb(xk), minimum(vmag(c, xk)), kcl_mismatch(c, fk, xk), spread(xk)], ",")); end
    end
    say(@sprintf("  V%d(end) = %.3f, min V(end) = %.3f, KCL(end) = %.1e", ev.fault_bus, vfb(sD.u[end]), minimum(vmag(c, sD.u[end])), kcl_mismatch(c, fk, sD.u[end])))
end

# ---------------------------------------------------------------- first implicit-Euler step after a re-initialization
# Newton residual history of the first step M (x - u_root) - h f(x) = 0 (algebraic rows g(x) = 0)
function first_step_history(u_root, h)
    Gh(x) = (du = zeros(eltype(x), n); f(du, x, p1, h); r = M*(x - u_root) - h*du; r[alg_idx] .= du[alg_idx]; r)
    y = copy(u_root); hist = Float64[]
    for k in 1:NIT
        g = Gh(y); J = ForwardDiff.jacobian(Gh, y); push!(hist, norm(g))
        (any(!isfinite, g) || any(!isfinite, J)) && break
        dx = -(J \ g); y .+= dx
        norm(dx) < 1e-7*(1 + norm(y)) && (push!(hist, norm(Gh(y))); break)
    end
    hist
end
for (label, u_root, prefix) in (("homotopy", uh, "post_hc_step"), ("flat-start", okF ? uF : nothing, "post_flat_step"))
    u_root === nothing && continue
    say("\n[First DAE step from the $(label) root]")
    for h in unique(vcat(collect(H_DIRECT), h_run))
        hist = first_step_history(u_root, h)
        open(joinpath(outdir, @sprintf("%s_h%g.csv", prefix, h)), "w") do o
            println(o, "k,G_2"); for (k, r) in enumerate(hist); println(o, join([k, r], ",")); end
        end
        say(@sprintf("  h=%g: %d iterations, residuals %s", h, length(hist) - 1, join([@sprintf("%.1e", r) for r in hist], " ")))
    end
end

# ---------------------------------------------------------------- runs from Newton's root and from the flat-start root (if they differ from the homotopy root)
function save_run(name, u_start, label)
    okX, rcX, tfX, sX = integrate(c, f, u_start, 1.0, h_run, T)
    say(@sprintf("\n[Integration from the %s root] h=%g T=%g: %s, t_end=%.3f", label, h_run, T, rcX, tfX))
    sX === nothing && return
    open(joinpath(outdir, name), "w") do o
        println(o, "t,V_faultbus,minV,kcl,angle_spread_deg")
        for (k, xk) in enumerate(sX.u)
            println(o, join([sX.time[k], vfb(xk), minimum(vmag(c, xk)), kcl_mismatch(c, fk, xk), spread(xk)], ","))
        end
    end
    say(@sprintf("  max angle spread = %.1f deg, min V over run = %.3f, max KCL = %.1e, V%d(end) = %.3f",
                 maximum(spread(xk) for xk in sX.u), minimum(minimum(vmag(c, xk)) for xk in sX.u),
                 maximum(kcl_mismatch(c, fk, xk) for xk in sX.u), ev.fault_bus, vfb(sX.u[end])))
end
if rN.converged && norm(uN[alg_idx] - uh[alg_idx], Inf) > 1e-6
    say(@sprintf("\nNewton root: V%d = %.4f, min V = %.4f, KCL = %.2e (differs from homotopy root)", ev.fault_bus, vfb(uN), minimum(vmag(c, uN)), kcl_mismatch(c, fk, uN)))
    save_run("lv_run.csv", uN, "Newton (low-voltage)")
end
if okF && norm(uF[alg_idx] - uh[alg_idx], Inf) > 1e-6
    say(@sprintf("\nflat-start root: V%d = %.4f, min V = %.4f, KCL = %.2e (differs from homotopy root)", ev.fault_bus, vfb(uF), minimum(vmag(c, uF)), kcl_mismatch(c, fk, uF)))
    norm(uF[alg_idx] - uN[alg_idx], Inf) < 1e-6 && say("  identical to Newton's root")
end
okF && save_run("flat_run.csv", uF, "flat-start")
close(io)
println("\nsaved to ", outdir)
