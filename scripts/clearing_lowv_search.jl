# After fault clearing, does the solver settle on a low-voltage solution while a
# high-voltage solution exists? And which re-initializations recover it?
#
#   julia --project=. scripts/clearing_lowv_search.jl mode=rms form=power lf=0.99 [xf=...] [tfault=...] [run=...]
#
# For every clearing event (fault removed, or cleared by tripping the faulted line):
#   1. inherited state = end of the fault-on trajectory (homotopy application + integration)
#   2. what the simulator does: Newton re-init from the inherited state, and the
#      first step of direct integration at every h in H_DIRECT
#   3. candidate solutions from independent starts (starts=, default hc,flat):
#      homotopy in the clearing parameter (hc), Newton from a flat start V = 1∠0
#      (flat), optionally from |V| = 1∠θ_inherited (flat_angle) and from the
#      pre-fault voltages (prefault_v).
#      Each candidate must satisfy KCL (physical). HV = the physical candidate with
#      the highest minimum bus voltage.
#   4. low-voltage flag: the simulator's state differs from HV (bus voltages) and its
#      minimum voltage is lower by more than DV_LOW. Spurious if its KCL mismatch is
#      large, otherwise a genuine low-voltage root.
#   5. persistence: integrate from Newton's root and from HV for T_HOLD (hold=) and
#      compare the voltage at the low-voltage bus. CLEAN = Newton's run stays at
#      V < 0.05 on that bus while the HV run stays synchronised (angle spread < 180
#      deg) and physical (KCL mismatch < KCL_TOL) for the whole hold.
#   6. recovery: which of HC / flat / flat-angle / pre-fault reach HV.
# Machinery (case build, events, homotopy form of the events, physicality check) is
# shared with scripts/reinit_search.jl. Output: outputs/clearing_lowv/<mode>_<form>_<run>_<time>/cases.csv

push!(ARGS, "families=fault,faultline")            # clearing events only (fault family also yields apply, filtered below)
include("reinit_search.jl")
using Printf, Dates, LinearAlgebra

const DV_LOW = parse(Float64, get(args, "dvlow", "0.05"))
const STARTS = split(get(args, "starts", "hc,flat"), ",")     # recovery candidates: hc, flat, flat_angle, prefault_v
const T_HOLD, H_HOLD = MODE == :rms ? (parse(Float64, get(args, "hold", "0.3")), 1e-2) : (parse(Float64, get(args, "hold", "0.05")), 5e-4)

vbus(c, u) = vmag(c, u)
vdist(c, u1, u2) = norm(vbus(c, u1) - vbus(c, u2), Inf)
function with_voltages(c, u, vd, vq)
    x = copy(u); x[c.address["balance_d"]] = vd; x[c.address["balance_q"]] = vq; x
end
flat0(c, u) = with_voltages(c, u, ones(length(c.nsb)), zeros(length(c.nsb)))
function flat_angle(c, u)
    θ = atan.(u[c.address["balance_q"]], u[c.address["balance_d"]])
    with_voltages(c, u, cos.(θ), sin.(θ))
end
prefault_v(c, u) = with_voltages(c, u, c.u0[c.address["balance_d"]], c.u0[c.address["balance_q"]])

spread_deg(c, x) = (d = x[c.address["delta"]]; (maximum(d) - minimum(d))*180/pi)
# Integrate the post-event model from u for T_HOLD; return the metrics used by the
# clean-case filter (angle spread and KCL mismatch along the whole run).
function hold_metrics(c, f, fk, u, ib)
    ok, rc, tf, s = integrate(c, f, u, 1.0, H_HOLD, T_HOLD)
    s === nothing && return (; ok=false, rc, tend=NaN, vend=NaN, maxspread=NaN, maxkcl=NaN)
    (; ok, rc, tend=s.time[end], vend=vmag(c, s.u[end])[ib],
       maxspread=maximum(spread_deg(c, x) for x in s.u), maxkcl=maximum(kcl_mismatch(c, fk, x) for x in s.u))
end

function lowv_case(c, ev, lf_tag)
    c.models.fault.bus[1] = ev.fault_bus; c.models.fault.x_fault[1] = ev.xf
    f = event_model(c, ev); fk = event_model(c, ev; form=:kcl); phys(u) = kcl_mismatch(c, fk, u)
    fa = event_model(c, Event(family=:fault, label="", fault_bus=ev.fault_bus, xf=ev.xf, fault_mode=:apply))
    ok, ua, _, _ = homotopy_reinit(c, fa, c.u0); ok || return nothing
    ok, _, _, sol = integrate(c, fa, ua, 1.0, H_FAULTON, T_FAULT); ok || return nothing
    u_inh = copy(sol.u[end])
    kb = findfirst(==(ev.fault_bus), c.nsb)

    # candidate solutions of the post-clearing algebraic system
    cand = Dict{String,Any}()
    nconv, uN, nit = newton_reinit(c, f, u_inh); cand["newton"] = (nconv, uN, nit)
    if "hc" in STARTS
        hconv, uH, hl, _ = homotopy_reinit(c, f, u_inh); cand["hc"] = (hconv, uH, 0)
    end
    for (name, mk) in (("flat", flat0), ("flat_angle", flat_angle), ("prefault_v", prefault_v))
        name in STARTS || continue
        cv, ux, it = newton_reinit(c, f, mk(c, u_inh)); cand[name] = (cv, ux, it)
    end
    info = Dict(k => (conv = v[1], u = v[2], it = v[3], kcl = v[1] ? phys(v[2]) : NaN,
                      minv = v[1] ? minimum(vbus(c, v[2])) : NaN) for (k, v) in cand)
    physical = [k for k in vcat(STARTS, "newton") if haskey(info, k) && info[k].conv && info[k].kcl < KCL_TOL]
    hv_src = isempty(physical) ? "" : physical[argmax([info[k].minv for k in physical])]
    uHV = hv_src == "" ? nothing : info[hv_src].u
    reaches_hv(k) = uHV !== nothing && haskey(info, k) && info[k].conv && info[k].kcl < KCL_TOL && vdist(c, info[k].u, uHV) < 1e-3
    is_low(u, kcl) = uHV !== nothing && vdist(c, u, uHV) > 1e-3 && minimum(vbus(c, u)) < minimum(vbus(c, uHV)) - DV_LOW

    # direct integration: state after the first step at each h
    dcols = String[]; direct_low = false
    for h in H_DIRECT
        okd, rc, tf, sd = integrate(c, f, u_inh, 1.0, h, N_DIRECT*h)
        if sd === nothing || length(sd.u) < 2
            push!(dcols, "fail", "NaN", "NaN"); continue
        end
        x1 = sd.u[2]; k1 = phys(x1); low1 = is_low(x1, k1); direct_low |= low1
        push!(dcols, okd ? (low1 ? (k1 < KCL_TOL ? "lowV" : "spurious") : "ok") : string(rc),
              @sprintf("%.3f", minimum(vbus(c, x1))), @sprintf("%.1e", k1))
    end

    n = info["newton"]
    newton_low = n.conv && is_low(n.u, n.kcl)
    nm = hm = nothing; clean = false
    lv_bus = newton_low ? c.nsb[argmin(vbus(c, n.u))] : 0
    # persistence: integrate from Newton's root and from HV
    persist = ""; vb_n = vb_h = NaN; kcl_n_end = NaN
    if newton_low
        ib = findfirst(==(lv_bus), c.nsb)
        okn, rcn, _, sn = integrate(c, f, n.u, 1.0, H_HOLD, T_HOLD)
        okh, _, _, sh = integrate(c, f, uHV, 1.0, H_HOLD, T_HOLD)
        vb_n = sn === nothing ? NaN : vbus(c, sn.u[end])[ib]; kcl_n_end = sn === nothing ? NaN : phys(sn.u[end])
        nm = hold_metrics(c, f, fk, n.u, ib); hm = hold_metrics(c, f, fk, uHV, ib)
        # clean: Newton's run stays on the spurious/low branch while the HV run stays
        # synchronised (spread < 180 deg) and physical for the whole hold
        clean = nm.ok && nm.vend < 0.05 && hm.ok && hm.maxspread < 180 && hm.maxkcl < KCL_TOL
        vb_h = sh === nothing ? NaN : vbus(c, sh.u[end])[ib]
        persist = sn === nothing ? "newton_run_failed($rcn)" :
                  (isfinite(vb_h) && vb_n < vb_h - 0.1) ? "stays_low" : "recovers"
        okh || (persist *= "; HV run failed")
    end
    root_type = !newton_low ? "" : n.kcl < KCL_TOL ? "physical_lowV" : "spurious"
    vals = [lf_tag, @sprintf("%.4f", c.lf), string(MODE), string(FORM), csvq(ev.label), ev.fault_bus, ev.xf, T_FAULT,
            @sprintf("%.3f", vbus(c, u_inh)[kb]), @sprintf("%.3f", minimum(vbus(c, u_inh))), dcols..., direct_low]
    for k in vcat("newton", STARTS)
        i = info[k]; push!(vals, i.conv, i.it, @sprintf("%.3f", i.minv), @sprintf("%.1e", i.kcl), reaches_hv(k))
    end
    push!(vals, hv_src != "", hv_src, uHV === nothing ? "NaN" : @sprintf("%.3f", minimum(vbus(c, uHV))),
          newton_low, root_type, lv_bus, persist, @sprintf("%.3f", vb_n), @sprintf("%.3f", vb_h), @sprintf("%.1e", kcl_n_end))
    for m in (nm, hm)
        m === nothing ? push!(vals, "", "NaN", "NaN", "NaN", "NaN") :
            push!(vals, m.ok, @sprintf("%.3f", m.tend), @sprintf("%.3f", m.vend), @sprintf("%.0f", m.maxspread), @sprintf("%.1e", m.maxkcl))
    end
    push!(vals, clean)
    return (; row = join(vals, ","), newton_low, direct_low, persist, root_type, clean,
              recovered_by = [k for k in STARTS if reaches_hv(k)])
end

function lowv_main()
    run = get(args, "run", "")
    outdir = joinpath("outputs", "clearing_lowv", join(filter(!isempty, [string(MODE), string(FORM), run, Dates.format(now(), "yyyymmdd_HHMMSS")]), "_"))
    mkpath(outdir)
    lfmax = any(!=("base"), LF_SPEC) ? lf_max_search() : NaN
    hdr = "lf,lf_abs,mode,form,label,fault_bus,xf,tfault,Vfb_preclear,minV_preclear," *
          join(["direct_h$(h)_state,direct_h$(h)_minV,direct_h$(h)_kcl" for h in H_DIRECT], ",") * ",direct_low," *
          join(["$(k)_conv,$(k)_it,$(k)_minV,$(k)_kcl,$(k)_reaches_hv" for k in vcat("newton", STARTS)], ",") *
          ",hv_exists,hv_source,hv_minV,newton_low,newton_root,lv_bus,persist,V_lvbus_end_newton,V_lvbus_end_hv,kcl_end_newton" *
          ",N_ok,N_tend,N_vend,N_maxspread,N_maxkcl,HV_ok,HV_tend,HV_vend,HV_maxspread,HV_maxkcl,clean"
    io = open(joinpath(outdir, "cases.csv"), "w"); println(io, hdr); flush(io)
    open(joinpath(outdir, "config.txt"), "w") do cf
        println(cf, "mode=$MODE form=$FORM loads=$LOAD_KW xf=$XF_LIST tfault=$T_FAULT kcl_tol=$KCL_TOL dv_low=$DV_LOW starts=$STARTS")
        println(cf, "h_direct=$H_DIRECT, hold=$T_HOLD s at h=$H_HOLD, LF_max=$lfmax, lf=$LF_SPEC")
    end
    for tag in LF_SPEC
        c = build_case(tag == "base" ? 1.0 : parse(Float64, tag)*lfmax)
        c === nothing && continue
        evs = filter(e -> e.fault_mode == :clear, make_events(c))
        cnt = Dict{String,Int}()
        for (k, ev) in enumerate(evs)
            r = try lowv_case(c, ev, tag) catch e; nothing end
            r === nothing && (cnt["skipped"] = get(cnt, "skipped", 0) + 1; continue)
            println(io, r.row); flush(io)
            key = r.newton_low ? "newton_low($(r.root_type), $(r.persist))" : r.direct_low ? "direct_low_only" : "no_low"
            cnt[key] = get(cnt, key, 0) + 1
            r.newton_low && (@printf "  [%s] %s: Newton -> %s, %s; HV recovered by %s%s\n" tag ev.label r.root_type r.persist string(r.recovered_by) (r.clean ? "  <-- CLEAN" : ""); flush(stdout))
        end
        @printf "LF %s: %d clearing events, %s\n" tag length(evs) string(cnt); flush(stdout)
    end
    close(io); println("done: $outdir")
end

lowv_main()
