# Event search for re-initialization failures (IEEE 39-bus first).
#
#   julia --project=. scripts/reinit_search.jl mode=rms families=line,gen,line2,topo,load,fault,faultline lf=base,0.9,0.97,0.99
#   reproduce one saved case: add  match="<label>"  and  lf=<tag>  (optional: xf=, tfault=, loadstep=, limit=)
#
# Pipeline for every (loading level, event):
#   1. Pre-event state. Steady state from power flow + static init at the loading
#      factor, or, for clearing events, the end of the fault-on trajectory.
#   2. Direct integration of the post-event model from the inherited state, with no
#      re-initialization, at every step size in the mode's h list.
#   3. If any step size fails: re-initialize from the inherited state with
#      (a) Newton on the post-event algebraic equations and
#      (b) natural-parameter homotopy in the event parameter s (0 = pre-event,
#          1 = post-event).
#   4. For every re-init that converged, integrate the post-event model from its
#      solution at every step size where the direct run failed (same horizon)
#      and over the longer H_POST checks. If it proceeds, the case is SAVED with all parameters needed to
#      reproduce it (saved_cases.csv). TARGET marks the cases the paper is after:
#      direct integration fails at every h, homotopy re-init converges and the
#      integration proceeds from it.
#
# Every event is written as a homotopy in s, so the same closure serves the direct
# run (s = 1), Newton (s = 1) and the continuation (s: 0 -> 1):
#   fault application  lambda_fault = s          (geometric x_eff in solve_fault!)
#   fault clearing     lambda_fault = 1 - s
#   line trip          (1-s)*branch equation - s*i = 0   (s = 1: i = 0)
#   generator trip     (1-s)*stator/swing equations + s*i_gen (s = 1: no current, frozen rotor)
#   load step          P, Q scaled by (1 + s*step)
#
# Loads are ZIP (0.7, 0.1, 0.2) with the PSS/E low-voltage characteristic and
# PQBRAK = 0.7. Bus equations (form=):
#   current  repo balance! (C dv/dt = i_net - conj(S/V)), the KCL form
#   power    power mismatch V conj(i_net) - S_load (as in PowerSAS), defined here
# Every state the search accepts (end of a direct run, a re-init solution, end of a
# post-re-init run) must also be physical: KCL current mismatch below KCL_TOL.
# In the power form Newton can converge to V = 0 at a bus (a spurious root);
# such runs count as failures, not successes.
#
# Network forms (mode=):
#   rms  only delta and omega are differential (algebraic lines and bus voltages)
#   emt  line currents differential (inductance), bus voltages algebraic
#   cap  line currents and bus voltages differential (C_eq on the buses)
# The mass matrix and the algebraic index set are built here, independent of
# build_mass_matrix and of the current _setup_alg: the solvers receive an address
# whose keys make _setup_alg return exactly the algebraic rows (asserted).
#
# Outputs: outputs/reinit_search/<mode>_<form>_<run>_<timestamp>/{config.txt, all_cases.csv, saved_cases.csv}

using Barq, MyDiffEq, LinearAlgebra, SparseArrays, Printf, Dates

# ---------------------------------------------------------------- configuration
args = Dict(String(split(a, "="; limit=2)[1]) => String(split(a, "="; limit=2)[2]) for a in ARGS if occursin("=", a))
const DATA     = get(args, "data", "cases/Fault_Cases/ieee39_fault.xlsx")
const MODE     = Symbol(get(args, "mode", "rms"))
const FAMILIES = Symbol.(split(get(args, "families", "line,gen,line2,topo,load,fault,faultline"), ","))
const LF_SPEC  = split(get(args, "lf", "base,0.9,0.97,0.99"), ",")    # "base" = 1.0, numbers = fraction of LF_max
const XF_LIST  = parse.(Float64, split(get(args, "xf", "0.001,0.005,0.02"), ","))
const T_FAULT  = parse(Float64, get(args, "tfault", "0.1"))           # fault duration before clearing
const LOAD_STEPS = parse.(Float64, split(get(args, "loadstep", "0.5,1.0"), ","))
const LIMIT    = parse(Int, get(args, "limit", "0"))                 # max events per loading level (0 = all)
const MATCH    = get(args, "match", "")                              # keep events whose label contains one of "a|b|..."
const METHOD   = Symbol(get(args, "method", "Euler"))
const FORM     = Symbol(get(args, "form", "current"))                # current | power
const KCL_TOL  = parse(Float64, get(args, "kcltol", "1e-3"))         # max KCL mismatch (pu) of an accepted state
const V_EPS    = parse(Float64, get(args, "veps", "1e-12"))          # |V| smoothing in the script-local load model
MODE in (:rms, :emt, :cap) || error("mode must be rms, emt or cap (got $(MODE)); pass each key=value as its own argument")
FORM in (:current, :power) || error("form must be current or power (got $(FORM))")
const LOAD_KW  = (; zip=Tuple(parse.(Float64, split(get(args, "zip", "0.7,0.1,0.2"), ","))), low_voltage=true,
                   pqbrak=parse(Float64, get(args, "pqbrak", "0.7")),
                   characteristic=parse(Int, get(args, "characteristic", "1")))
# step sizes: direct test (n steps each), post-re-init check (h => duration), fault-on integration
# (1e-2, 1e-3, 1e-4)
const H_DIRECT, N_DIRECT, H_POST, H_FAULTON =
    MODE == :rms ? ((1e-3, 1e-4, 1e-5), 20, ((1e-2, 0.3), (1e-3, 0.05)), 1e-2) :
                   ((5e-4, 5e-5, 5e-6), 20, ((5e-4, 0.05), (5e-5, 0.01)), 5e-4)

quiet(f) = redirect_stdout(f, devnull)
# Jacobians are ~200x200: one BLAS thread per process is fastest and lets several
# searches run side by side without oversubscribing the cores.
BLAS.set_num_threads(parse(Int, get(args, "blas", "1")))

# ---------------------------------------------------------------- bus equations
# Net current injected into every bus (generators, fault, lines, charging) and the
# ZIP load power, shared by the power-form balance and the KCL check.
function bus_terms(u, p)
    T = eltype(u); address, models, incidence_matrix, C_eq, nsb, _ = p
    bus, gen, fault, load = models.bus, models.generator, models.fault, models.load
    vd = Vector{T}(bus.vd); vq = Vector{T}(bus.vq)
    vd[nsb] = u[address["balance_d"]]; vq[nsb] = u[address["balance_q"]]
    δ = u[address["delta"]]; gid = u[address["gen_id"]]; giq = u[address["gen_iq"]]
    id = zeros(T, length(bus.idx)); iq = zeros(T, length(bus.idx))
    id[gen.bus] += @. gid*sin(δ) + giq*cos(δ); iq[gen.bus] += @. giq*sin(δ) - gid*cos(δ)
    id[fault.bus] -= u[address["fault_id"]]; iq[fault.bus] -= u[address["fault_iq"]]
    id .+= incidence_matrix*u[address["line_id"]]; iq .+= incidence_matrix*u[address["line_iq"]]
    w = 2pi*60; id .+= @. w*C_eq*vq; iq .-= @. w*C_eq*vd
    k_z, k_i, k_p = LOAD_KW.zip
    pl = zeros(T, length(bus.idx)); ql = zeros(T, length(bus.idx))
    for (k, b) in enumerate(load.bus)
        # |V| smoothed by V_EPS: at an exact V = 0 (reached by underflow on the spurious
        # power-form branch) d|V|/dV is 0/0 and the Jacobian becomes NaN.
        v = sqrt(vd[b]^2 + vq[b]^2 + V_EPS^2); v0 = bus.v[b]
        kP, kI = LOAD_KW.low_voltage ? Barq.Models.BusModel._openipsl_load_factors(v, LOAD_KW.pqbrak, LOAD_KW.characteristic) : (one(v), one(v))
        lf = k_z*(v/v0)^2 + kI*k_i*(v/v0) + kP*k_p
        pl[b] += load.p[k]*lf; ql[b] += load.q[k]*lf
    end
    return vd, vq, id, iq, pl, ql, nsb, address
end
function balance_power!(du, u, p)          # S mismatch: V conj(i_net) - S_load
    vd, vq, id, iq, pl, ql, nsb, a = bus_terms(u, p)
    du[a["balance_d"]] = (@. id*vd + iq*vq - pl)[nsb]
    du[a["balance_q"]] = (@. id*vq - iq*vd - ql)[nsb]
end
function balance_kcl!(du, u, p)            # current mismatch, V = 0 guarded (a spurious root shows as a large mismatch)
    vd, vq, id, iq, pl, ql, nsb, a = bus_terms(u, p)
    v2 = @. max(vd^2 + vq^2, 1e-12)
    du[a["balance_d"]] = (@. id - (pl*vd + ql*vq)/v2)[nsb]
    du[a["balance_q"]] = (@. iq - (pl*vq - ql*vd)/v2)[nsb]
end

# ---------------------------------------------------------------- system at a loading factor
function build_case(lf)
    models = load_data(DATA)
    models.load.p .*= lf; models.load.q .*= lf; models.generator.p_m .*= lf
    sys = build_system(models)
    pf = quiet(() -> solve_power_flow_continuation!(sys; verbose=false))
    pf.converged || return nothing
    quiet(() -> run_static_init!(sys))
    address = build_dynamic_address(sys)
    u0 = build_initial_conditions(sys, address)
    n = length(u0); nsb = sys.non_slack_buses
    M = spzeros(n, n)
    for i in address["delta"]; M[i, i] = 1.0; end
    for (j, i) in enumerate(address["omega"]); M[i, i] = models.generator.M[j]; end
    if MODE in (:emt, :cap)
        for (j, i) in enumerate(address["line_id"]); M[i, i] = models.line.L[j]; end
        for (j, i) in enumerate(address["line_iq"]); M[i, i] = models.line.L[j]; end
    end
    if MODE == :cap
        for (j, i) in enumerate(address["balance_d"]); M[i, i] = sys.C_eq[nsb[j]]; end
        for (j, i) in enumerate(address["balance_q"]); M[i, i] = sys.C_eq[nsb[j]]; end
    end
    # All differential rows after omega form one contiguous block (address order is
    # delta, omega, line_id, line_iq, balance_d, balance_q, fault, gen). Pass it under
    # "line_id" and leave the other keys empty, so any version of _setup_alg that sums
    # these keys returns exactly the algebraic rows.
    n_extra = MODE == :rms ? 0 : MODE == :emt ? 2length(address["line_id"]) : 2length(address["line_id"]) + 2length(nsb)
    w_end = last(address["omega"]); none = 1:0
    alg_addr = Dict("delta" => address["delta"], "omega" => address["omega"], "line_id" => (w_end+1):(w_end+n_extra),
                    "line_iq" => none, "balance_d" => none, "balance_q" => none)
    _, alg_idx = Barq.DynamicSim._setup_alg(u0, alg_addr)
    @assert collect(alg_idx) == findall(iszero, diag(M)) "solver algebraic set does not match the mass matrix"
    pb = (address, sys.models, sys.incidence_matrix, sys.C_eq, nsb)
    return (; lf, models, sys, address, u0, M, alg_addr, alg_idx, pb, nsb,
              p_load0 = copy(models.load.p), q_load0 = copy(models.load.q))
end

# ---------------------------------------------------------------- events
Base.@kwdef struct Event
    family::Symbol
    label::String
    lines::Vector{Int} = Int[]
    gens::Vector{Int} = Int[]
    loads::Vector{Int} = Int[]
    load_step::Float64 = 0.0
    fault_bus::Int = 0
    xf::Float64 = NaN
    fault_mode::Symbol = :none        # :none, :apply, :clear
end

function event_model(c, ev; form=FORM)
    a = c.address; models = c.models
    lid, liq = a["line_id"], a["line_iq"]; gid, giq = a["gen_id"], a["gen_iq"]
    (du, u, p, t) -> begin
        s = p[end]
        λf = ev.fault_mode == :apply ? s : ev.fault_mode == :clear ? 1 - s : 0.0
        pf = (p[1:end-1]..., λf)
        for i in ev.loads
            models.load.p[i] = c.p_load0[i]*(1 + s*ev.load_step); models.load.q[i] = c.q_load0[i]*(1 + s*ev.load_step)
        end
        solve_generator!(du, u, pf); solve_line!(du, u, pf, t); solve_fault!(du, u, pf, t)
        form == :power ? balance_power!(du, u, pf) : form == :kcl ? balance_kcl!(du, u, pf) : balance!(du, u, pf; LOAD_KW...)
        for l in ev.lines
            du[lid[l]] = (1 - s)*du[lid[l]] - s*u[lid[l]]
            du[liq[l]] = (1 - s)*du[liq[l]] - s*u[liq[l]]
        end
        for g in ev.gens
            du[gid[g]] = (1 - s)*du[gid[g]] + s*u[gid[g]]
            du[giq[g]] = (1 - s)*du[giq[g]] + s*u[giq[g]]
            du[a["delta"][g]] *= (1 - s); du[a["omega"][g]] *= (1 - s)
        end
    end
end

function connected_without(c, out_lines)
    nb = length(c.models.bus.idx); adj = [Int[] for _ in 1:nb]
    for l in eachindex(c.models.line.bus1_idx)
        l in out_lines && continue
        b1, b2 = c.models.line.bus1_idx[l], c.models.line.bus2_idx[l]
        push!(adj[b1], b2); push!(adj[b2], b1)
    end
    seen = falses(nb); stack = [setdiff(1:nb, c.nsb)[1]]; seen[stack[1]] = true
    while !isempty(stack)
        b = pop!(stack)
        for nb2 in adj[b]; seen[nb2] || (seen[nb2] = true; push!(stack, nb2)); end
    end
    return all(seen)
end

function make_events(c)
    m = c.models; nl = length(m.line.idx); evs = Event[]
    lines_at(b) = [l for l in 1:nl if m.line.bus1_idx[l] == b || m.line.bus2_idx[l] == b]
    lname(l) = "$(m.line.bus1_idx[l])-$(m.line.bus2_idx[l])"
    if :line in FAMILIES
        for l in 1:nl; connected_without(c, [l]) && push!(evs, Event(family=:line, label="trip $(lname(l))", lines=[l])); end
    end
    if :line2 in FAMILIES     # N-2: pairs of lines sharing a bus
        pairs = Set{Tuple{Int,Int}}()
        for b in m.bus.idx, (i, l1) in enumerate(lines_at(b)), l2 in lines_at(b)[i+1:end]; push!(pairs, (l1, l2)); end
        for (l1, l2) in sort(collect(pairs))
            connected_without(c, [l1, l2]) && push!(evs, Event(family=:line2, label="trip $(lname(l1)) + $(lname(l2))", lines=[l1, l2]))
        end
    end
    if :topo in FAMILIES      # substation-scale change: all lines at a bus except one
        for b in m.bus.idx
            ls = lines_at(b); length(ls) >= 3 || continue
            for keep in ls
                out = setdiff(ls, [keep])
                if connected_without(c, out)
                    push!(evs, Event(family=:topo, label="bus $b: trip $(length(out)) lines, keep $(lname(keep))", lines=out)); break
                end
            end
        end
    end
    if :gen in FAMILIES
        for g in eachindex(m.generator.bus); push!(evs, Event(family=:gen, label="trip gen at bus $(m.generator.bus[g])", gens=[g])); end
    end
    if :load in FAMILIES
        for (i, b) in enumerate(m.load.bus), st in LOAD_STEPS
            abs(c.p_load0[i]) > 1e-6 && push!(evs, Event(family=:load, label="load +$(Int(100st))% at bus $b", loads=[i], load_step=st))
        end
    end
    if :fault in FAMILIES
        for b in c.nsb, xf in XF_LIST
            push!(evs, Event(family=:fault, label="fault bus $b xf=$xf apply", fault_bus=b, xf=xf, fault_mode=:apply))
            push!(evs, Event(family=:fault, label="fault bus $b xf=$xf clear@$(T_FAULT)", fault_bus=b, xf=xf, fault_mode=:clear))
        end
    end
    if :faultline in FAMILIES   # fault at the sending end of a line, cleared by tripping that line
        for l in 1:nl, xf in XF_LIST[1:min(2, end)]
            b = m.line.bus1_idx[l]
            (b in c.nsb && connected_without(c, [l])) || continue
            push!(evs, Event(family=:faultline, label="fault bus $b xf=$xf cleared by trip $(lname(l))", lines=[l], fault_bus=b, xf=xf, fault_mode=:clear))
        end
    end
    pats = split(MATCH, "|"); evs = filter(e -> any(p -> occursin(p, e.label), pats), evs)
    return LIMIT > 0 ? evs[1:min(LIMIT, end)] : evs
end

# ---------------------------------------------------------------- numerics
function integrate(c, f, u, s, h, T)
    try
        sol = quiet(() -> MyDiffEq.Solve(MyDiffEq.ODEProblem(f, u, (0.0, T), (c.pb..., s), c.M), h,
                                         method=METHOD, adaptive=false, tstops=[], always_new=true))
        return sol.retcode == :Success, sol.retcode, sol.time[end], sol
    catch e
        return false, :exception, NaN, nothing
    end
end
alg_resid(c, f, u, s) = (du = zeros(length(u)); f(du, u, (c.pb..., s), 0.0); norm(du[c.alg_idx], Inf))
function kcl_mismatch(c, fk, u)             # fk = event_model(...; form=:kcl), evaluated post-event
    du = zeros(length(u)); fk(du, u, (c.pb..., 1.0), 0.0)
    r = du[vcat(collect(c.address["balance_d"]), collect(c.address["balance_q"]))]
    return all(isfinite, r) ? norm(r, Inf) : Inf
end
vmag(c, u) = hypot.(u[c.address["balance_d"]], u[c.address["balance_q"]])

function newton_reinit(c, f, u)
    ur = copy(u)
    r = try quiet(() -> solve_newton!(ur, (c.pb..., 1.0), c.alg_addr; max_iter=50, always_new=true, model! = f)) catch; nothing end
    ok = r !== nothing && r.converged && all(isfinite, ur)
    return ok, ur, r === nothing ? -1 : r.iters
end
function homotopy_reinit(c, f, u)
    ur = copy(u)
    # r = try quiet(() -> solve_homotopy!(ur, c.pb, c.alg_addr; Δλ=0.01, always_new=true, model! = f)) catch; nothing end
    r = try quiet(() -> solve_adaptive_homotopy!(ur, c.pb, c.alg_addr; always_new=true, model! = f)) catch; nothing end
    ok = r !== nothing && r.converged && all(isfinite, ur)
    return ok, ur, r === nothing ? NaN : something(r.λ_failed, 1.0), r === nothing ? -1 : r.total_iters
end
function post_check(c, f, u, failed_h, phys)
    # Recovery must hold at every step size where the direct run failed (same
    # horizon), plus the longer checks in H_POST.
    checks = vcat([(h, N_DIRECT*h) for h in failed_h], collect(H_POST))
    oks = Bool[]; notes = String[]
    for (h, T) in checks
        ok, rc, tf, sol = integrate(c, f, u, 1.0, h, T)
        k = ok ? phys(sol.u[end]) : NaN
        ok &= k < KCL_TOL
        push!(oks, ok); push!(notes, ok ? @sprintf("h=%g/%g ok", h, T) :
                              isfinite(tf) && rc == :Success ? @sprintf("h=%g/%g spurious(kcl=%.1e)", h, T, k) : @sprintf("h=%g/%g %s@%.4g", h, T, rc, tf))
    end
    return all(oks), any(oks), join(notes, "; ")
end

# ---------------------------------------------------------------- one case
const HEADER = "lf,lf_abs,family,label,lines,gens,loads,load_step,fault_bus,xf,fault_mode,pre_resid," *
               join(["direct_h$(h)" for h in H_DIRECT], ",") * ",direct_all_ok,direct_all_fail,direct_spurious," *
               "newton_ok,newton_it,newton_resid,newton_kcl,newton_minV,post_newton," *
               "hom_ok,hom_lambda_end,hom_it,hom_resid,hom_kcl,hom_minV,post_hom,root_dist,class,saved,target,seconds"
csvq(x) = "\"" * replace(string(x), "\"" => "'") * "\""

function run_case(c, ev, lf_tag)
    t0 = time()
    try
        c.models.fault.bus[1] = ev.fault_bus > 0 ? ev.fault_bus : c.nsb[1]
        c.models.fault.x_fault[1] = isnan(ev.xf) ? 0.01 : ev.xf
        f = event_model(c, ev)
        fk = event_model(c, ev; form=:kcl); phys(u) = kcl_mismatch(c, fk, u)
        u_inh = copy(c.u0)
        if ev.fault_mode == :clear          # fault-on trajectory: homotopy application, then integrate
            fa = event_model(c, Event(family=:fault, label="", fault_bus=ev.fault_bus, xf=ev.xf, fault_mode=:apply))
            ok, ua, _, _ = homotopy_reinit(c, fa, c.u0)
            ok || return (; row = "", note = "fault application re-init failed")
            ok, _, _, sol = integrate(c, fa, ua, 1.0, H_FAULTON, T_FAULT)
            ok || return (; row = "", note = "fault-on integration failed")
            u_inh = copy(sol.u[end])
            fka = event_model(c, Event(family=:fault, label="", fault_bus=ev.fault_bus, xf=ev.xf, fault_mode=:apply); form=:kcl)
            kcl_mismatch(c, fka, u_inh) < KCL_TOL || return (; row = "", note = "fault-on state not physical")
        end
        pre = alg_resid(c, f, u_inh, 0.0)
        direct = [integrate(c, f, u_inh, 1.0, h, N_DIRECT*h) for h in H_DIRECT]
        dk = [d[1] ? phys(d[4].u[end]) : NaN for d in direct]          # KCL mismatch at the end of each run
        dphys = [d[1] && k < KCL_TOL for (d, k) in zip(direct, dk)]
        dstr = [p ? "ok" : d[1] ? @sprintf("spurious(kcl=%.1e)", k) : @sprintf("%s@%.2e", d[2], d[3]) for (d, k, p) in zip(direct, dk, dphys)]
        all_ok = all(dphys); all_fail = !any(dphys)
        n_spur = count(d[1] && !p for (d, p) in zip(direct, dphys))
        failed_h = [h for (h, p) in zip(H_DIRECT, dphys) if !p]
        nk, nit, nres, nkcl, nmin, pn = false, 0, NaN, NaN, NaN, ""
        hk, hl, hit, hres, hkcl, hmin, ph, dist = false, NaN, 0, NaN, NaN, NaN, "", NaN
        post_n_ok = post_h_ok = false
        if !all_ok
            nconv, un, nit = newton_reinit(c, f, u_inh)
            if nconv
                nres = alg_resid(c, f, un, 1.0); nkcl = phys(un); nmin = minimum(vmag(c, un))
                nk = nkcl < KCL_TOL                          # converged AND physical
                nk && ((post_n_ok, _, pn) = post_check(c, f, un, failed_h, phys))
                nk || (pn = @sprintf("converged to non-physical root (kcl=%.1e)", nkcl))
            end
            hconv, uh, hl, hit = homotopy_reinit(c, f, u_inh)
            if hconv
                hres = alg_resid(c, f, uh, 1.0); hkcl = phys(uh); hmin = minimum(vmag(c, uh))
                hk = hkcl < KCL_TOL
                hk && ((post_h_ok, _, ph) = post_check(c, f, uh, failed_h, phys))
                hk || (ph = @sprintf("converged to non-physical root (kcl=%.1e)", hkcl))
                nconv && (dist = norm(un[c.alg_idx] - uh[c.alg_idx], Inf))
            end
        end
        cls = all_ok ? "A_direct_ok" :
              (!nk && !hk) ? "C_no_reinit" :
              (hk && post_h_ok && (!nk || dist > 1e-4)) ? "B_homotopy_only" :
              ((nk && post_n_ok) || (hk && post_h_ok)) ? "reinit_recovers" : "reinit_then_fails"
        saved = !all_ok && ((nk && post_n_ok) || (hk && post_h_ok))
        target = all_fail && hk && post_h_ok
        row = join([lf_tag, @sprintf("%.4f", c.lf), ev.family, csvq(ev.label), csvq(ev.lines), csvq(ev.gens), csvq(ev.loads),
                    ev.load_step, ev.fault_bus, ev.xf, ev.fault_mode, @sprintf("%.2e", pre), csvq.(dstr)..., all_ok, all_fail, n_spur,
                    nk, nit, @sprintf("%.2e", nres), @sprintf("%.2e", nkcl), @sprintf("%.3f", nmin), csvq(pn),
                    hk, hl, hit, @sprintf("%.2e", hres), @sprintf("%.2e", hkcl), @sprintf("%.3f", hmin), csvq(ph), @sprintf("%.2e", dist),
                    cls, saved, target, @sprintf("%.1f", time() - t0)], ",")
        return (; row, note = cls, saved, target)
    finally
        c.models.load.p .= c.p_load0; c.models.load.q .= c.q_load0
    end
end

# ---------------------------------------------------------------- loading levels
function lf_max_search(lo=1.0, hi=4.0, tol=0.005)
    build_case(lo) === nothing && error("base case power flow failed")
    while build_case(hi) !== nothing; lo, hi = hi, 2hi; end
    while hi - lo > tol
        mid = (lo + hi)/2
        build_case(mid) === nothing ? (hi = mid) : (lo = mid)
    end
    return lo
end

function main()
    run = get(args, "run", "")
    outdir = joinpath("outputs", "reinit_search", join(filter(!isempty, [string(MODE), string(FORM), run, Dates.format(now(), "yyyymmdd_HHMMSS")]), "_"))
    mkpath(outdir)
    lfmax = any(!=("base"), LF_SPEC) ? lf_max_search() : NaN
    levels = [(tag, tag == "base" ? 1.0 : parse(Float64, tag)*lfmax) for tag in LF_SPEC]
    open(joinpath(outdir, "config.txt"), "w") do io
        println(io, "data=$DATA mode=$MODE form=$FORM kcl_tol=$KCL_TOL method=$METHOD families=$FAMILIES")
        println(io, "loads=$LOAD_KW xf=$XF_LIST tfault=$T_FAULT loadsteps=$LOAD_STEPS")
        println(io, "h_direct=$H_DIRECT x $N_DIRECT steps, h_post=$H_POST, h_faulton=$H_FAULTON")
        println(io, "LF_max (power-flow continuation) = $lfmax; levels = $levels")
    end
    all_io = open(joinpath(outdir, "all_cases.csv"), "w"); println(all_io, HEADER)
    sav_io = open(joinpath(outdir, "saved_cases.csv"), "w"); println(sav_io, HEADER); flush(all_io); flush(sav_io)
    @printf "mode=%s LF_max=%.4f levels=%s -> %s\n" MODE lfmax string(levels) outdir
    for (tag, lf) in levels
        c = build_case(lf)
        c === nothing && (println("LF $tag ($lf): power flow failed, skipped"); continue)
        evs = make_events(c)
        @printf "LF %s = %.4f: %d events, pre-event residual %.1e\n" tag lf length(evs) alg_resid(c, event_model(c, Event(family=:none, label="")), c.u0, 0.0)
        counts = Dict{String,Int}()
        for (k, ev) in enumerate(evs)
            r = run_case(c, ev, tag)
            counts[r.note] = get(counts, r.note, 0) + 1
            if r.row != ""
                println(all_io, r.row); flush(all_io)
                r.saved && (println(sav_io, r.row); flush(sav_io))
                (r.saved || r.target) && @printf "  [%s] %s: %s%s\n" tag ev.label r.note (r.target ? "  <-- TARGET" : "")
            end
            k % 25 == 0 && (@printf "  %d/%d done %s\n" k length(evs) string(counts); flush(stdout))
        end
        @printf "LF %s summary: %s\n" tag string(counts); flush(stdout)
    end
    close(all_io); close(sav_io)
    println("done: $outdir")
end

if abspath(PROGRAM_FILE) == @__FILE__   # include() for diagnostics without running the search
    main()
end
