# Run from the repository root:
#   julia --project=. script/plot_converter_case1.jl [output_directory]
# Produces PDF/PNG figures and numerical records. No existing script is executed
# or modified: only GFM!, GFL!, and the parameter tuple p are loaded from it.
ENV["GKSwstype"] = "100"
using ForwardDiff, LinearAlgebra, Printf, MyDiffEq, Plots

const ROOT = normpath(joinpath(@__DIR__, ".."))
const SOURCE = joinpath(ROOT, "scripts", "converter_homotopy.jl")
const OUT = isempty(ARGS) ? joinpath(ROOT, "figures", "review1", "converter_case1") : abspath(ARGS[1])
const STEPS = (1e-4, 1e-5, 1e-6)
const METHODS = (:Euler, :Trap)
const MAXITER = 80
const ALG_TOL = 1e-10
const END_TIME = 0.01
const LAMBDAS = collect(range(0.0, 1.0; length=6))
const COLORS = ["#0072B2", "#D55E00", "#009E73"]
const TEAL = "#008C95"
const ORANGE = "#D55E00"

module ConverterEquations end
let source = read(SOURCE, String), position = 1, found = Set{Symbol}()
    while position <= lastindex(source)
        expression, position = Meta.parse(source, position; raise=true)
        expression === nothing && break
        if expression isa Expr && expression.head == :function
            signature = expression.args[1]
            if signature isa Expr && signature.head == :call && signature.args[1] in (:GFM!, :GFL!)
                Core.eval(ConverterEquations, expression)
                push!(found, signature.args[1])
            end
        elseif expression isa Expr && expression.head == :(=) && expression.args[1] == :p
            Core.eval(ConverterEquations, expression)
            push!(found, :p)
        end
    end
    @assert found == Set([:GFM!, :GFL!, :p]) "Converter definitions or parameters were not found."
end
const P = ConverterEquations.p
with_lambda(p, lambda) = (p[1:9]..., lambda*p[10], p[11])
rhs(u, p) = (r = similar(u); ConverterEquations.GFL!(r, u, p, 0.0); r)
algebraic(y, p, eq) = rhs(vcat(eq, y), p)[2:end]
# Fixed row scaling removes omega/X from the two line equations. All algebraic
# variables use unit per-unit scales. Scaling is diagnostic only, not a solver change.
const ROW_SCALE = [P[4]/P[5], P[4]/P[5], 1.0, 1.0, 1.0, 1.0]

function newton_trace(f, initial; maxiter=MAXITER, tol=ALG_TOL)
    y = copy(initial)
    records = NamedTuple[]
    for k in 0:maxiter
        r = f(y)
        if !all(isfinite, r) || !all(isfinite, y)
            return (; y, records, success=false, updates=k, status="nonfinite")
        end
        J = ForwardDiff.jacobian(f, y)
        s = svdvals(Diagonal(ROW_SCALE)*J)
        entry = (; iteration=k, residual=norm(r), scaled_residual=norm(ROW_SCALE.*r, Inf),
                 sigma_min=minimum(s), sigma_max=maximum(s), condition=maximum(s)/minimum(s),
                 raw_condition=cond(J), y=copy(y))
        push!(records, entry)
        norm(r, Inf) <= tol && return (; y, records, success=true, updates=k, status="converged")
        k == maxiter && break
        delta = try
            -(J\r)
        catch error
            error isa SingularException || rethrow()
            return (; y, records, success=false, updates=k, status="singular")
        end
        y += delta
    end
    return (; y, records, success=false, updates=maxiter, status="iteration_limit")
end

function reduced_constraints(vd, vq, p, eq)
    vr, E, R, X, omega, tq, ki, kp, iL, gf, imax = p
    V = hypot(vd, vq)
    ild = (R*(vd-E) + X*vq)/(R^2+X^2)
    ilq = (-X*(vd-E) + R*vq)/(R^2+X^2)
    a = iL/max(V, 0.7) + gf
    icd, icq = ild+a*vd, ilq+a*vq
    return icq-eq-kp*(vr-V), icd^2+icq^2-imax^2
end

function step_residual(u, previous, h, method)
    f = rhs(u, P)
    r = copy(f)
    r[1] = u[1]-previous[1]-h*(method == :Euler ? f[1] : (f[1]+rhs(previous, P)[1])/2)
    return r
end

function simulate(initial, h, method, tend)
    mass = zeros(7, 7)
    mass[1, 1] = 1.0
    problem = MyDiffEq.ODEProblem(ConverterEquations.GFL!, copy(initial), (0.0, tend), P, mass)
    return MyDiffEq.Solve(problem, h; method, adaptive=false, always_new=true,
                         max_iter=MAXITER, tol=1e-9)
end

function savefigure(p, name; width=1100, height=420)
    plot!(p; size=(width, height))
    savefig(p, joinpath(OUT, name*".pdf"))
    savefig(p, joinpath(OUT, name*".png"))
end

function write_records(path, rows)
    isempty(rows) && return
    open(path, "w") do io
        println(io, join(string.(keys(first(rows))), ","))
        for row in rows
            println(io, join(values(row), ","))
        end
    end
end

function geometry_figure(pre, direct, continuation)
    eq = pre[1]
    vdroots = [s.y[3] for s in continuation]
    vqroots = [s.y[4] for s in continuation]
    xrange = extrema(vcat(pre[4], vdroots))
    yrange = extrema(vcat(pre[5], vqroots))
    padx = max(0.12, 0.25*(xrange[2]-xrange[1]))
    pady = max(0.12, 0.25*(yrange[2]-yrange[1]))
    xs = range(xrange[1]-padx, xrange[2]+padx; length=450)
    ys = range(yrange[1]-pady, yrange[2]+pady; length=450)
    panels = Any[]
    for homotopy in (false, true)
        p = plot(; xlabel="v_d (p.u.)", ylabel="v_q (p.u.)", aspect_ratio=:equal,
                 xlims=extrema(xs), ylims=extrema(ys), legend=:bottomleft, legendfontsize=7,
                 title=homotopy ? "(b) Continuation solutions" : "(a) Direct Newton: voltage projection")
        stages = homotopy ? [1, 3, length(LAMBDAS)] : [length(LAMBDAS)]
        for j in stages
            lambda = LAMBDAS[j]
            for (constraint, color) in ((1, COLORS[1]), (2, ORANGE))
                z = [reduced_constraints(x, y, with_lambda(P, lambda), eq)[constraint] for y in ys, x in xs]
                contour!(p, xs, ys, z; levels=[0.0], color, colorbar=false,
                         linewidth=lambda == 1 ? 2.2 : 1.2, linestyle=lambda == 1 ? :solid : :dash,
                         alpha=lambda == 1 ? 1.0 : 0.4, label="")
            end
        end
        plot!(p, [NaN], [NaN]; color=COLORS[1], label="Reactive-current constraint")
        plot!(p, [NaN], [NaN]; color=ORANGE, label="Current-circle constraint")
        scatter!(p, [pre[4]], [pre[5]]; marker=:diamond, color=:black, markersize=6, label="Inherited pre-event point")
        scatter!(p, [last(vdroots)], [last(vqroots)]; marker=:star5, color=TEAL, markersize=9, label="Target solution")
        if homotopy
            plot!(p, vdroots, vqroots; color=TEAL, marker=:circle, markersize=4, arrow=true, label="Converged stages")
            for j in eachindex(LAMBDAS)
                annotate!(p, vdroots[j]+0.012, vqroots[j]+0.017,
                          text(@sprintf("%.1f", LAMBDAS[j]), 8, TEAL))
            end
            plot!(p, [pre[4], first(vdroots)], [pre[5], first(vqroots)];
                  color=:gray, linestyle=:dot, label="Auxiliary initialization (endpoints)")
        else
            if length(direct.records) > 1
                start = direct.records[1].y[3:4]
                next = direct.records[2].y[3:4]
                direction = next-start
                limits = (extrema(xs), extrema(ys))
                fraction = 1.0
                for axis in 1:2
                    if direction[axis] > 0
                        fraction = min(fraction, (limits[axis][2]-start[axis])/direction[axis])
                    elseif direction[axis] < 0
                        fraction = min(fraction, (limits[axis][1]-start[axis])/direction[axis])
                    end
                end
                endpoint = start + (fraction < 1 ? 0.94fraction : fraction)*direction
                plot!(p, [start[1], endpoint[1]], [start[2], endpoint[2]];
                      color="#CC3311", arrow=true, label="Direction of first Newton update")
                if fraction < 1
                    annotate!(p, sum(extrema(xs))/2, last(ys)-0.025,
                              text("First iterate lies outside this view", 8, "#CC3311"))
                end
            end
        end
        push!(panels, p)
    end
    savefigure(plot(panels...; layout=(1,2)), "01_constraint_geometry"; height=530)
end

function main()
    mkpath(OUT)
    gr()
    default(; fontfamily="Times", linewidth=1.8, markersize=3, guidefontsize=11,
            tickfontsize=9, legendfontsize=8, titlefontsize=11, framestyle=:box,
            gridalpha=0.15, dpi=220, margin=4Plots.mm)
    # Analytic pre-event equilibrium of the source model; checked against all 7 equations.
    vr, E, R, X, omega, tq, ki, kp, iL, gf, imax = P
    ild = R*(vr-E)/(R^2+X^2)
    ilq = -X*(vr-E)/(R^2+X^2)
    pre = [0.0, ild, ilq, vr, 0.0, ild+iL*vr/max(vr,0.7), ilq]
    check = zeros(7)
    ConverterEquations.GFM!(check, pre, P, 0.0)
    @assert norm(check, Inf) < 1e-12
    eq = pre[1]
    initial = pre[2:end]
    direct = newton_trace(y -> algebraic(y, P, eq), initial)
    continuation = []
    y = copy(initial)
    traces = NamedTuple[]
    for lambda in LAMBDAS
        stage = newton_trace(z -> algebraic(z, with_lambda(P, lambda), eq), y)
        @assert stage.success "Continuation stage lambda=$lambda failed; inspect the equations/parameters."
        push!(continuation, stage)
        y = copy(stage.y)
        for rec in stage.records
            push!(traces, (; run="homotopy", lambda, iteration=rec.iteration,
                           residual_2=rec.residual, scaled_residual_inf=rec.scaled_residual,
                           sigma_min=rec.sigma_min, sigma_max=rec.sigma_max,
                           scaled_condition=rec.condition, raw_condition=rec.raw_condition,
                           ild=rec.y[1], ilq=rec.y[2], vd=rec.y[3], vq=rec.y[4], icd=rec.y[5], icq=rec.y[6]))
        end
    end
    for rec in direct.records
        push!(traces, (; run="direct", lambda=1.0, iteration=rec.iteration,
                       residual_2=rec.residual, scaled_residual_inf=rec.scaled_residual,
                       sigma_min=rec.sigma_min, sigma_max=rec.sigma_max,
                       scaled_condition=rec.condition, raw_condition=rec.raw_condition,
                       ild=rec.y[1], ilq=rec.y[2], vd=rec.y[3], vq=rec.y[4], icd=rec.y[5], icq=rec.y[6]))
    end
    restart = vcat(eq, y)
    @assert restart[1] == pre[1]
    @assert norm(algebraic(y, P, eq), Inf) <= ALG_TOL
    @assert norm(collect(reduced_constraints(y[3], y[4], P, eq)), Inf) < 1e-9
    geometry_figure(pre, direct, continuation)

    failure_plot = plot(; xlabel="Newton iteration", ylabel="Step residual (2-norm)", yscale=:log10,
                        title="(a) First step without reinitialization", legend=:bottomright)
    restart_plot = plot(; xlabel="Newton iteration", ylabel="Step residual (2-norm)", yscale=:log10,
                        title="First step after reinitialization", legend=:topright, xticks=0:5)
    voltage_plot = plot(; xlabel="Time after event (ms)", ylabel="PCC voltage (p.u.)",
                        title="(d) Resumed DAE integration", legend=:topright)
    plot!(voltage_plot, [-0.5, 0.0], [hypot(pre[4],pre[5]), hypot(pre[4],pre[5])]; color=:black, label="Pre-event")
    plot!(voltage_plot, [0.0, 0.0], [hypot(pre[4],pre[5]), hypot(restart[4],restart[5])]; color=:gray, linestyle=:dot, label="Event jump")
    summaries, steptraces, trajectories = NamedTuple[], NamedTuple[], NamedTuple[]
    for method in METHODS, (j,h) in enumerate(STEPS)
        style = method == :Euler ? :solid : :dash
        label = "$(method == :Euler ? "Euler" : "Trap"), $(round(Int,h*1e6)) μs"
        first_direct = simulate(pre, h, method, h)
        first_restart = simulate(restart, h, method, h)
        full = simulate(restart, h, method, END_TIME)
        @assert first_restart.retcode == :Success && full.retcode == :Success
        @assert isapprox(full.time[end], END_TIME; atol=h*1e-6)
        maxalg = maximum(norm(rhs(u, P)[2:end], Inf) for u in full.u)
        maxstep = maximum(norm(step_residual(full.u[k], full.u[k-1], full.time[k]-full.time[k-1], method), Inf) for k in 2:length(full.u))
        @assert maxalg < 1e-7 && maxstep < 1e-7
        for (kind, sol, plt, previous) in (("direct",first_direct,failure_plot,pre), ("restart",first_restart,restart_plot,restart))
            rr = sol.newton_log.residual_norm
            candidate = sol.retcode == :Success ? sol.u[end] : sol.newton_log.u_final
            final_residual = step_residual(candidate,previous,h,method)
            plot!(plt, 0:length(rr), max.(vcat(rr,norm(final_residual)),1e-16);
                  color=COLORS[j], linestyle=style, label)
            push!(summaries, (; run=kind, method=string(method), step_s=h, status=string(sol.retcode),
                               iterations=sol.newton_log.iters, final_step_residual_inf=norm(step_residual(candidate,previous,h,method),Inf),
                               max_trajectory_algebraic_residual=kind=="restart" ? maxalg : NaN,
                               max_trajectory_step_residual=kind=="restart" ? maxstep : NaN,
                               final_voltage=kind=="restart" ? hypot(full.u[end][4],full.u[end][5]) : NaN))
            for (k,r) in enumerate(rr)
                push!(steptraces, (; run=kind, method=string(method), step_s=h, iteration=k-1, residual_2=r,
                                    correction_2=sol.newton_log.correction_norm[k]))
            end
            push!(steptraces, (; run=kind, method=string(method), step_s=h,
                                iteration=length(rr), residual_2=norm(final_residual), correction_2=NaN))
        end
        plot!(voltage_plot, 1e3 .* full.time, [hypot(u[4],u[5]) for u in full.u]; color=COLORS[j], linestyle=style, label)
        for (t,u) in zip(full.time, full.u)
            push!(trajectories, (; method=string(method), step_s=h, time_s=t, eq=u[1], voltage=hypot(u[4],u[5]),
                                vd=u[4], vq=u[5], current=hypot(u[6],u[7])))
        end
    end

    algebraic_plot = plot([r.iteration for r in direct.records], max.([r.residual for r in direct.records],1e-16);
                         color=ORANGE, label="Direct Newton ($(direct.status))", yscale=:log10,
                         xlabel="Cumulative Newton updates", ylabel="Algebraic residual (2-norm)",
                         title="(b) Algebraic reinitialization", legend=:bottomright)
    offset = 0
    for (j,stage) in enumerate(continuation)
        xx = offset .+ [r.iteration for r in stage.records]
        plot!(algebraic_plot, xx, max.([r.residual for r in stage.records],1e-16); color=TEAL,
              label=j==1 ? "Homotopy: stage residual" : "")
        vline!(algebraic_plot, [offset]; color=:gray, linewidth=0.6, linestyle=:dot, label="")
        annotate!(algebraic_plot, offset+0.4, isodd(j) ? 30.0 : 2000.0,
                  text(@sprintf("%.1f", LAMBDAS[j]), 7, TEAL))
        offset += stage.updates
    end
    effort_plot = bar(LAMBDAS, [s.updates for s in continuation]; color=TEAL, bar_width=0.11, label="",
                      xlabel="Continuation parameter lambda", ylabel="Newton updates",
                      title="(c) Continuation effort (total = $offset)", xticks=LAMBDAS,
                      ylims=(0, maximum(s.updates for s in continuation)+2))
    annotate!(effort_plot, 0.0, continuation[1].updates+0.6, text("Auxiliary solve",8))
    savefigure(plot(failure_plot, algebraic_plot, effort_plot, voltage_plot; layout=(2,2)), "02_numerical_validation"; height=850)
    savefigure(restart_plot, "03_restart_newton"; width=660, height=440)
    savefigure(voltage_plot, "04_voltage_response"; width=660, height=440)
    diagnostic = plot(LAMBDAS, [last(s.records).sigma_min for s in continuation]; color=TEAL, marker=:circle,
                      xlabel="Continuation parameter lambda", ylabel="Scaled minimum singular value", label="",
                      title="Converged continuation points",
                      ylims=(0, 1.1maximum(last(s.records).sigma_min for s in continuation)))
    savefigure(diagnostic, "05_jacobian_diagnostic"; width=660, height=380)
    write_records(joinpath(OUT,"algebraic_newton_records.csv"), traces)
    write_records(joinpath(OUT,"integration_summary.csv"), summaries)
    write_records(joinpath(OUT,"step_newton_records.csv"), steptraces)
    write_records(joinpath(OUT,"voltage_trajectories.csv"), trajectories)
    write_records(joinpath(OUT,"continuation_effort.csv"), [(; lambda=LAMBDAS[j], updates=s.updates,
                  final_residual_inf=norm(algebraic(s.y,with_lambda(P,LAMBDAS[j]),eq),Inf)) for (j,s) in enumerate(continuation)])
    open(joinpath(OUT,"README.md"),"w") do io
        println(io, """
        # Converter case-study figures

        Generated by `script/plot_converter_case1.jl` from GFM!, GFL!, and p in
        `scripts/converter_homotopy.jl`. Parameters: $P.

        Only e_q is differential. Reinitialization preserves e_q = $eq exactly.
        The event imposes the current-limited mode and adds g_f; this script does
        not establish that the mode change is triggered automatically by a limiter.
        The load is constant current with impedance conversion below 0.7 p.u.

        Continuation ramps g_f while the current-limited equations are active.
        Lambda=0 therefore requires an auxiliary solve from the inherited state;
        it is not the pre-event manifold. Total Newton updates: $offset, including
        $(continuation[1].updates) for the auxiliary solve.
        Direct algebraic Newton: $(direct.status), $(direct.updates) updates.

        Algebraic Newton checks raw infinity-norm residual <= $ALG_TOL and allows
        $MAXITER updates. MyDiffEq uses correction tolerance 1e-9 and $MAXITER iterations.
        All restarted trajectories are independently checked for step and algebraic
        residuals below 1e-7 through $END_TIME seconds. Trapezoidal is restarted with
        the post-event vector field. CSV files retain actual statuses and diagnostics.
        Step residual plots include an independently evaluated final residual after
        the last Newton update; its CSV correction is NaN because no further step is taken.

        Geometry uses exact elimination of the algebraic line and KCL equations.
        The red arrow is the voltage projection of the first full six-variable Newton
        update, not a reduced Newton step. If that iterate is outside the root-focused
        viewport, the arrow is shortened and annotated; complete iterates are in CSV.
        Homotopy labels are lambda values. Gray dotted connection shows only the
        inherited/auxiliary endpoints, not a computed solution branch.
        Dashed contours show lambda=0 and 0.4; solid contours show lambda=1.

        In panel (b), each teal segment uses its own H(y,lambda); boundaries reset
        the residual, and separate segments are deliberately not connected.
        Log plots clip values below 1e-16 for display only. Diagnostic row scaling
        is diag(X/omega,X/omega,1,1,1,1), with unit variable scales, and never changes
        the Newton equations. Singular values are sampled, not a certified path bound.
        Timings and physical validation of the converter controller are outside scope.
        """)
    end
    println("Figures and records written to $OUT")
end

main()
