using Pkg; Pkg.activate(".")

using Barq

using MyDiffEq

using JSON
using LinearAlgebra

const OMEGA = 2 * pi * 60

"""Open `lines` in place: scale R and X by `mult`, keeping L = X/OMEGA and dropping the charging B."""
function _open_lines!(models, lines::AbstractVector; mult=1e4)
    models.line.R[lines] .*= mult
    models.line.X[lines] .*= mult
    models.line.L[lines] .= models.line.X[lines] ./ OMEGA
    models.line.B[lines] .= 0.0
    models.line.C[lines] .= 0.0
    return nothing
end

"""Build the ODE problem for one simulation phase."""
function _phase_problem(u0, address, mass_matrix, sys, lambda, t_end)
    p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
    return MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, t_end), p, mass_matrix)
end

"""Solve `prob` non-adaptively, returning the solution or the thrown exception."""
function _solve(prob, dt; method=:Euler, max_iter=50)
    try
        return MyDiffEq.Solve(prob, dt,
                              method=method,
                              adaptive=false,
                              tstops=[],
                              always_new=true,
                              max_iter=max_iter)
    catch e
        return e
    end
end

"""L2 norm of the algebraic constraint residual of `prob` evaluated at state `u`."""
function _algebraic_residual(prob, u)
    du = zeros(length(u))
    prob.f!(du, u, prob.p, 0.0)
    return norm(du[prob.nx+1:end], 2)
end

"""Return (retcode, first_step_failure) for a phase result, using the exception name if it threw."""
function _classify(result)
    result isa Exception && return (nameof(typeof(result)), false)
    return (result.retcode, result.retcode == :MaxIter && length(result.time) == 1)
end

"""Trip one line as a breaker would and report whether Newton fails on the first post-trip step."""
function _run_case(data_file, line; dt, method, max_iter, mult, settle_steps, post_steps)
    models = load_data(data_file)
    # models.load.p[:] .*= 1.2
    sys = build_system(models)
    solve_power_flow!(sys)
    run_static_init!(sys)
    address = build_dynamic_address(sys)
    mass_pre = build_mass_matrix(sys, address)
    u0 = build_initial_conditions(sys, address)

    # lambda = 0 holds the fault branch open, so this phase is a pure pre-trip settle
    prob_pre = _phase_problem(u0, address, mass_pre, sys, 0.0, settle_steps * dt)
    sol_pre = _solve(prob_pre, dt; method=method, max_iter=max_iter)
    if sol_pre isa Exception || sol_pre.retcode != :Success
        return nothing
    end
    u_pre = sol_pre.u[end]

    id_pre = u_pre[address["line_id"][line]]
    iq_pre = u_pre[address["line_iq"][line]]

    # breaker action: open the branch, then interrupt its current
    _open_lines!(models, [line]; mult=mult)
    sys_post = build_system(models)
    mass_post = build_mass_matrix(sys_post, address)
    u_post = copy(u_pre)
    u_post[address["line_id"][line]] = 0.0
    u_post[address["line_iq"][line]] = 0.0

    prob_post = _phase_problem(u_post, address, mass_post, sys_post, 0.0, post_steps * dt)
    resid0 = _algebraic_residual(prob_post, u_post)
    sol_post = _solve(prob_post, dt; method=method, max_iter=max_iter)
    retcode, first_step = _classify(sol_post)

    threw = sol_post isa Exception
    return Dict("line" => line,
                "bus1" => Int(models.line.bus1_idx[line]),
                "bus2" => Int(models.line.bus2_idx[line]),
                "i_pre" => sqrt(id_pre^2 + iq_pre^2),
                "resid0" => resid0,
                "retcode" => String(retcode),
                "first_step_failure" => first_step,
                # an exception aborts the solve, so the step it happened on is unknown
                "steps_completed" => threw ? nothing : length(sol_post.time) - 1,
                "iters" => threw ? Int[] : copy(sol_post.iters),
                "correction_norm" => threw ? Float64[] : copy(sol_post.newton_log.correction_norm),
                "exception" => threw ? sprint(showerror, sol_post) : nothing)
end

# data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
# data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
# data_file = "cases/Simple_Cases/wecc_full_slack.xlsx"

"""Sweep every single-line outage (N-1) and record the trips whose first step fails to converge."""
function main(; data_file="cases/Fault_Cases/ieee39_fault.xlsx",
                # on ieee39 the first post-trip step converges for every line up to dt = 1e-3;
                # first-step MaxIter failures start appearing at dt = 2e-3
                dt=5e-4,
                method=:Euler,
                # max_iter is what separates "hard" from "failed" -- a healthy first step
                # takes 4-6 iterations here, so anything near the cap is diverging
                max_iter=50,
                mult=1e4,
                settle_steps=200,
                post_steps=3,
                out_dir="results_topology")

    line_idx = load_data(data_file).line.idx
    records = Dict{String,Any}[]
    skipped = Int[]

    for line in eachindex(line_idx)
        println("================================================")
        println("Tripping line $line of $(length(line_idx))")
        println("================================================")
        record = _run_case(data_file, line;
                           dt=dt, method=method, max_iter=max_iter,
                           mult=mult, settle_steps=settle_steps, post_steps=post_steps)
        if record === nothing
            println("pre-trip settle failed -- line $line not attributable, skipping")
            push!(skipped, line)
            continue
        end
        push!(records, record)
        flag = record["first_step_failure"] ? " <<< first-step failure" : ""
        println("line $line: retcode = $(record["retcode"]), " *
                "resid0 = $(round(record["resid0"], sigdigits=4)), " *
                "i_pre = $(round(record["i_pre"], sigdigits=4))$flag")
    end

    meta = Dict("data_file" => data_file,
                "dt" => dt,
                "method" => String(method),
                "max_iter" => max_iter,
                "mult" => mult,
                "settle_steps" => settle_steps,
                "post_steps" => post_steps,
                "skipped_lines" => skipped)

    mkpath(out_dir)
    out_file = joinpath(out_dir, replace(basename(data_file), ".xlsx" => "_n1.json"))
    open(out_file, "w") do io
        JSON.print(io, Dict("meta" => meta, "cases" => records), 2)
    end

    first_step = filter(r -> r["first_step_failure"], records)
    other = filter(r -> !r["first_step_failure"] && r["retcode"] != "Success", records)
    sort!(first_step, by=r -> -r["resid0"])

    println("================================================")
    println("Cases run: $(length(records))  (skipped: $(length(skipped)))")
    println("Successful: $(count(r -> r["retcode"] == "Success", records))")
    println("First-step Newton failures: $(length(first_step))")
    for r in first_step
        println("    line $(r["line"]) ($(r["bus1"])-$(r["bus2"])): " *
                "resid0 = $(round(r["resid0"], sigdigits=4)), i_pre = $(round(r["i_pre"], sigdigits=4))")
    end
    println("Other failures: $(length(other))")
    for r in other
        where_failed = r["steps_completed"] === nothing ? "an unknown step" :
                       "$(r["steps_completed"]) step(s)"
        println("    line $(r["line"]) ($(r["bus1"])-$(r["bus2"])): " *
                "retcode = $(r["retcode"]) after $where_failed")
    end
    println("Results written to $out_file")
    println("================================================")

    return records
end

records = main()
