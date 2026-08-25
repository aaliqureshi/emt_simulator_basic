using Pkg

# Examples:
#   julia --project=. scripts/compare_power_flow_solvers.jl
#
# Reproduce an NR-fails/PGD-converges case at a strict residual tolerance:
#   julia --project=. scripts/compare_power_flow_solvers.jl \
#     --scales=1.0 --initial-voltages=0.6 --initial-angles=0.0 \
#     --tol=1e-8 --nr-max=100 --pgd-max=100000 --pgd-alpha=0.2

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT)

include(joinpath(PROJECT_ROOT, "src", "data_loader.jl"))
include(joinpath(PROJECT_ROOT, "src", "system.jl"))
include(joinpath(PROJECT_ROOT, "src", "power_flow.jl"))

using .DataLoader
using .SystemModel
using .PowerFlow
using CSV
using DataFrames
using LinearAlgebra
using MyDiffEq
using Printf

const DEFAULT_CASE = joinpath(PROJECT_ROOT, "cases", "Fault_Cases",
                              "ieee14_fault_barq.xlsx")

"""Return a MyDiffEq power-flow problem without changing the system state."""
function power_flow_problem(base_models, load_scale;
                            scale_reactive=false,
                            initial_voltage=1.0,
                            initial_angle=0.0)
    models = deepcopy(base_models)
    models.load.p .*= load_scale
    scale_reactive && (models.load.q .*= load_scale)
    sys = build_system(models)

    nq = length(sys.v_update_buses)
    np = length(sys.non_slack_buses)
    address = Dict(
        "balance_reactive_power" => 1:nq,
        "balance_real_power" => (nq + 1):(nq + np),
    )
    u0 = vcat(fill(initial_voltage, nq), fill(initial_angle, np))
    p = (address, models, sys.Y, sys.non_slack_buses, sys.v_update_buses)
    prob = MyDiffEq.NonLinearProblem(PowerFlow.powerflow, u0, p)
    return prob, address, sys
end

"""
Run one solver and independently validate the final residual.

`tol` is made an absolute 2-norm target for both solvers. MyDiffEq currently
uses `rtol * max(1, norm(F(u0)))`, so `rtol` is adjusted per problem.
"""
function run_solver(prob, address, method;
                    tol=1e-6, nr_max_iter=100, pgd_max_iter=150_000,
                    pgd_alpha=1.0)
    initial_residual = PowerFlow.powerflow(prob.u0, prob.p)
    initial_norm_2 = norm(initial_residual, 2)
    adjusted_rtol = tol / max(1.0, initial_norm_2)
    max_iter = method == :NewtonRaphson ? nr_max_iter : pgd_max_iter
    alpha = method == :NewtonRaphson ? 1.0 : pgd_alpha

    sol = nothing
    elapsed = @elapsed redirect_stdout(devnull) do
        sol = MyDiffEq.NLsolve(
            prob;
            method=method,
            max_iter=max_iter,
            atol=tol,
            rtol=adjusted_rtol,
            alpha=alpha,
            verbose=false,
        )
    end

    final_residual = PowerFlow.powerflow(sol.u_final, prob.p)
    final_norm_2 = norm(final_residual, 2)
    final_norm_inf = norm(final_residual, Inf)
    finite = all(isfinite, sol.u_final) && all(isfinite, final_residual)
    validated = finite && final_norm_2 <= tol

    voltage = sol.u_final[address["balance_reactive_power"]]
    return (
        solver=String(method),
        reported_converged=sol.converged,
        validated_converged=validated,
        retcode=String(sol.retcode),
        iterations=sol.iters,
        initial_norm_2=initial_norm_2,
        final_norm_2=final_norm_2,
        final_norm_inf=final_norm_inf,
        min_voltage=isempty(voltage) ? NaN : minimum(voltage),
        max_voltage=isempty(voltage) ? NaN : maximum(voltage),
        elapsed_seconds=elapsed,
    )
end

function failure_result(method, err)
    return (
        solver=String(method),
        reported_converged=false,
        validated_converged=false,
        retcode="exception: $(typeof(err))",
        iterations=0,
        initial_norm_2=NaN,
        final_norm_2=NaN,
        final_norm_inf=NaN,
        min_voltage=NaN,
        max_voltage=NaN,
        elapsed_seconds=NaN,
    )
end

function parse_options(args)
    options = Dict{String,String}()
    for arg in args
        startswith(arg, "--") || error("Expected --name=value, got: $arg")
        pieces = split(arg[3:end], "="; limit=2)
        length(pieces) == 2 || error("Expected --name=value, got: $arg")
        options[pieces[1]] = pieces[2]
    end
    return options
end

parse_bool(value) = lowercase(value) in ("true", "yes", "1")

function main(args=ARGS)
    options = parse_options(args)
    case_file = get(options, "case", DEFAULT_CASE)
    scales = parse.(Float64, split(get(options, "scales", "1.0,4.4,4.465"), ','))
    initial_voltages = parse.(Float64, split(get(options, "initial-voltages", "1.0"), ','))
    initial_angles = parse.(Float64, split(get(options, "initial-angles", "0.0"), ','))
    scale_reactive = parse_bool(get(options, "scale-q", "false"))
    tol = parse(Float64, get(options, "tol", "1e-6"))
    nr_max_iter = parse(Int, get(options, "nr-max", "100"))
    pgd_max_iter = parse(Int, get(options, "pgd-max", "150000"))
    pgd_alpha = parse(Float64, get(options, "pgd-alpha", "1.0"))
    output_file = get(options, "output",
                      joinpath(PROJECT_ROOT, "power_flow_solver_comparison.csv"))

    println("Case: $case_file")
    println("Load scaling: P$(scale_reactive ? " and Q" : " only")")
    @printf("Validation target: ||F||₂ ≤ %.3e\n\n", tol)

    base_models = load_data(case_file)
    rows = NamedTuple[]
    for scale in scales, initial_voltage in initial_voltages, initial_angle in initial_angles
        prob, address, sys = power_flow_problem(
            base_models, scale;
            scale_reactive=scale_reactive,
            initial_voltage=initial_voltage,
            initial_angle=initial_angle,
        )
        for method in (:NewtonRaphson, :GradientDescent)
            result = try
                run_solver(
                    prob, address, method;
                    tol=tol,
                    nr_max_iter=nr_max_iter,
                    pgd_max_iter=pgd_max_iter,
                    pgd_alpha=pgd_alpha,
                )
            catch err
                failure_result(method, err)
            end
            row = merge((load_scale=scale, initial_voltage=initial_voltage,
                         initial_angle=initial_angle), result)
            push!(rows, row)
            @printf(
                "scale=%7.4f  V0=%6.3f  θ0=%7.3f  %-15s  converged=%-5s  iters=%7d  ||F||₂=%10.3e  Vmin=%8.4f  time=%8.3fs  %s\n",
                scale, initial_voltage, initial_angle, result.solver,
                string(result.validated_converged), result.iterations,
                result.final_norm_2, result.min_voltage,
                result.elapsed_seconds, result.retcode,
            )
        end
    end

    results = DataFrame(rows)
    CSV.write(output_file, results)
    println("\nWrote $(abspath(output_file))")

    separations = DataFrame()
    for scale in scales, initial_voltage in initial_voltages, initial_angle in initial_angles
        pair = results[(results.load_scale .== scale) .&
                       (results.initial_voltage .== initial_voltage) .&
                       (results.initial_angle .== initial_angle), :]
        nr_ok = only(pair[pair.solver .== "NewtonRaphson", :validated_converged])
        pgd_ok = only(pair[pair.solver .== "GradientDescent", :validated_converged])
        if !nr_ok && pgd_ok
            append!(separations, pair)
        end
    end
    if nrow(separations) == 0
        println("No NR-fails/PGD-converges case was found in this sweep.")
    else
        println("NR-fails/PGD-converges configurations:")
        display(unique(separations[:, [:load_scale, :initial_voltage, :initial_angle]]))
    end
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
