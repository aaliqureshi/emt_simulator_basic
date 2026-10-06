# Run from the repository root:
# julia --project=. scripts/power_flow_continuation.jl
using Barq, LinearAlgebra, NonlinearSolve
using Barq.PowerFlow: PowerFlowCache, powerflow!, _powerflow_jacobian

function compare_power_flow_starts(filename)
    sys = build_system(load_data(filename))
    flat = vcat(ones(length(sys.v_update_buses)), zeros(length(sys.non_slack_buses)))
    data = vcat(sys.models.bus.v[sys.v_update_buses], sys.models.bus.theta[sys.non_slack_buses])
    cache = PowerFlowCache(sys, flat)
    # Same equations and undamped Newton method as solve_power_flow!, using
    # an analytic sparse Jacobian to keep this large-case comparison inexpensive.
    jac! = (J, u, p) -> copyto!(J, _powerflow_jacobian(u, p))
    f = NonlinearFunction{true, NonlinearSolve.SciMLBase.FullSpecialize}(
        powerflow!; jac=jac!, jac_prototype=_powerflow_jacobian(flat, cache))
    direct(u) = solve(NonlinearProblem(f, u, cache), NewtonRaphson();
                      abstol=1e-9, reltol=1e-9, maxiters=20)
    from_data = direct(data)
    from_flat = direct(flat)
    println("Case: $filename; unknowns: $(length(flat))")
    println("Newton/data: $(from_data.retcode), residual=$(norm(from_data.resid, Inf))")
    println("Newton/flat: $(from_flat.retcode), residual=$(norm(from_flat.resid, Inf))")

    result = solve_power_flow_continuation!(sys)
    println("Continuation/flat: $(result.retcode), lambda=$(result.lambda), " *
            "residual=$(result.residual_norm), corrector iterations=$(result.iterations)")
    println("Accepted steps: $(count(h -> h.accepted, result.history)); " *
            "rejected steps: $(count(h -> !h.accepted, result.history))")
    result.converged || error("Continuation did not reach the target power flow")
    NonlinearSolve.SciMLBase.successful_retcode(from_data) || error("Data-start reference failed")
    difference = norm(result.u - from_data.u, Inf)
    println("Maximum state difference from data-start solution: $difference")
    @assert difference < 1e-7
    residual = similar(flat)
    powerflow!(residual, result.u, PowerFlowCache(sys, result.u))
    @assert norm(residual, Inf) <= 1e-9

    output_dir = joinpath(@__DIR__, "..", "outputs", "power_flow_continuation")
    mkpath(output_dir)
    output = joinpath(output_dir, splitext(basename(filename))[1] * "_history.csv")
    open(output, "w") do io
        println(io, "lambda,step,accepted,iterations,homotopy_residual_inf")
        for h in result.history
            println(io, "$(h.lambda),$(h.step),$(h.accepted),$(h.iterations),$(h.residual_norm)")
        end
    end
    println("History written to: $output")
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    filename = isempty(ARGS) ? joinpath(@__DIR__, "..", "cases", "Fault_Cases", "case3012wp_barq.xlsx") : ARGS[1]
    compare_power_flow_starts(filename)
end
