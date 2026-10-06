using Pkg; Pkg.activate(".")
using Revise
using Barq
using MyDiffEq
using LinearAlgebra, Printf

# Fault-on re-initialization on the 39-bus system with a WECC motor D share at
# every load bus. The remainder of each load is ZIP with the PSS/E PQBRAK
# low-voltage characteristic. For each motor D fraction the script records:
#   - direct Newton from the pre-fault point (converged?, iterations),
#   - integrator steps of decreasing size from the pre-fault point (retcode),
#   - natural-parameter continuation (converged?, total Newton iterations, fold),
#   - voltage at the faulted bus and the number of load buses in run state II.

data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
models = load_data(data_file)
sys = build_system(models)

models.fault.bus = [24]
x_fault_list = [0.015, 0.02, 0.025, 0.03, 0.04, 0.05]

solve_power_flow!(sys)
run_static_init!(sys)
address = build_dynamic_address(sys)
mass_matrix = build_mass_matrix(sys, address)
u0 = build_initial_conditions(sys, address)

p_direct = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, 1.0)
p_base = p_direct[1:end-1]
vd_fault_idx = address["balance_d"][models.fault.bus[1]]
vq_fault_idx = address["balance_q"][models.fault.bus[1]]

# reviewer-requested load model for the non-motor share
load_kwargs = (zip=(0.7, 0.1, 0.2), pqbrak=0.7, characteristic=1, low_voltage=true)
fractions = vcat(0.0:0.1:0.4, 0.45:0.01:0.7)
dt_list = [5e-4, 5e-5, 5e-6]
outdir = "outputs/motor_d_sweep_39"
mkpath(outdir)

"""
Return the DAE residual with a motor D share `f` of every load, remainder ZIP.
"""
function make_model(f)
    motor = f > 0 ? motor_d_load(f) : nothing
    return (du, u, p, t) -> begin
        solve_generator!(du, u, p)
        solve_line!(du, u, p, t)
        solve_fault!(du, u, p, t)
        balance!(du, u, p; load_kwargs..., motor=motor)
    end
end

"""
Voltage magnitudes of all buses from a state vector (slack bus from the power flow).
"""
function bus_voltages(u)
    v = copy(models.bus.v)
    v[sys.non_slack_buses] = hypot.(u[address["balance_d"]], u[address["balance_q"]])
    return v
end

"""
Integrate two steps of size `dt` from `u` with the fault applied; return retcode
(`:Diverged` if the step Newton produced a non-finite Jacobian).
"""
function step_retcode(model!, u, dt)
    prob = MyDiffEq.ODEProblem(model!, copy(u), (0.0, 2dt), p_direct, mass_matrix)
    sol = try
        MyDiffEq.Solve(prob, dt, method=:Euler, adaptive=false, tstops=[], always_new=true)
    catch e
        e isa ArgumentError || rethrow()
        return :Diverged
    end
    return sol.retcode
end

n, alg_idx = Barq.DynamicSim._setup_alg(u0, address)
prm = Barq.Models.BusModel.MOTOR_D
vstallbrk = Barq.Models.BusModel._motor_d_vstallbrk(prm)
load_buses = unique(models.load.bus)
@printf("motor D run state II band: (%.3f, %.3f)\n", vstallbrk, prm.vbrk)

"""
Faulted-bus voltage, number of load buses in run state II, number stalled, and
cond(J) at a converged point `u`; NaNs if not converged.
"""
function root_summary(u, converged, model!)
    converged || return (NaN, -1, -1, NaN)
    v = bus_voltages(u)
    n_run2 = count(b -> vstallbrk < v[b] < prm.vbrk, load_buses)
    n_stall = count(b -> v[b] <= vstallbrk, load_buses)
    _, J = Barq.DynamicSim._eval_g_jac(u, p_direct, n, alg_idx, model!)
    kappa = all(isfinite, J) ? cond(J) : NaN
    return (hypot(u[vd_fault_idx], u[vq_fault_idx]), n_run2, n_stall, kappa)
end

rows = Any[]
for xf in x_fault_list, f in fractions
    models.fault.x_fault[1] = xf
    model! = make_model(f)

    # pre-fault consistency of the load split
    du = zeros(n)
    model!(du, u0, (p_base..., 0.0), 0.0)
    res0 = norm(du[alg_idx], Inf)

    # direct Newton from the inherited point
    u_n = copy(u0)
    r_n = solve_newton!(u_n, p_direct, address; max_iter=50, always_new=true, model! = model!)
    eta = r_n.correction_norm[1]

    # step-size reduction from the inherited point
    codes = [step_retcode(model!, u0, dt) for dt in dt_list]

    # continuation from the inherited point
    u_h = copy(u0)
    r_h = redirect_stdout(devnull) do
        solve_adaptive_homotopy!(u_h, p_base, address; always_new=true, Δλ_init=0.01, Δλ_max=0.02,
                                 vd_idx=vd_fault_idx, vq_idx=vq_fault_idx, model! = model!)
    end
    open(joinpath(outdir, @sprintf("path_xf%.3f_f%.2f.csv", xf, f)), "w") do io
        println(io, "lambda,v_fault")
        for (lam, vd, vq) in zip(r_h.λ_hist, r_h.vd_hist, r_h.vq_hist)
            println(io, lam, ',', hypot(vd, vq))
        end
    end

    vn, run2_n, stall_n, cond_n = root_summary(u_n, r_n.converged, model!)
    vh, run2_h, stall_h, cond_h = root_summary(u_h, r_h.converged, model!)
    root_gap = (r_n.converged && r_h.converged) ? norm(u_n - u_h, Inf) : NaN
    lam_fail = r_h.λ_failed === nothing ? NaN : r_h.λ_failed

    @printf("xf=%.3f f=%.2f res0=%.0e | NR %s %2d it eta=%.1f V24=%.3f run2=%2d stall=%2d cond=%.0e | steps %s | HC %s %3d it lam_fail=%.3f V24=%.3f run2=%2d stall=%2d cond=%.0e | gap=%.1e\n",
            xf, f, res0, r_n.converged ? "ok  " : "FAIL", r_n.iters, eta, vn, run2_n, stall_n, cond_n,
            join(string.(codes), ","), r_h.converged ? "ok  " : "FAIL",
            r_h.total_newton_iters, lam_fail, vh, run2_h, stall_h, cond_h, root_gap)
    push!(rows, (xf, f, res0, r_n.converged, r_n.iters, eta, vn, run2_n, stall_n, cond_n, codes...,
                 r_h.converged, r_h.total_newton_iters, lam_fail, vh, run2_h, stall_h, cond_h, root_gap))
end

header = ["x_fault", "motor_d", "prefault_residual", "newton_converged", "newton_iters", "first_step_norm",
          "newton_v_fault", "newton_n_run2", "newton_n_stall", "newton_cond_J",
          ["retcode_dt_$(dt)" for dt in dt_list]..., "homotopy_converged", "homotopy_iters",
          "lambda_failed", "homotopy_v_fault", "homotopy_n_run2", "homotopy_n_stall", "homotopy_cond_J",
          "root_gap"]
open(joinpath(outdir, "summary.csv"), "w") do io
    println(io, join(header, ','))
    for r in rows
        println(io, join(r, ','))
    end
end
println("written ", joinpath(outdir, "summary.csv"))
