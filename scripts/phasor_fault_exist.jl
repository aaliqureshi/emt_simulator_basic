# Does the fault-on steady state exist? All-algebraic (phasor) network: only delta, omega differential.
using Barq, LinearAlgebra, Printf
models = load_data("cases/Fault_Cases/ieee39_fault.xlsx"); sys = build_system(models)
fb = parse(Int, ARGS[1]); models.fault.bus = [fb]; models.fault.x_fault[1] = parse(Float64, ARGS[2])
solve_power_flow!(sys); run_static_init!(sys)
address = build_dynamic_address(sys); u0 = build_initial_conditions(sys, address)
nsb = sys.non_slack_buses; k = findfirst(==(fb), nsb); bd, bq = address["balance_d"], address["balance_q"]
pb = (address, sys.models, sys.incidence_matrix, sys.C_eq, nsb)
alg_addr = Dict("delta"=>address["delta"], "omega"=>address["omega"], "line_id"=>1:0, "line_iq"=>1:0, "balance_d"=>1:0, "balance_q"=>1:0)
mk(kw) = (du, u, p, t) -> (solve_generator!(du, u, p); solve_line!(du, u, p, t); solve_fault!(du, u, p, t); balance!(du, u, p; kw...))
for (label, kw) in (("default kwargs", (;)), ("ZIP no low-voltage", (; zip=(0.7,0.1,0.2), low_voltage=false)),
                    ("ZIP + PQBRAK 0.7", (; zip=(0.7,0.1,0.2), low_voltage=true, pqbrak=0.7, characteristic=1)), ("pure Z", (; zip=(1.0,0.0,0.0), low_voltage=false)))
    f = mk(kw)
    un = copy(u0); rn = solve_newton!(un, (pb..., 1.0), alg_addr; max_iter=50, always_new=true, model! = f)
    uh = copy(u0); rh = solve_homotopy!(uh, pb, alg_addr; Δλ=0.01, always_new=true, model! = f)
    @printf "bus %d xf=%s %-20s | Newton conv=%s it=%d V%d=%.3f | natural homotopy conv=%s λ_failed=%s V%d at last λ=%.3f\n" fb ARGS[2] label rn.converged rn.iters fb hypot(un[bd[k]], un[bq[k]]) rh.converged string(rh.λ_failed) fb hypot(uh[bd[k]], uh[bq[k]])
end
