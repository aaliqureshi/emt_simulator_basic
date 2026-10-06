# Power flow from flat initial conditions

`solve_power_flow_continuation!` solves the existing Barq power-flow equations
from a flat start, without using the PQ voltage or non-slack angle guesses
stored in the case file.

```julia
using Barq

models = load_data("cases/Fault_Cases/case3012wp_barq.xlsx")
sys = build_system(models)
pf = solve_power_flow_continuation!(sys)
pf.converged || error("Power flow failed: $(pf.retcode), lambda=$(pf.lambda)")

# Successful continuation has already populated the bus solution.
run_static_init!(sys)
```

This call can replace `solve_power_flow!(sys)` in a simulation script. The
existing direct Newton solver remains available and its behavior is unchanged.

## Equations and initialization

Let `x = [V_PQ; theta_non_slack]` and let `F(x)` be the existing residual
`[Q_PQ; P_non_slack]`, including network injections, loads, and generation.
Define

```text
x_flat = [ones(n_PQ); zeros(n_non_slack)]
r_flat = F(x_flat)
H(x, lambda) = F(x) - (1 - lambda) r_flat = 0.
```

At `lambda = 0`, `x_flat` is an exact root. At `lambda = 1`, this is precisely
the original power flow. Intermediate stages interpolate specified net
injections from those consistent with the flat voltage profile to the target
injections. Those intermediate injections need not represent a practical
generation dispatch; they provide a numerical path to the target solution.

PV and slack voltage magnitudes, and the slack angle, remain prescribed.
Thus “flat” applies to the unknowns; setting every bus magnitude to 1 would
change the voltage controls. The current implementation uses zero non-slack
angles, matching the existing direct solver's flat initialization.

Simply scaling loads and generation to zero would not guarantee that this
flat profile solves the starting equations: line charging, shunts, transformer
taps, and unequal prescribed voltages can produce nonzero network injections.
Subtracting `r_flat` accounts for all effects already represented in `F`.

## Tracking the solution

At an accepted point, the tangent satisfies

```text
J(x) t = -r_flat,             J = dF/dx
x_predict = x + delta_lambda t.
```

Damped Newton then corrects the prediction at the new fixed `lambda`.
The implementation uses an analytic sparse Jacobian, backtracking to reduce
the residual, and rejects nonpositive PQ voltage magnitudes. A failed
correction is discarded and retried from the previous accepted point with
half the continuation step.

The default initial step is 0.05, maximum step 0.2, minimum step 1e-6, with
20 Newton iterations per corrector and at most 200 continuation attempts.
Steps grow by 1.5 after corrections taking at most 4 iterations and shrink
after corrections taking at least 8 iterations. The absolute infinity-norm
residual tolerance is 1e-9 p.u.

```julia
pf = solve_power_flow_continuation!(sys;
    initial_step=0.05, max_step=0.2, min_step=1e-6,
    abstol=1e-9, maxiters=20, maxsteps=200, verbose=true)
```

Check `pf.converged` before static initialization or simulation. On failure,
the system is unchanged and `pf.u` is the last accepted continuation state;
`pf.lambda` identifies that state's parameter. `pf.residual_norm` always
measures the original target residual `F`, while each history entry reports
the attempted stage's `H` residual. `pf.retcode` is a Symbol, rather than a
NonlinearSolve return-code object.

## Verified case3012wp result

Tested on the local `cases/Fault_Cases/case3012wp_barq.xlsx` with 5,725 unknowns:

| Method | Result | Final target mismatch, infinity norm |
| --- | --- | ---: |
| Newton, data-file guess | Success | 9.76e-12 p.u. |
| Newton, flat guess, 20 iterations | MaxIters | 5.46e-1 p.u. |
| Continuation, flat guess | Success | 1.36e-10 p.u. |

Both direct baselines use NonlinearSolve's undamped Newton with an analytic
sparse Jacobian for this comparison. The continuation run accepted all seven
steps, at lambda values 0.05, 0.125, 0.2375, 0.40625, 0.60625, 0.80625, and 1.
It used 22 total corrector iterations. Its maximum state difference from the
data-start Newton solution was 2.62e-12 (magnitudes in p.u., angles in radians).

Reproduce the comparison from the repository root:

```sh
julia --project=. scripts/power_flow_continuation.jl
julia --project=. tests/test_power_flow_continuation.jl
```

The comparison writes its history to
`outputs/power_flow_continuation/case3012wp_barq_history.csv`. The tests cover
the existing residual and direct solver, the analytic Jacobian against
automatic differentiation, independence from file guesses, step rejection,
failure without modifying the system, and successful IEEE 14/39-bus solves.

This result concerns the case as modeled by Barq; it is not a validation of
the workbook conversion against the original MATPOWER case. Generator
reactive limits and PV-to-PQ switching are not enforced by the existing
formulation or this continuation wrapper.

Natural-parameter continuation is not guaranteed to reach every target. If
the path encounters a fold and the Jacobian becomes singular, shrinking the
step may be insufficient; pseudo-arclength continuation would be a further
extension. See the [MATPOWER continuation power-flow description](https://matpower.app/manual/matpower/ContinuationPowerFlow.html)
and [parameterization discussion](https://matpower.app/manual/matpower/Parameterization.html)
for the general predictor-corrector framework.
