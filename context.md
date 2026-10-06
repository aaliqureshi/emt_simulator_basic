# Context: re-initialization after discontinuous events (Newton vs homotopy)

Working notes, last updated 2026-09-29 01:20. Purpose: let a new chat pick up the
investigation into Newton failure at DAE re-initialization after a fault, the
homotopy alternatives, and the strategy for revising the paper. Section "Current
state" is authoritative where it disagrees with the older findings below.

## Current state of the working tree (2026-09-29)

Edits the user made after the first draft of these notes:

- `src/dynamic_sim.jl`, `build_mass_matrix`: the line inductance entries for
  `line_id` / `line_iq` are commented out. `MyDiffEq.ODEProblem` classifies a row
  as algebraic when its mass-matrix diagonal is zero, so the line equations in
  `solve_line!` (`V1/tap - V2 - (R + jX) I = 0`) are now enforced as quasi-static
  phasor branch equations by the integrator. Combined with classical generators
  behind `x_d'` and the ZIP bus balance, the in-repo model is now the algebraic
  network that the `/tmp` phasor prototype (see below) approximated.
- `src/algebraic_solvers.jl`, `_setup_alg`: NOT yet updated. It still excludes
  `line_id`, `line_iq` from the algebraic set (the alternative line with
  `n_mech = n_delta + n_omega` is present but commented out). With the new mass
  matrix this is inconsistent: re-init freezes the line currents and leaves the
  line residuals unsatisfied, and the first integrator step then has to solve
  them. First thing to do: switch `_setup_alg` to the `delta + omega` only form.
  `alg_idx = n_mech+1:n` relies on the address order `delta, omega, line_id,
  line_iq, balance_d, balance_q, fault_id, fault_iq, gen_id, gen_iq`, which holds.
- `src/models/bus.jl`, `balance!` defaults: `zip=(0.7, 0.1, 0.2)` (sums to 1, the
  earlier 0.9 issue is fixed), `pqbrak=0.7e-6`, `characteristic=2`,
  `low_voltage=true`. Note that for characteristic 2 the constant-current factor
  `kI` is independent of `pqbrak` and kicks in below `v = 0.5`, so the small
  `pqbrak` only disables the constant-power reduction. PSS/E default `pqbrak=0.7`
  has not been adopted.
- `src/main_convergence_failure_testing.jl`: `ieee39_fault.xlsx`, fault at bus 24,
  `x_fault = 0.101` (was 0.015). Computes `r1` (direct Newton), `r2`
  (`solve_homotopy!`, `dlambda=0.01`), `r3` (`solve_adaptive_homotopy!`), then
  the post-re-init simulation and the plotting blocks. The `r4` residual-homotopy
  block was removed from the driver; `build_residual_homotopy` remains in
  `algebraic_solvers.jl` and is exported.

Everything below about "frozen line currents" and the localization bound
describes the EMT (inductive line) configuration and no longer applies once the
lines are algebraic. Expect the 39-bus fault case to behave like the phasor
prototype: the fault pulls V24 down far less, Newton converges for most load
models, and the interesting regime has to be created by pushing the post-event
algebraic system toward its existence limit.

## Repository pointers

- `src/algebraic_solvers.jl`: `solve_newton!`, `solve_homotopy!` (fixed-step
  natural parameter), `solve_adaptive_homotopy!` (secant predictor + Newton
  corrector + PI step control), `_newton_corrector!`, `build_residual_homotopy`
  (section 8). All solvers take a `model!(du, u, p, t)` keyword (default
  `solve_dynamic_sim!`) and read the homotopy parameter from `p[end]`.
  `src/algebraic_solvers_OG.jl` is an untracked copy of the pre-session version.
- `src/models/fault.jl`: the natural homotopy parameter enters only here,
  `x_eff = 1e10^(1-lambda) * x_fault^lambda`, a geometric sweep of the fault
  susceptance over ~12 decades.
- `src/models/bus.jl`: `balance!(du, u, p; zip, pqbrak, characteristic,
  low_voltage)` writes the bus equations in power-mismatch form with a ZIP load
  and the OpenIPSL/PSS/E low-voltage characteristic (`_openipsl_load_factors`).
  The previous implementation is kept in a block comment above it.
- `src/models/line.jl`: `solve_line!` residual is the R-L branch equation; with
  zero mass it is the phasor branch equation.
- `src/dynamic_sim.jl`: `solve_dynamic_sim!` = generator, line, fault, balance.
- `src/homotopy.jl`: `pseudo_arclength!(u1, lam1, u2, lam2, p_base, mass_matrix,
  address; ...)` and `solve_algebraic!`, exported from `Barq.jl` (used in
  `src/main.jl`, not in the driver).
- `src/power_flow_continuation.jl`: `solve_power_flow_continuation!` (loading
  factor continuation for power flow; separate from the re-init work).
- Heavier cases for a loading sweep: `cases/Fault_Cases/case118_gc.xlsx`,
  `case2383wp_gc.xlsx`, `case3012wp_barq.xlsx`, `case13659pegase_barq.xlsx`.
- `MyDiffEq` lives at `/Users/aali27/Work/Research/DAE_Reinitialization/
  DAEReinitialization/MyDiffEq/MyDiffEq` (dev dependency). Implicit Euler and
  trapezoid steppers solve `M (u - u_prev) - dt f(u) = 0` with Newton; `sol.
  newton_log.residual_norm` holds the per-iteration residuals used in the plots.

## Code changes made by the assistant (session of 2026-09-28)

1. `src/algebraic_solvers.jl`: added `build_residual_homotopy(u0, p_direct,
   address; model!)`, a closure implementing the Newton (residual) homotopy
   `H(u, lam) = g(u) - (1 - lam) g(u0)` on the algebraic rows, with `g` evaluated
   at the fixed target parameters `p_direct`. Its Jacobian equals that of `g`.
2. `src/dynamic_sim.jl`, `src/Barq.jl`: export `build_residual_homotopy`.
3. The `r4` block added to the driver has since been removed by the user.

All test scripts were written under `/tmp` and deleted; recipes are at the end.

## Review of the existing homotopy code

`solve_homotopy!` is a correct natural-parameter continuation with these defects:
no rollback of `u` when a stage fails (returns the iterate after `max_iter`), no
`isfinite` guards on `g` and the Newton step, `0.0:dlambda:lambda_target` only
lands on `lambda_target` when `dlambda` divides it, and the post-update residual
of a converged stage is not recorded (plotted curves are offset by one iteration
relative to `solve_newton!`).

`solve_adaptive_homotopy!`: the assignments `dlambda_max=0.2` and `n_target=2`
right after the signature silently override the keyword arguments; `n_target=2`
makes the PI controller shrink the step whenever Newton needs 3 iterations (this
is why failing runs spent ~600-700 iterations bisecting toward a fold).
`total_jac_evals += r.iters` is wrong when `always_new=false`. Still not fixed.

Parameterization: with `x_fault=0.015` the geometric sweep spans 11.8 decades;
the network only responds once `b = 1/x_eff` is of order 1 pu, i.e. for
`lambda > ~0.8`. Roughly 80% of uniform stages re-solve the pre-fault point.

## Residual (Newton) homotopy: interpretation and behavior

With `g0 = g(u_prefault)` under fault-on equations, `g0` is dominated by the two
fault rows `[b vd0, b vq0]`. Substituting into `H = 0`, the fault branch becomes
a reactance `x_f` connected to a voltage source `(1 - lam) V0` shrinking to
zero. With stiff `b` the fault-bus voltage tracks `(1 - lam) V0` linearly, so
the path has no dead zone and `dH/dlam = -g0` is constant. It works well for
stiff faults; for a weak fault (`x_f = 0.15`) the residual path folded while the
natural homotopy converged, so it is not uniformly more robust.

## Findings for the EMT configuration (inductive lines; superseded by current tree)

39-bus, fault at bus 24, `x_fault = 0.015`, static-init start. Because line
currents were frozen during re-init, the problem reduced to two scalar equations
at the faulted bus. With `|I_lines| = 3.16 pu` and `b = 66.7 pu`, the reactive
balance forces `|V24| <~ 0.047`, so deliverable active power at bus 24 is at
most ~0.15 pu (empirical existence limit: constant-P share `k_p ~ 0.025`).
Constant-P shares of 0.1 or 0.2 have no solution. The Newton divergence in the
original script was non-existence, not a solver weakness.

| load model | solution exists | NR (power form) | NR (current form) | homotopies |
|---|---|---|---|---|
| pure Z (1,0,0) | yes | 2 it. | - | converge, match NR to 1e-15 |
| ZIP (.7,.1,.2), pqbrak off | no | diverges | diverges | fold at lam*=0.877 (resid), 0.918 (natural) |
| k_p >= 0.026, pqbrak off | no | diverges | diverges | fold |
| k_p <= 0.024 | yes | 8-10 it. | 6-10 it. | same root |
| ZIP (.7,.1,.2), pqbrak=0.7 (char 1 or 2) | yes | converges to V24=0 (char1) / 0.002 (char2), fault current ~0 | 5 it., physical root | physical root V24=0.049, I_f=3.28 |
| ZI (.5,.5,0) or (.3,.7,0), pqbrak off | yes | converges to V24=0 | 6 it., physical root | physical root |
| ZIP (.7,.1,.2), x_f=0.15 | yes | converges | - | natural OK, residual folds |

The `V=0` root is an artifact of the power-mismatch form: at `V=0` both `P` and
`Q` mismatches vanish regardless of current mismatch. It is not a root of the
current-mismatch equations. This artifact persists in the current `balance!`.

Conclusion for the EMT model: the inductive network localizes the inconsistency
to the faulted bus; the local problem is nearly quadratic and Newton is robust
whenever a solution exists. No "region B" (exists but Newton fails) was found.

## Phasor (algebraic-network) prototype: this is what the repo now implements

A `/tmp` script built the same event with an algebraic Y-bus network, classical
generators with frozen `E' = e_q_prime*exp(j delta)` behind `x_d_prime`, slack
bus voltage fixed, ZIP loads with `_openipsl_load_factors`, fault as shunt
`-j b` at bus 24, current-injection residual in `[vd; vq]` of the non-slack
buses, Newton via ForwardDiff, plus a residual homotopy with step halving.

Findings: at `x_f = 0.015` the phasor network only pulls V24 to 0.48 and Newton
converges in 4-6 iterations for every load model at both fault application and
clearing. Pure constant-P without low-voltage characteristic has no fault-on
solution for any `x_f` (fold at lam ~ 0.72-0.78).

Sweeping `pqbrak` toward the existence limit with constant-P loads
(`x_f = 0.005`): Newton iterations 6, 6, 9, 25, 19, 9, 11, 14, then fail at
`pqbrak = 0.32` while the solution exists (found by homotopy), then
non-existence at 0.30 and below. At `x_f = 0.002`: 6, 6, 8, 9, 11, 17, 87, 51,
90, fail. The natural-parameter residual homotopy folded at `pqbrak = 0.36,
0.34, 0.30` where Newton succeeded (S-shaped branch; pseudo-arclength needed).

## Strategy recommended for the paper revision

The paper's claim: after a mode transition the inherited point is displaced from
the post-event constraint manifold by an amount step-size reduction cannot
shrink, and whether Newton recovers depends on the nonlinearity of the
post-event algebraic map, not on the displacement itself.

1. Concede that the original 39-bus EMT example is a non-existence case; present
   the continuation's fold detection (`lam*`, `V24(lam*)`) as a certificate that
   the fault-on system has no consistent initialization for that load model.
2. Report the localization result: in the inductive network only the faulted bus
   is displaced, the displacement is huge (1.0 -> 0.05 pu) yet Newton recovers
   because the local Newton map is benign. This supports "rather than on the
   inconsistency itself".
3. Build the demonstration where the theory says it must live: a post-event
   algebraic system parametrized toward its existence limit, so three regimes
   appear: A (exists, Newton converges), B (exists, Newton fails at any step
   size, continuation converges), C (no solution, continuation certifies the
   fold). Region B exists in the phasor formulation, which the repo now runs.
   Physically defensible knob: system loading factor with `pqbrak = 0.7` on a
   heavier case (case118, case2383wp); `pqbrak` itself is the cleanest numerical
   knob but harder to defend.
4. Figure: x-axis = the parameter; displacement norm (flat), Newton iterations
   from the inherited point (erratic, hits the cap), continuation cost (bounded
   until the fold). Second panel: Newton outcome identical for h spanning three
   decades.
5. Use `pseudo_arclength!` for the continuation so it never loses to Newton on
   S-shaped branches. Use current-form bus equations for the Newton baseline, or
   state explicitly that the power form admits a degenerate V=0 root.

## Motor D (WECC single-phase A/C) study, 2026-09-29

Motivation: the reviewer requires constant-P loads to carry PQBRAK, which makes
the ZIP map Z-like below 0.7 pu and removes the load curvature. The WECC motor D
performance model is algebraic, standard, and strongly curved in 0.52-0.86 pu:
run I (V > 0.86) P = P0, Q = Q0' + 6 (V-0.86)^2; run II (Vstallbrk < V < 0.86)
P = P0 + 12 (0.86-V)^3.2, Q = Q0' + 11 (0.86-V)^2.5; stall (V < Vstallbrk)
impedance Rstall = Xstall = 0.1. Vstallbrk = 0.523 with the defaults. PNNL's
reference CLM implementation states the network solution cannot converge with
run II enforced and replaces it by a lagged impedance.

Code: `src/models/bus.jl` gained `MOTOR_D` (defaults from the WECC spec),
`_motor_d_vstallbrk`, `_motor_d_pq`, `_motor_d_load`, and the exported
`motor_d_load(fraction; params=MOTOR_D)`, passed to `balance!` through the
`motor` keyword (default `nothing` = pure ZIP, no motor load). The motor's Q at the initial voltage is
removed from the ZIP share so the power flow point is preserved (checked:
pre-fault residual 1e-10). `src/algebraic_solvers.jl`: the hard-coded
`dlambda_max=0.2; n_target=2` overrides are now the signature defaults (same
behaviour for existing callers, but overridable), and `solve_newton!` /
`_newton_corrector!` stop on a non-finite Jacobian instead of throwing.
Script: `scripts/motor_d_sweep_39.jl` sweeps `x_fault` in {0.015, 0.02, 0.025,
0.03, 0.04, 0.05} and motor fraction 0 to 0.7 on the 39-bus fault at bus 24, with
ZIP (0.7,0.1,0.2) + PQBRAK 0.7 (characteristic 1) for the remainder; records
Newton from the pre-fault point, integrator retcodes for h = 5e-4, 5e-5, 5e-6,
natural-parameter continuation (dlambda_max = 0.02), roots, cond(J), and the
V24(lambda) path. Output: `outputs/motor_d_sweep_39/`. Run outside the sandbox
(`Pkg.activate` writes to `~/.julia/logs`), ~2 min.

Results:
- x_f >= 0.025 (V24 >= 0.59): nothing happens up to f = 0.7. Newton 4-6
  iterations, continuation agrees, cond(J) 4e3-7e3. The run II curvature is weak
  above 0.6 pu (at V = 0.7 the P increase is 3%, Q increase 11% of motor P).
- x_f = 0.02 (V24 0.55 -> 0.455): Newton 7-8 iterations up to f = 0.68 with
  cond(J) rising 8e3 -> 2e4, then at f = 0.69 the running branch disappears
  (continuation folds at lambda = 0.9999), Newton fails, all three step sizes
  fail. No window between "Newton slows" and "no solution".
- x_f = 0.015 (V24 0.48 -> 0.40, bus 24 itself below Vstallbrk): running branch
  exists to f = 0.53 (Newton 8 it.). For f >= 0.54 the natural continuation folds
  at lambda = 0.997-0.9997 and the problem becomes multi-root: Newton from the
  pre-fault point fails at f = 0.54, 0.62, 0.65, 0.69 and otherwise converges in
  10-34 iterations to different partially stalled roots (7, 8, 22, 23, 27, 28
  stalled load buses, V24 = 0.12-0.31); for f = 0.56, 0.57 Newton and the
  continuation reach different roots. Integrator outcome vs step size is
  non-monotone (f = 0.54: h = 5e-4 ok, 5e-5 MaxIter, 5e-6 ok; f = 0.60: ok,
  Diverged, MaxIter; f = 0.63: Diverged, ok, ok).

Interpretation: with WECC defaults the model behaves like the earlier load
studies, Newton is benign until the running branch folds, and the failures live
beyond the fold in a multi-root regime created by the stall switch (Q is
discontinuous at Vstallbrk). The region-B evidence available is: (i) erratic,
step-size-independent Newton outcomes beyond the fold while stalled roots exist,
and (ii) the need for pseudo-arclength to follow the running branch around the
fold to the stalled branch; the natural-parameter continuation cannot do it.
Whether a stalled fault-on solution at t = 0+ is acceptable physically is a
question for the paper (per spec the impedance branch applies instantly below
Vstallbrk; the Tstall timer only latches it).

## Open items, in order

1. Done: `_setup_alg` now includes line currents (user edit); confirmed V24 = 0.48
   at `x_f = 0.015`, Newton 4 it., matching the phasor prototype.
2. Decide how to use the motor D regime beyond the fold: run `pseudo_arclength!`
   (in `src/homotopy.jl`, takes two points on the path and `mass_matrix`) from
   the last two natural-continuation points to follow the running branch around
   the fold at `x_f = 0.015, f = 0.54-0.60`, and compare with Newton's erratic
   outcomes. Alternatively explore `MOTOR_D.lf < 1` (larger curvature on the
   motor base) or a lower `vstall` to widen the smooth run II band.
3. Fault-clearing re-init (Energies 2024 case): clearing time as the parameter,
   not yet built. `balance!` default `pqbrak=0.7e-6` is still the repo default;
   the sweep script passes `pqbrak=0.7` explicitly.
4. `solve_homotopy!` still lacks rollback and finite checks.
5. Post-clearing re-init has not been tested in the repo model.

## How to reproduce quick checks (outside the repo, e.g. /tmp)

Run with `julia --project=.` (no `Pkg.activate` needed in the sandbox):
`load_data`, `build_system`, set `models.fault.bus=[24]`, `x_fault[1]=...`,
`solve_power_flow!`, `run_static_init!`, `build_dynamic_address`,
`build_initial_conditions` (this `u0` is a valid pre-fault point; the driver's
two-step `sol_pf` is not required). `p_direct = (address, sys.models,
sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, 1.0)`; `p_base` is the same
tuple without the last element. To vary the load model without editing the repo,
pass a closure as `model!` that calls `solve_generator!`, `solve_line!`,
`solve_fault!`, then `balance!(du, u, p; zip=..., pqbrak=..., characteristic=...,
low_voltage=...)`. `Barq.Models.BusModel._openipsl_load_factors` is accessible
for custom load formulations. Newton iteration counts from the integrator are in
`sol.newton_log.residual_norm`; retcode `:MaxIter` marks a failed step.

## Fault clearing with frozen line currents (EMT lines), 2026-09-29

Quick check (scratch script, not in repo): 39-bus, bus 24, `x_f = 0.015`, pure Z
loads. Fault-on re-init Newton 2 it. (V24 = 0.048, |If| = 3.2); after 0.1 s of
fault-on integration (Trap, h = 5e-4) V24 = 0.497, |If| = 33.2, net line current
into bus 24 = 32.7 pu. Clearing re-init (lambda = 0, line currents frozen):
Z loads converge in 2 it. to V24 = 10.8 pu; ZI (0.7, 0.3, 0) 6 it., V24 = 14.8;
ZIP (0.7, 0.1, 0.2) without PQBRAK no convergence in 50 it. The ~33 pu of
trapped line current can only go into bus 24's load and quasi-static shunt, so
the consistent point is an inductive kick. It is a result of the ideal switch
combined with inductive lines and algebraic bus voltages, not a physical TRV.
Implication: with EMT lines the clearing re-init exists without the P part but
is physically implausible. With the dynamic (R-L) fault idea, application needs
no re-init and clearing still gives this kick.

## Inductive (R-L) fault cleared through a snubber resistor, 2026-09-29

Scratch script (not in repo): fault as R-L branch with mass `L_f = x_f/omega`,
applied at t = 0 from the pre-fault point with no re-init; implicit Euler,
h = 5e-5 fault-on, 1e-5 after clearing. At clearing the branch resistance
switches to `R_s` (current continuous, so no algebraic jump and no re-init).
Pure Z loads: fault-on V24 dips to 0.089 (DC offset) and settles at 0.489,
|If| = 32.6 at t_c = 0.1 s. Peak V24 after clearing vs R_s: 1e4 -> 8.8 pu,
10 -> 8.6, 3 -> 8.1, 1 -> 7.0, 0.3 -> 4.7, 0.1 -> 2.4; V24 at +30 ms about
1.03 in all cases. A snubber small enough to limit the voltage stays a large
permanent shunt (steady current about V/R_s, i.e. 1-10 pu at R_s = 1-0.1),
so the fault is not really cleared and a second interruption is still needed.
ZIP (0.7, 0.1, 0.2) + PQBRAK 0.7 char 1: fault-on OK (V24(tc) = 0.477), but the
first post-clearing step at R_s = 1e4 threw NaN in the Jacobian inside MyDiffEq
(not diagnosed; likely Newton heading toward the V = 0 root of the P part).
Next idea: time-varying (arc-like) resistance R(t) ramped to open over T_r;
the peak should scale like L_line * dI / T_r.

## Reproduction of Yao et al., TPWRS 2020 (HE dynamic simulation), Fig. 10/11, 2026-09-29

PowerSAS.m (the authors' open-source toolbox, github.com/ANL-CEEESA/powersas.m) is
cloned at `~/Work/repos/powersas.m`. MATLAB R2024a is the Intel build: run with
`arch -x86_64 /Applications/MATLAB_R2024a.app/bin/matlab -batch ...`, outside the
sandbox. `runPowerSAS` calls `initpowersas` -> `savepath`, which fails harmlessly
(install dir not writable). Run one MATLAB session at a time: PowerSAS temp files
are timestamp-named and `clearAllTempFiles` deletes other sessions' files.
Public data `d_039_mod.m`: 10 six-order machines, no AVR/TG, 19 PQ loads; the
paper's 18 ZIP loads and 18 induction motors are NOT in the public file. PQ loads
become the constant-power part of PowerSAS's ZIP model (the PSAT vmin/z columns
are ignored), i.e. pure constant P/Q with no low-voltage switch.
Faults are specified on lines: `[line, position, r, x]`; bus 1 = line 1 (1-2) at 0.

Scripts and data: `scripts/powersas39/` (sweep39.m, isolate39.m, export39.m, a
shadow copy of solveAlgebraicNR.m that dumps the switch problem to JSON, and
clear_reinit.jl which solves it with Barq's solvers through `model!`).

Fig. 11 reproduced (fault bus 1, 0.5-0.75 s, dt = 0.01, |Zf| 0.010:0.001:0.020):
HE, ME-HE, RK4-HE normal everywhere; ME-NR and RK4-NR low-voltage or diverged for
|Zf| <= 0.019, normal at 0.020. TRAP crashes in public PowerSAS (`tgIdx` bug).
Isolation (method set per event): NR only at clearing -> V1 = 0; NR only at fault
application -> normal; NR only in the time steps -> normal. The failure is the
clearing re-init.

Cause: PowerSAS NR uses the power-mismatch residual
`conj(S) + I|V| + Ig conj(V) - (Y V) conj(V)`. Bus 1 has no load or generator,
so its row vanishes identically at V1 = 0. The "low-voltage solution" is V1 = 0
with every other bus satisfying KCL to 1e-8 and 56 pu of unbalanced current at
bus 1 (a phantom bolted fault). From the same inherited point, Newton on the
current-mismatch form converges in 4 iterations to the HE root (|Zf| = 0.010,
0.015, 0.019), even though V1 moves 0.37-0.54 -> 1.06-1.08. The Y- -> Y+
natural homotopy (50 steps) reaches the HE root in both forms. So the published
NR failure is a spurious root of the power formulation, not non-existence and
not intrinsic Newton difficulty. Our own `balance!` is power-form and has the
same artifact (seen earlier at V24 = 0).
Divergence ("X" in Fig. 11): all 13 diverged runs in our sweep already have
V1 = 0.000 at t = 0.76 s, i.e. they start with the spurious clearing root. The
power-form NR in later steps stays on V1 = 0, so bus 1 acts as a fault that
never cleared: ME-NR at |Zf| = 0.014 has the rotor-angle spread growing 46 ->
106 deg over 0.76-1.41 s, and PowerSAS stops at 1.41 s (flag -1, "Simulation
terminated earlier") when NR fails to converge. ME-HE on the same case: spread
peaks at 53 deg and decays, stable. Whether a run shows as L or X depends on
whether the separation gets large enough before t_end. No case found where the
clearing NR itself diverges (public data, no motors/ZIP). Script:
`scripts/powersas39/diverge39.m`.

## Polish 2383-bus case (Yao et al. first Polish case), 2026-09-29

Data `d_2383wp_mod2_ind_zip_syn.m` matches the paper: 327 six-order machines,
1827 ZIP loads, 1542 induction motors, Exc/Tg empty (controllers disarmed, the
paper's first Polish case). Fault at bus 1396 = line row 1674 (1396-1140) at
position 0, Zf = j0.01, 0.5-0.95 s. Export via `scripts/powersas39/export_case.m`
(the shadow solveAlgebraicNR now dumps Y as sparse triplets), results in
`scripts/powersas2383/fault1396/{he,appNR,clrNR}`; solve with
`julia --project=. scripts/powersas39/reinit_sparse.jl scripts/powersas2383/fault1396 clr 0.01`.

HE everywhere: V1396 1.152 -> 0.362 at application, 0.281 -> 0.896 at clearing,
run finishes. NR only at application: same root, fine. NR only at clearing:
PowerSAS NR converges (loop 6, flag 0, residual 5e-10) to V1396 = 0, then the
simulation terminates at t = 0.9506 s. Same spurious power-form root as the
39-bus case: all buses satisfy KCL except 1396 (40 pu mismatch). Barq power-form
Newton reproduces it (5 it.); current-form Newton from the same inherited point
converges in 3 iterations to the HE root (dist 8e-7), displacement 0.62 pu.
The paper's second Polish case (AVR/TG on, faults at 42 and 540 then 1396) cannot
be run: Exc/Tg data are not in the public file.
The V = 0 root only exists at a bus with no load/generator injection, so faults
at load buses are where a genuine Newton failure could show up.

## In-simulation re-initialization (fault ramped in over a window), 2026-09-29

Script: `scripts/insim_reinit_39.jl`. 39-bus, EMT lines, ZIP (0.7, 0.1, 0.2) with
`low_voltage=false` (the case with no frozen-current re-init solution), fault
bus 24, x_f = 0.015. From the pre-fault point, the fault susceptance is ramped
over a window T_c while implicit Euler integrates (h = T_c/200), then integration
continues with b = 1/x_f to t = 0.1 s. Ramps: linear in susceptance
(b = s/x_f) and the existing geometric one (x_eff = 1e10^(1-s) x_f^s).

- Linear: T_c = 1e-5, 1e-4, 5e-4 fail (MaxIter) at s = 0.105, 0.150, 0.645;
  T_c = 1, 2, 5 ms succeed (V24 at end of window 0.405, 0.527, 0.484).
- Geometric: T_c = 1e-5 ... 1e-3 fail at s = 0.915, 0.920, 0.940, 0.970;
  T_c = 2, 5 ms succeed (V24 dips to 0.137, 0.228 at end of window).
- Short-window failures sit at the frozen-current fold: geometric s = 0.915-0.92
  (frozen continuation folds at 0.918) and linear s = 0.105 are the same fault
  susceptance, b ~ 7 pu. Longer windows push the failure later, then remove it.
- After a successful window, Euler h = 5e-4 and 5e-5 and Trap h = 5e-4 all run to
  0.1 s; V24(0.1) = 0.476-0.477 (Euler) and 0.479-0.498 (Trap) for every T_c and
  shape, i.e. the end state does not depend on the window.
- Step-size control at fixed T_c: lin T_c = 1 ms succeeds with N = 50, 200, 2000
  (same V24 = 0.405). lin T_c = 0.5 ms fails at s = 0.645 (N = 200) and s = 0.396
  (N = 2000) but "succeeds" with N = 50: near the threshold the fine-step
  trajectory hits an impasse and coarse steps jump over it. geom T_c = 0.1 ms
  fails at s = 0.920 for N = 200 and 2000. The knob is T_c, not h.
Critical window ~0.5-1 ms (linear), 1-2 ms (geometric), in line with the
earlier estimate L*dI/dV ~ 0.1-0.6 ms for the line currents to redistribute.
Brown-style restart with original currents (`scripts/brown_restart_39.jl`):
y* = algebraic state at the end of the 1 ms linear window (line currents moved
by up to 15.8 pu during the window), restart from (x0, y*) with the pre-fault
delta, omega, line currents. Algebraic residual at (x0, y*) is 10 pu. Euler and
Trap: h = 2, 1, 0.5 ms succeed (first-step V24 = 0.445 / 0.374 / 0.268 Euler,
depends on h), h = 1e-4, 5e-5, 5e-6 all fail (MaxIter). y* only improves the
first step's Newton guess; solvability of that step needs h >~ 0.5 ms, the same
threshold as the window length. Smaller steps make it fail.

## Spurious V = 0 root at clearing: Newton vs homotopy on the same equations, 2026-09-29

- PowerSAS 39-bus (algebraic lines, constant PQ, no PQBRAK), fault bus 1:
  power-form Newton at clearing -> V1 = 0 (regular root, residual 1e-10), carried
  by every later step (phantom bolted fault), angle spread 46 -> 106 deg, NR
  failure and termination at 1.41 s. Y- -> Y+ natural homotopy (50 steps) on the
  same power-form equations from the same point -> physical root (matches HE to
  4e-6) for |Zf| = 0.010, 0.015, 0.019.
- Polish 2383, fault bus 1396: power-form Newton -> V1396 = 0 (5 it.), PowerSAS
  terminates at 0.9506 s. Power-form Y-homotopy (50 steps, 151 Newton it., 1028 s
  with dense ForwardDiff Jacobians) -> V1396 = 0.8964, matches HE to 8e-7.
  Current-form Newton -> same physical root in 3 it.
- Why it never recovers: at a bus with no injection the power-mismatch row is
  V_i * conj(sum of branch currents), so V_i = 0 is an exact root; its Jacobian is
  generically nonsingular, so the root persists as the other states evolve and
  the simulator follows a smooth non-physical branch.
- With PQBRAK, kP(V) -> 0 as V -> 0, so V = 0 also becomes a root at LOAD buses
  in the power form (seen earlier: V24 = 0 with pqbrak = 0.7). PQBRAK adds
  spurious roots rather than removing them.
- Barq's own 39-bus in RMS form (`scripts/spurious_root_39.jl`, classical
  machines, local power-form balance because src/models/bus.jl is being edited
  to a current-balance form): fault at zero-injection buses 1, 2, 5, 6, 9 with
  xf = 0.005-0.03, 0.25 s fault: Newton at clearing always reaches the same
  physical root as homotopy (3-7 it.). Buses 5 and 6 with xf <= 0.01 are
  transiently unstable (both roots lose synchronism). Milder fault-on stress
  than PowerSAS (angle spread ~20 deg at clearing vs 45 deg). Sweep of buses
  10-22 in progress.
  Full sweep done (buses 1, 2, 5, 6, 9, 10, 11, 13, 14, 17, 19, 22; xf = 0.005,
  0.01, 0.02): Newton at clearing matches the homotopy root in every case (3-8
  it.). No spurious V = 0 root in Barq's model. The post-clearing MaxIter cases
  fail identically from the Newton and homotopy roots (loss of synchronism, a
  stability outcome, not a re-init artifact).

## Capacitor bus model (C dv/dt = i_net - conj(S/V)), fault at bus 20, 2026-09-29

User switched to dynamic bus voltages: `build_mass_matrix` puts C_eq on
balance_d/q, `_setup_alg` treats delta, omega, line currents and bus voltages as
differential (re-init solves only fault and generator currents, linear), and
`balance!` is now the current form (verified: matches an independent
i_net - conj(S/V) reference exactly, pre-fault residual 5e-11). C_eq = B/2 line
charging + artificial 1e-6 on every bus; buses 12, 20, 31 carry load with only
the 1e-6. Driver: fault bus 20, x_f = 0.009.

`scripts/cap_model_fault20.jl`: re-init converges (2 it.) but the first step
fails at h = 5e-4, 5e-5, 5e-6 for the default kwargs and for ZIP without
low-voltage; pure Z and ZIP + PQBRAK 0.7 run at every h (V20 -> ~0.25). With
C20 raised to 1e-2 the collapse is visible: h = 5e-6 fails at step 171
(~0.86 ms, V20 = 0.045). `scripts/phasor_fault_exist.jl`: the phasor fault-on
steady state exists for all load models (V20 = 0.247 ZIP, Newton 4 it., natural
homotopy reaches lambda = 1). So the failure is transient: frozen line currents
at t = 0+, fast bus-20 voltage (tau ~ 1e-8 s) has no equilibrium, constant-P
current P/V blows up before the lines (~ms) redistribute. In-simulation ramp
of the fault admittance: T_c = 0.1, 0.5 ms fail (s = 0.065, 0.110); T_c = 1, 2,
5 ms succeed and the continuation at h = 5e-4 and 5e-5 reaches V20 = 0.245-0.252
at 50 ms, matching the phasor steady state.

## Event search for re-init failures: `scripts/reinit_search.jl`, 2026-09-29

Pipeline per (loading level, event): direct integration of the post-event model
from the inherited state at every h; if any h fails, Newton re-init and natural
homotopy re-init (event parameter s: 0 pre-event, 1 post-event); integrate from
each converged re-init; SAVED = direct fails at some h and a re-init recovers,
TARGET = direct fails at every h and the homotopy recovers. Events as homotopies
in s: fault apply/clear (lambda_fault = s or 1-s), line trip ((1-s) branch eq
- s i), generator trip (stator/swing blended to zero current, frozen rotor), load
step (P, Q x (1 + s step)). Families: line (N-1, no islanding), line2 (pairs at a
bus), topo (all lines at a bus except one), gen, load (+50%, +100%), fault (all
non-slack buses, xf = 0.001/0.005/0.02, apply and clear after 0.1 s), faultline
(fault at the sending bus cleared by tripping the line). Loading: base and 0.9,
0.97, 0.99 x LF_max (power-flow continuation, loads and dispatch scaled;
LF_max = 2.3125 for the 39-bus). Loads ZIP (0.7, 0.1, 0.2) + PQBRAK 0.7 char 1
(`pqbrak=` arg). Modes: rms (only delta, omega differential), emt (line currents
differential), cap (line currents + bus voltages). The script builds its own mass
matrix and passes the solvers an address whose keys make any version of
`_setup_alg` return the algebraic rows (asserted). BLAS pinned to one thread
(two default-threaded runs in parallel were ~60x slower). Output:
outputs/reinit_search/<mode>_<run>_<timestamp>/{config.txt, all_cases.csv,
saved_cases.csv}; reproduce a case with `match="<label>" lf=<tag>`.
Branch checks: emt bus 24 xf 0.015 with pqbrak=1e-6 -> C (homotopy folds at
0.92); cap bus 20 xf 0.009 pqbrak=1e-6 -> reinit_then_fails; rms bus 20 -> A.
First full run: 8 processes (rms and emt x 4 loading levels), logs in
outputs/reinit_search/log_<mode>_lf<tag>.txt.
First full 39-bus search results (PQBRAK 0.7 char 1, current-form balance!,
LF = 1.0 and 0.9/0.97/0.99 x LF_max = 2.3125, 426 events per level):
- rms: 1704 cases, no re-init failure. 12 cases fail only at h = 0.01 after
  15-19 steps (t = 0.15-0.19 s, fault clearing at buses 5, 6, 10, 11, 13 with
  xf <= 0.005); Newton re-init converges in 4-5 it. to the homotopy root and the
  run from it fails identically. Step-size failures mid-trajectory, not re-init.
  No fold anywhere: with PQBRAK the loads turn Z-like below 0.7 pu, so the
  post-event algebraic system always has a solution and Newton finds it.
- emt: 52 failures, all at small h (5e-6, some 5e-5) a few steps after line or
  generator trips / faults cleared by a trip. Newton re-init converges in 1-2 it.
  The first run's 27 "reinit_recovers" were an artifact of the post check not
  including h = 5e-6; with the corrected criterion (recover at every h where the
  direct run failed, same horizon) all are reinit_then_fails. Cause: the tripped
  line/gen current decays with time constant ~L (~7 us) under the status model,
  and small steps resolve that interruption while the other line currents stay
  frozen (inductive kick again).
- Conclusion: standard RMS model + PQBRAK + current (KCL) form gives a benign
  re-init for every standard event on the 39-bus up to 0.99 LF_max. Remaining
  routes: power-mismatch form (spurious V = 0 roots, which PQBRAK adds at load
  buses), larger stressed systems (2383 needs sparse Jacobians), IBR limiters.

Power-form RMS search (`form=power`, KCL physicality check on every accepted
state, 1704 cases): no re-init failure and no TARGET. The same 12 events that
failed at h = 0.01 in the current form (faults at zero-injection buses 5, 6, 10,
13 near LF_max, cleared by removal or by a line trip) now "succeed" at h = 0.01
but end on a spurious root (KCL mismatch 0.7-17 pu); Newton and homotopy re-init
both give the same physical root. `scripts/trace_case.jl` (lf 0.99, fault bus 6
xf 0.001 clear@0.1): the post-clearing system slips poles (angle spread 55 ->
323 deg), bus 7 sits near the electrical centre and its voltage physically goes
toward 0 at t = 0.16-0.19 s. Current form: singular Jacobian there (char 1 keeps
the constant-current part, whose direction is undefined at V = 0) -> MaxIter.
Power form at h = 0.01: jumps onto V7 = 0 and continues silently (KCL mismatch
12 -> 55 pu); at h = 1e-3, 1e-4: exception. Integration-phase phenomenon, not
re-init. Parser bug fixed: values containing "=" (e.g. match="... xf=0.001 ...")
were truncated at the second "=", so earlier match= runs selected extra cases.
Idea not yet tested: secondary events during the post-fault swing (line/gen trip
at t_clear + 0.05-0.2 s, as in protection/cascading), where the inherited state
is far from equilibrium with depressed voltages.

## Post-clearing low-voltage search: `scripts/clearing_lowv_search.jl`, 2026-09-29

Question: after clearing, does the solver settle on a low-voltage solution while a
high-voltage (HV) one exists, and which re-inits recover HV? Per clearing event
(fault removed / cleared by line trip; xf 0.001/0.005/0.02; 0.1 s fault; LF base,
0.9/0.97/0.99 x LF_max; PQBRAK 0.7 char 1): Newton re-init and the first direct
step vs. physical candidates from HC, flat (1∠0), flat magnitude with inherited
angles, and pre-fault voltages; low = differs from HV with min V lower by > 0.05.
2944 cases (184 per LF x 4 LF x {rms, emt} x {power, current}):
- rms power / rms current: 0 low-voltage outcomes. Newton always reaches HV.
- emt current: Newton never low; first direct step lands on a genuine physical
  low-voltage root in 153 cases (small h; frozen line currents make each load bus
  a scalar ZIP equation with a high and a low root).
- emt power: Newton -> spurious V = 0 at the faulted load bus in 294/736 cases
  (all LF, buses 3, 4, 7, 8, 12, 15, 16, 18, 20, 21, 23-29; 6-12 it., KCL mismatch
  15-68 pu); HC, flat, flat-angle and pre-fault all recover HV in 294/294. With
  PQBRAK, V = 0 is an exact power-form root at load buses. Integration from
  Newton's root throws (d|V|/dV undefined at V = 0 -> NaN Jacobian). Direct runs
  at small h (5e-5, 5e-6) land on the spurious root and continue for the 20-step
  check; `scripts/trace_lowv.jl` (lf 0.97, fault bus 3, xf 0.001, h = 5e-5):
  crashes before 0.02 s. Caveat: the HV root is the inductive kick, e.g. V3 = 3.5
  pu at clearing (trapped fault current), decaying to ~0.65 pu within 2 ms.

PowerSAS-like RMS runs (mode=rms: only delta, omega differential; form=power;
starts = hc, flat; outputs/clearing_lowv/rms_power_<cfg>_lf*):
| config | cases | Newton -> low | clean (stays low 0.3 s, HV run ok) | direct low at every h | HC / flat recover |
| ZIP + PQBRAK 0.7, 0.25 s | 736 | 5 (spurious V = 0) | 0 | 0 | 4 / 5 |
| ZIP no PQBRAK, 0.1 s | 537 | 109 (physical low-V) | 70 | 53 | 109 / 109 |
| ZIP no PQBRAK, 0.25 s | 331 | 42 | 26 | 21 | 42 / 42 |
| PQ no PQBRAK, 0.1 s | 36 | 4 | 4 | 4 | 4 / 4 |
| PQ no PQBRAK, 0.25 s | 29 | 2 | 2 | 1 | 1 / 2 |
(cases < 736 without PQBRAK: fault-on application/integration infeasible, skipped.)
Without PQBRAK, Newton from the inherited post-fault state converges to the lower
PV-branch root at the faulted bus (KCL satisfied to 1e-13, V ~ 0.006-0.2 while the
HV root has V ~ 0.66-1.0); integration from it stays low (e.g. fault bus 26,
xf 0.001, cleared by trip, base: V ~ 0.01 at 0.3 s vs 1.05 from HV) and the direct
run lands on it at h = 0.01, 1e-3 and 1e-4. This is the HE paper's "impractical
low-voltage solution" in our RMS model. HC and flat both recover HV in every case.
With PQBRAK 0.7 only 5 near-nose cases (faults at 8, 10, 28 cleared by a line trip,
LF >= 0.9): power-form Newton -> spurious V = 0 (KCL 9-14 pu), flat recovers 5/5,
HC 4/5; integration from Newton's root throws, direct runs fail or go spurious.
PQBRAK removes the lower PV root (load turns Z-like at low V), so the reviewer's
requirement removes most of the phenomenon.

Fault-duration sweep with stability filter (2026-09-30): RMS, power form, ZIP +
PQBRAK 0.7, LF 0.90/0.93/0.95/0.97/0.99, xf 0.001/0.003/0.005, tfault 0.10-0.20,
2 s hold; CLEAN = Newton's run stays at V < 0.05 on the LV bus while the HV run
keeps angle spread < 180 deg and KCL < 1e-3 throughout (columns N_*, HV_*, clean;
outputs/clearing_lowv/rms_power_dur_t*_lf*). 4600 cases: Newton -> low in 12, all
the same event (fault bus 28 xf 0.001 cleared by trip 28-29, tfault >= 0.15 s),
spurious V28 = 0 (KCL 10-12 pu) held for 2 s; HC and flat recover HV 12/12
(min V 0.55-0.71); direct runs spurious at every h in 11/12. No clean case: the HV
run loses synchronism in all 12 (spread 564-3700 deg or run fails). In this
classical-machine model the V = 0 trap needs large angle excursions at clearing,
which only occur beyond the critical clearing time. `scripts/trace_three_roots.jl`
traces Newton / HC / flat runs and saves CSVs. |V| is smoothed (veps=1e-12) in the
script-local load model; without it exact V = 0 (reached by underflow) gave NaN
Jacobians. Beware zsh: an unquoted "$var" holding several key=value pairs is passed
as ONE argument; scripts now reject malformed mode=/form=.

PowerSAS cases with PQBRAK (2026-09-30, `scripts/powersas39/pqbrak_reinit.jl`):
PowerSAS has no PQBRAK, so the exported clearing problems were re-solved with the
constant-power load part scaled by kP(|V|) (OpenIPSL char 1, PQBRAK 0.7); the
inherited pre-clearing state is still PowerSAS's (fault-on trajectory without
PQBRAK). The HE root keeps all |V| > 0.89, so it is unchanged (KCL 1e-4).
39-bus, |Zf| = 0.010/0.015/0.019: power-form Newton from the inherited state ->
V1 = 0 (6/7/25 it., KCL 46-56 pu); current-form Newton -> HE root (4 it.); flat
1∠0 power-form Newton -> HE root (5 it.); power-form HC (50 steps) -> HE root.
Polish 1396: power Newton -> V1396 = 0 (5 it., KCL 41 pu); current Newton (3 it.)
and flat power Newton (5 it.) -> HE root. PQBRAK does not change the outcome:
the faulted buses (1, 1396) carry no load, so their V = 0 root is independent of
the load model. Caveat: the fault-on trajectory was not recomputed with PQBRAK.
PowerSAS low-voltage flag: the PSAT PQ.con format carries Vmax, Vmin and an
"allow conversion to impedance" flag (d_039_mod: 1.2, 0.8, 1), but PowerSAS never
reads those columns (regulateSystemData uses P, Q, status; no low-voltage load
logic anywhere in internal/, util/, interface/). Empirical check on the exports:
Polish clearing state has 14 load buses below 0.8 pu (bus 1140 at 0.445, ...) still
exported with their full constant-power S, i.e. no conversion. 39-bus: all load
buses are >= 0.90 pu at clearing, so the flag would be inactive anyway. In both
cases the faulted bus (1, 1396) has no load, so its V = 0 root does not depend on
any load model.
ZIP instead of PQ on the PowerSAS 39-bus clearing problem (`pqbrak_reinit.jl ... zip=kz,ki,kp`;
constant-power loads re-split on the base |V_HE|, so the HE root stays a root; inherited
state still from PowerSAS's PQ trajectory): |Zf| = 0.010 and 0.015: power-form Newton
-> V1 = 0 for every mix tried, ZIP (0.7,0.1,0.2) with and without PQBRAK 0.7,
(0.4,0.3,0.3), pure I, pure Z; flat, HC and current-form Newton -> HE root.
|Zf| = 0.019: with ZIP + PQBRAK or pure Z, power-form Newton reaches the HE root
(13 / 12 it.; with constant PQ it went to V1 = 0 in 25 it.), so the load model only
moves the edge of the region. The Polish case already has ZIP loads (1827) + motors.

## Homotopy at fault clearing (λ: 1 -> 0), 2026-09-30

`solve_homotopy!` and `solve_adaptive_homotopy!` take `λ_start` (default 0.0), so
clearing re-init is `λ_start=1.0, λ_target=0.0` from the fault-on state at t_c.
`solve_homotopy!` now uses `range(λ_start, λ_target; length=...)` (always lands on
the target). Adaptive: direction `σ`, predictor uses the signed step, exact
landing on `λ_target`. Default 0 -> 1 behaviour unchanged.
Smoke test (39-bus, bus 24, x_f = 0.015, algebraic lines, `balance!` defaults,
Euler h = 5e-4 fault-on to 0.1 s): clearing Newton from the inherited point
DIVERGES (|V24| -> 1e159 in 39 it.); natural homotopy (Δλ = 0.05, 21 stages,
hardest stage λ = 0.95 with 25 it.) and adaptive (31 steps, 86 it.) both reach
V24 = 1.037. With the geometric x_eff the work sits in λ ∈ [0.75, 1].
Driver note: `main_convergence_failure_testing.jl` post-clearing run passes `u3`
(fault re-init point) instead of `u4 = sol_post.u[end]`.
Note: the V24 = 1.037 result above was with algebraic lines. On 2026-10-01 the tree
has EMT lines again (L in mass matrix, line currents out of `_setup_alg`) and the
keyword `balance!` is current form with I_load = conj(S_load/V) restored. Same
smoke test: fault-application homotopy folds (last converged λ = 0.95, |V24| =
0.19; ZIP with pqbrak 0.7e-6, no frozen-current solution); fault-on integration
from that iterate still succeeds (first Euler step h = 5e-4 lets currents move);
clearing Newton diverges, both homotopies 1 -> 0 converge to |V24| = 14.3 pu (the
frozen-current inductive kick). Algebraic Jacobian at the pre-fault point is
well conditioned (cond 94, σ_min 0.073): the 12 zero-injection buses all carry
real line charging (ωC_eq 0.07-0.79 pu), and the three smallest singular values
equal ωC_eq at buses 10, 11, 13. C_eq = 1e-6 only at buses 12, 20, 30-38, which
have loads or generators. So the EMT current form is index 1 on this case.

## case118, fault bus 12, x_f = 0.03: re-init converges but integration fails, 2026-10-01

EMT lines, current-form `balance!` defaults (ZIP 0.7/0.1/0.2, pqbrak 0.7e-6).
Frozen-current re-init converges (Newton 6 it. = homotopy, |V12| = 0.078), but
implicit Euler from it fails at any h: h = 1e-4 fails at step 1; h = 2e-5 ... 2e-7
all stop at t* ~ 7 µs. Cause: impasse point at bus 14 (load P = 0.14, Q = 0.01,
lines only to 12 and 15, ωC_eq = 0.035). On µs scales the inductive lines act as
current sources, so V14 solves I_net = I_ZIP(V); with k_z P0 V + k_p P0 / V the
deliverable current has a minimum 2 sqrt(k_z k_p) P0 + k_i P0 ~ 0.12 pu at
V = sqrt(k_p/k_z) = 0.535. Line 12-14 (L = 1.9e-4) drains bus 14 at ~5e3 pu/s,
I_net falls below the minimum in ~7 µs, V14 has fallen 0.98 -> 0.536 and
σ_min(g_y) drops with null vector on bus 14. Same for direct integration.
With pqbrak = 0.7, direct integration and integration from re-init both
succeed to 0.1 s (re-init 4 it.). Lesson: consistency at t+ is not enough; the
algebraic branch must stay regular along the flow. Scripts: scratchpad c118*.jl.
Smooth fault insertion (C∞ smoothstep, b = γ(t/Tc)/x_f, window h = Tc/200), same
case, pqbrak 0.7e-6 (scratchpad c118_ramp.jl): Tc = 1e-5 window ok (V14 min 0.537,
max |Δi_line| 0.10) but continuation fails at once (h = 5e-4, 5e-5); Tc = 1e-4,
1e-3, 3e-3, 1e-2 fail inside the window at γ = 0.26, 0.11, 0.39, 0.9995 (V14 min
0.55, 0.67, 0.60, 0.52). Non-monotone in Tc. No ramp works: the obstacle is the
fault-on DAE (bus-14 fold), not the switching discontinuity. Fast-ramp limit =
frozen-current homotopy in time (di/dτ = (Tc/L)(Δv - R i) -> 0).

## Letter case study 2b: Newton non-convergence at clearing (same model as the LV-root example), 2026-10-01

Question: besides the LV-root example (Newton converges to the wrong root), find an
event where Newton does not converge at all while the homotopy re-init recovers.
Source: the existing no-PQBRAK clearing searches (`outputs/clearing_lowv/
rms_power_zipNoBrak_*`) already contain ~40 such rows (`newton_conv=false`,
`newton_it=50`, direct runs fail at every h, `hc_conv=true`). Most are at 0.9-0.99
LF_max or with 0.25 s faults, where the HV run then loses synchronism (beyond
the CCT), so they are not clean. Two clean ones at BASE loading, 0.1 s fault:

- fault bus 16, xf = 0.001, cleared by tripping line 16-21 (chosen)
- fault bus 10, xf = 0.001, cleared by tripping line 10-11 (direct h = 0.01
  replicate converges in 50 it. to a genuine low-voltage root, V10 = 0.47, min V
  0.02; less clean)

Model: RMS (only delta, omega differential), power-mismatch bus equations, ZIP
(0.7, 0.1, 0.2), low_voltage=true with pqbrak = 1e-6 (i.e. no PQBRAK), char 1,
classical machines; identical to the LV-root case study.
Scripts: `scripts/trace_divergence.jl` (per-iteration CSVs, reuses reinit_search.jl
machinery) and `scripts/plot_divergence_case.jl` (figure). Output:
`outputs/clearing_lowv/divergence_fault_bus_16_xf_0_001_cleared_by_trip_16_21_lfbase_power/`,
figure `figures/review1/divergence_39_clearing{,_2panel}.pdf`.
Run: `julia --project=. scripts/trace_divergence.jl mode=rms form=power pqbrak=1e-6
tfault=0.1 lf=base "match=fault bus 16 xf=0.001 cleared by trip 16-21"`.

Bus 16 / trip 16-21 numbers:
- inherited (fault-on, t = 0.1 s) state: V16 = 0.064, angle spread 28 deg,
  post-clearing algebraic residual 63 (inf norm).
- Newton re-init from the inherited state: no convergence in 50 iterations. Iterates
  bounce (V16 = 0.06 -> 0.25 -> 0.06 -> 0.03 -> 0.12 -> 0.28 -> ...) and after ~12
  iterations settle in a wandering regime near V16 = 0.011 with |g| ~ 1e-2 and
  Newton step norm ~1 that never shrinks; cond(J) 1e4-1e5 there (1.7e3 at the
  start). Not a blow-up and not a root: bounded non-convergence near the lower
  PV branch / V = 0 region.
- Direct integration (implicit Euler, no re-init): MaxIter at h = 10, 1, 0.1 ms;
  the first-step Newton iterates reproduce the re-init behaviour almost exactly
  (the algebraic rows dominate), so step-size reduction does nothing.
- Homotopy in the clearing parameter (dlambda = 0.01): 234 Newton iterations
  total, V16 rises 0.06 -> 0.58 (lambda 0.1) -> 0.98 (0.2) -> 1.02 (1.0); root
  V16 = 1.022, min V = 0.983, KCL 2e-15. Flat start (1 angle 0) also converges to
  the same root in 5 iterations.
- Run from the homotopy root, h = 10 ms, 2 s: stable, max angle spread 61 deg, min
  V 0.86, V16(2 s) = 1.05. So the HV solution is the physical post-clearing state.
- Current-form (KCL) balance, same event (`form=current`, output dir `..._current`):
  Newton stalls the same way (V16 -> 0.013, |dy| ~ 3), so the stall is not a
  power-form artifact; the HC root exists there too (V16 = 0.96). The script's
  KCL check and the 2 s run are not meaningful in that mode (repo `balance!` is
  mid-edit), ignore them.
Interpretation for the letter: the inherited point (V16 = 0.06) lies in a region
where the post-clearing Newton map is ill-conditioned and non-contractive; the
displacement is not the issue (the flat start, much farther from the root, needs
5 iterations), the basin geometry is.

## 118-bus clearing search (same model as the 39-bus case studies), 2026-10-01

`scripts/clearing_lowv_search.jl mode=rms form=power pqbrak=1e-6 tfault=0.1 lf=base
data=cases/Fault_Cases/case118_gc.xlsx xf=<0.001|0.005|0.02> run=c118_xf<xf>`, three
processes, ~1 h each. Output `outputs/clearing_lowv/rms_power_c118_xf*_20261001_111010/`,
logs `log_c118_xf*.txt`. Summary script (scratch): counts over cases.csv.

| system | clearing events | direct ok at all h | Newton stall (50 it.) | Newton -> physical LV root (stays low) | HC -> HV | flat -> HV | HV run unstable |
|---|---|---|---|---|---|---|---|
| 39-bus base, xf 0.001/0.005/0.02 | 170 | 129 | 4 | 37 | 41/41 | 41/41 | 0 |
| 118-bus base, same xf | 745 | 409 | 31 | 236 | 263/267 | 266/267 | 0 |

118 by xf: 0.001 -> 168 events run (123 skipped: fault-on solution does not exist
or fault-on integration fails with constant-P loads), 6 stalls, 129 LV roots;
0.005 -> 286 events, 25 stalls, 106 LV roots; 0.02 -> 291 events, 0 stalls, 1 LV root.
Stalls cluster at buses 5, 6, 11, 17, 32, 54-56, 92 (fault-on V 0.02-0.10).
Spurious V = 0 roots: none (all LV roots satisfy KCL). No spurious column because
the faulted buses are load buses and PQBRAK is off.
Two deviations from "HC always recovers": (i) 4 events (fault bus 2 and 104,
xf 0.001, fault-on V ~ 0.004-0.01) where the natural homotopy follows the LV
branch continuously from the fault-on state to the same LV root Newton finds;
flat start reaches HV. Path-following follows the branch it starts on when the LV
branch connects without a fold. (ii) 1 event (fault bus 48 xf 0.001 cleared by
trip 48-49) recovered by HC only, flat start fails.
Note: the "HV run unstable" column is the hold check (0.3 s) from the search, not a
2 s stability check; the trace script does the longer run.

118-bus divergence trace (`scripts/trace_divergence.jl ... data=cases/Fault_Cases/case118_gc.xlsx
"match=fault bus 54 xf=0.005 clear@0.1"`), output
`outputs/clearing_lowv/divergence_fault_bus_54_xf_0_005_clear_0_1_lfbase_power/`,
figure `figures/review1/divergence_118_clearing.pdf`:
- inherited state V54 = 0.096, angle spread 96 deg (the 118 case has large static
  angle differences), post-clearing residual 18.
- Newton re-init: no convergence in 50 it., iterates wander V54 = 0.02-0.15 with
  min V down to 0.006, |g| ~ 1e-1, step norm ~1.
- direct integration: MaxIter at h = 1, 0.1, 0.01 ms (H_DIRECT in reinit_search.jl
  is now (1e-3, 1e-4, 1e-5) for rms).
- homotopy: 180 Newton it., root V54 = 0.931, min V 0.921; flat start 6 it., same
  root. 2 s run from it: stable, max spread 102 deg, min V 0.90, V54(2 s) = 0.955.
- bus 32 xf 0.005 (also a stall) is NOT clean: the HV run loses synchronism.

118-bus wrong-root case where ONLY the homotopy recovers (`trace_divergence.jl ...
data=cases/Fault_Cases/case118_gc.xlsx "match=fault bus 48 xf=0.001 cleared by trip 48-49" T=3.0`),
output `outputs/clearing_lowv/divergence_fault_bus_48_xf_0_001_cleared_by_trip_48_49_lfbase_power/`,
figure `figures/review1/wrongroot_118_bus48.pdf`:
- inherited state V48 = 0.011 (nearly bolted fault, 0.1 s), post-clearing residual 10.
- Newton re-init converges quadratically in 5 it. (residuals 10, 5e-2, 3e-4, 1e-8,
  1e-15) to a genuine low-voltage root: V48 = 0.0138, KCL 1e-15. Direct integration
  at h = 1 ms, 0.1 ms, 10 us converges at every step and runs on this root. Flat
  start (1 angle 0) converges in 8 it. to the SAME low-voltage root.
- Homotopy (dlambda 0.01, 240 Newton it.) reaches the HV root V48 = 0.945, min V 0.924.
- 3 s runs (h = 10 ms): from the HV root stable, min V 0.91, V48(3 s) = 0.966, spread
  <= 98.5 deg (static spread of this case ~96 deg). From the LV root also "stable":
  V48 stays at 0.013 for 3 s, KCL 1e-14, spread <= 99.8 deg. Both are regular
  trajectories; the simulator has no internal signal that the LV one is wrong.
- The trace script now also writes lv_run.csv / flat_run.csv when those roots differ
  from the homotopy root; plot_divergence_case.jl switches panel (c) to the
  time-domain comparison when lv_run.csv exists.

## Letter figures (three two-panel figures), 2026-10-01

Traces for the figures: `trace_divergence.jl ... h=0.001 T=0.3 tag=fig` (output dirs end in
`_fig`); figures from `scripts/plot_letter_figs.jl` -> `figures/review1/letter_fig{1,2,3}_*.pdf`.
The trace script now also writes faulton_run.csv (pre-fault 50 ms + fault-on segment, t = 0 at
clearing), direct_run.csv (no re-init, h_run, T), post_hc_step_h<h>.csv (Newton residuals of the
first implicit-Euler step from the homotopy root), lv_run.csv / flat_run.csv.
- Case 1 (user's existing LV-root example): 39-bus, fault bus 6, xf = 0.002, 0.1 s, fault removed.
  Inherited V6 = 0.095. Newton 14 it. -> LV root V6 = 0.097, min V 0.026 (bus 6 has no load; the
  low voltage sits at a neighbouring load bus), KCL 5e-15. Direct run at h = 1 ms follows the LV
  root, MaxIter at 0.16 s. HC (185 it.) and flat (5 it.) -> HV root V6 = 1.02, min V 0.98; 0.3 s
  run stable (spread <= 47 deg). First DAE step after HC re-init: 3 Newton iterations.
  `outputs/clearing_lowv/divergence_fault_bus_6_xf_0_002_clear_0_1_lfbase_power_fig/`
- Case 2: 39-bus bus 16 / trip 16-21 (see above); with H_DIRECT = (1e-3, 1e-4, 1e-5) now.
- Case 3: 118-bus bus 48 / trip 48-49 (see above).
