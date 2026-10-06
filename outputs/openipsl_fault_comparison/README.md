# OpenIPSL low-voltage load comparison

Implemented both low-voltage characteristics in `src/models/bus.jl`, preserving the previous function inside a block comment. The network uses the original power-balance equations and differential line currents. No homotopy or state reset is used in these tests.

## Configuration

- Case: `cases/Fault_Cases/ieee39_fault.xlsx`; fault bus 20; fault reactance 0.015 p.u.
- Pre-event integration: implicit Euler, 500 microseconds, through 1 ms. All initial non-slack bus voltages exceed 0.7 p.u.
- Fault-on test: one step, using implicit Euler or trapezoidal at 500, 50, and 5 microseconds.
- Newton: fresh Jacobian every iteration, correction tolerance 1e-7, maximum 200 iterations. The repository solver normally allows 30 iterations.
- Post-solve verification: full first-step residual infinity norm <= 1e-6; equivalent current mismatch checked separately.
- Loads: 100% PQ and the paper's 70/10/20 ZIP mixture. Low-voltage breakpoint: 0.7 p.u.

## Results for the paper's ZIP mixture

| Method | Step (microseconds) | Unmodified ZIP | OpenIPSL characteristic 1 | OpenIPSL characteristic 2 |
|---|---:|---|---|---|
| Euler | 500 | Converged (6 iterations) | Converged (5 iterations) | Converged (5 iterations) |
| Euler | 50 | Failed within 200 iterations | Converged (7 iterations) | Converged (7 iterations) |
| Euler | 5 | Failed within 200 iterations | Converged (9 iterations) | Converged (9 iterations) |
| Trap | 500 | Failed within 200 iterations | Converged (6 iterations) | Converged (6 iterations) |
| Trap | 50 | Failed within 200 iterations | Converged (8 iterations) | Converged (8 iterations) |
| Trap | 5 | Failed within 200 iterations | Converged (9 iterations) | Converged (9 iterations) |

All 12 modified-ZIP solves also pass the current-balance check (maximum current mismatch 3.301e-13 p.u.). The minimum bus voltage at the converged first steps ranges from 0.096845 to 0.241102 p.u. These checks verify a numerical root, not its uniqueness or physical branch selection.

## Qualifications for 100% PQ

- Unmodified PQ fails in all six combinations within 200 iterations.
- Characteristic 1 returns a small power residual in all six combinations, but minimum voltages approach zero and equivalent current mismatches are about 7-13 p.u. These are not acceptable current-balanced network solutions, despite the solver success code.
- Characteristic 2 returns small power and current residuals in all six combinations. However, the lowest-voltage load multiplier is negative (about -0.0075 to -0.0098), inherited from the upstream trigonometric fit at very low voltages. These results cannot be treated as physically validated passive-load solutions.
- Characteristic 2 with PQ, Euler, and 50 microseconds takes 39 iterations: it exceeds the normal 30-iteration budget even though it converges with 200 allowed.

## Implementation choices and source behavior

- `balance!` defaults to `zip=(0.0, 0.0, 1.0)`, preserving the working file's previous load mixture, with `characteristic=1`, `pqbrak=0.7`, and `low_voltage=true`. To test the paper mixture, use `zip=(0.7, 0.1, 0.2)`. The comparison script forwards these options through its DAE right-hand side.
- OpenIPSL's voltage-dependent multipliers are applied to the existing ZIP components normalized at initial voltage. Its separate default P/Q transfer fractions are not imposed, so the requested load mixture is retained.
- Characteristic 1 preserves the upstream strict comparisons exactly: at voltage equal to zero or half of PQBRAK it falls through to a constant-power multiplier of one. This is documented in the code; no silent boundary correction was made.
- Source: https://github.com/OpenIPSL/OpenIPSL/blob/master/OpenIPSL/Electrical/Loads/PSSE/BaseClasses/baseLoad.mo
- Load equation: https://github.com/OpenIPSL/OpenIPSL/blob/master/OpenIPSL/Electrical/Loads/PSSE/Load.mo

## Reproduction and limits

Run `julia --project=. --startup-file=no scripts/test_openipsl_fault.jl` from the repository root.
The script passes 13 checks, including comparison with the preserved old balance function and directional finite-difference checks of automatic derivatives. It writes `summary.csv` and `newton_traces.csv` here. `elapsed_s` includes compilation for some runs and is not a controlled performance benchmark.
This experiment tests first-step Newton convergence only, at one bus and one fault reactance. It does not establish successful full fault-on trajectories or continuation recovery.

## Input fingerprints (SHA-256)

- `src/models/bus.jl`: `f013411525d6d5dabc6c85fbce4812e657787dd124bc4579bdc435a22cc8f8ed`
- `scripts/test_openipsl_fault.jl`: `9fd10759dcbe595fc5588b63f156e97c3299f68fcb288a96105331e9207a5f00`
- `src/models/line.jl`: `6983f4d1f3aca20e104cc2a55222f09ec40cfcf9787f713f167338a69eaae270`
- `src/dynamic_sim.jl`: `2da00b1373d315a4e54f46e0248a4849cfc71cfab9933eae82682eb673608370`
- `src/models/fault.jl`: `a94a57a65908527a473e4723b98c31e075078fc0e0df29a0be4b7e8e3db5d86d`
- `cases/Fault_Cases/ieee39_fault.xlsx`: `916523cee3ac889a05b89c6e6dc8b1bddfbcb9c3b96f6224690b0e17bce91211`
