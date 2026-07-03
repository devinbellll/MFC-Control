# Tests

Regression tests for the MFC SISO controllers in `../functions`.

## `test_estimators.m` — MATLAB (tests the real System objects)

Drives the actual `matlab.System` objects through `step()` in closed loop and
asserts tracking, steady-state `F` accuracy/sign, and a sign-flip regression.

```matlab
>> cd <repo root>      % the folder containing functions/ and tests/
>> test_estimators
```

Covers:
- `mfc_siso_non_algebraic` (`use_first_order=false`) — step tracking +
  steady-state `F` vs analytic value.
- **Sign-flip regression** — an open-loop-*unstable* plant (`ÿ = 4y + 2u`) with
  weak PD so that correct `F`-cancellation is the only thing keeping it bounded.
  The historical `-60/T^5` bug diverges here; the fixed `+60/T^5` stays bounded.
- `mfc_siso_non_algebraic` (`use_first_order=true`) — step tracking + steady-state `F`.
- `mfc_siso_core` (`use_first_order=false`, folded) / `mfc_siso_decoupled` — still
  track (validated estimators).
- `mfc_siso_core` (`use_first_order=true`) — the new algebraic 1st-order F
  estimator (raw F identification + explicit `kp` feedback) — step tracking +
  steady-state `F`.

Prints `N passed, N failed` and `error()`s if anything fails.

> Note: this harness needs `matlab.System`, so it only runs in MATLAB, not
> Octave. Its closed-loop logic mirrors the Octave check below (which was run and
> is green); the MATLAB file itself has not been executed in CI yet — run it once
> in MATLAB to confirm against your toolbox version.

## `octave_sanity.m` — Octave or MATLAB (tests the estimator math)

MATLAB-free guard on the estimator kernels and signs. Octave cannot instantiate
`matlab.System`, so this re-implements the `use_first_order=true/false`
branches of `mfc_siso_non_algebraic`'s `stepImpl` and checks it against plants
with an analytically known `F`.

```bash
octave --no-gui -q tests/octave_sanity.m   # exit code 0 = pass
```

Covers constant-`F` recovery, slow-sinusoid tracking, and **positive-sign
guards** that fail on the two bugs that were fixed:
- 2nd-order prefactor sign (`-60/T^5` returns `-F`),
- 1st-order `u`-term sign (`- u_kernel` returns `F + extra*u`).

> If you change an estimator formula in `functions/mfc_siso_non_algebraic.m`,
> mirror the change in `est1`/`est2` inside `octave_sanity.m` (they are
> intentionally a standalone copy of the math).
