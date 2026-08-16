# Tests

Regression tests for the MFC SISO controller in `../functions` and `../blocks`.

Two tiers, because the numerics and the Simulink object layer fail in different
ways:

| File | Runs in | Covers |
|---|---|---|
| `octave_sanity.m` | Octave + MATLAB | estimator/smoother math vs analytic answers (SISO §1-6, MIMO §7-8) |
| `test_golden.m` | Octave + MATLAB | the controller still reproduces reference traces |
| `test_composed.m` | Octave + MATLAB | the pipeline decomposes faithfully into stages |
| `test_golden_mimo.m` | Octave + MATLAB | the 2x2 matrix-alpha pair reproduces its reference traces, and reduces to the SISO trace when alpha is diagonal |
| `test_estimators.m` | **MATLAB only** | the `matlab.System` objects: ports, reset, masks |

```bash
octave --no-gui -q --path tests --path functions --eval octave_sanity
octave --no-gui -q --path tests --path functions --eval test_golden
octave --no-gui -q --path tests --path functions --eval test_composed
octave --no-gui -q --path tests --path functions --eval test_golden_mimo
```

```matlab
>> setup
>> test_estimators
```

## `test_golden.m` — the no-regression contract

`golden/*.csv` holds the exact per-sample outputs (`u, F_hat, sp_filt, err,
u_raw, valid`) of all six supported variants on the `val_mfc` drone plant,
captured from the tree **before** the block decomposition. `test_golden`
re-runs `mfc_golden_trace` and demands a **bit-identical** match.

No tolerance is used, deliberately. Moving code between files does not change
floating-point results, so any difference at all means something real changed.
If this ever seems to need a tolerance, the change is bigger than advertised.

Regenerating the reference is an explicit decision, not a way to make the test
pass:

```bash
octave --no-gui -q --path tests --path functions --eval golden_capture
```

`1st_coupled_alg` **diverges** on this plant, and that is part of the contract —
a true double integrator cannot be stabilized from a first-order ultra-local
model, which has no derivative room to fold. It must keep diverging the same
way. See [`../Knowledge/control-law-coupled-vs-decoupled.md`](../Knowledge/control-law-coupled-vs-decoupled.md).

## `test_golden_mimo.m` — the MIMO (matrix-alpha) pair's no-regression contract

`golden/mimo_diag.csv` and `golden/mimo_cross.csv` hold the closed-loop traces
of the 2x2 decoupled estimator + command pair (`mfc_fhat_alg2_decoupled_mimo_block`
/ `mfc_command_mimo_block`), captured the same way as the SISO golden traces —
see `mfc_golden_trace_mimo.m`. `test_golden_mimo` re-runs it and demands a
**bit-identical** match, then does one more thing the SISO test cannot: it
cross-checks `mimo_diag` (`alpha = eye(2)`, two identical channels) against
`golden/2nd_decoupled_alg.csv` itself, to a `1e-9` tolerance (not bit-exact,
since the matrix solve and the vector estimator arithmetic take a different —
if mathematically equivalent — path than the scalar SISO code). That is the
direct evidence that the matrix solve in `mfc_command_mimo` reduces to the
SISO result when the input gain is diagonal.

`mimo_cross` uses an off-diagonal `alpha` and two **different** plants, so
`F_hat` differs per channel while the estimator's `den_raw = t^2` stays
scalar and shared — see `octave_sanity.m` §7 for a direct structural pin of
that asymmetry, independent of the golden capture.

Regenerating the reference is an explicit decision, exactly as for `test_golden.m`:

```bash
octave --no-gui -q --path tests --path functions --eval golden_capture_mimo
```

## `test_composed.m` — is the decomposition honest?

Rebuilds the controller from the individual stages (`mfc_siso.ref_traj`, an
`mfc_fhat_*` estimator, `mfc_siso.feedback`, `mfc_siso.command`,
`mfc_siso.limit`) wired by hand, and checks it against the same golden traces.
It duplicates the wiring on purpose — that duplication *is* the test.

This proves the **stage math**. It does not exercise Simulink port wiring or
sample-time propagation; `test_estimators.m` section 5 and
`examples/val_mfc_composed.m` do that in MATLAB.

## `octave_sanity.m` — the math against closed-form answers

Feeds each estimator the exact sampled signals of a plant built to have a chosen
`F`, and demands it recovers that `F` with the right magnitude **and sign**.

This file used to carry its own hand-written copy of the estimator math, because
the algorithms lived inside `matlab.System` classes Octave cannot instantiate.
They are plain functions now, so it calls the real code — there is no second copy
to keep in sync.

Tolerances are not uniform, and the differences are meaningful:

- 1st-order sliding window: `1e-6`. The integrand is `y_kernel` (linear) ×
  `y` (linear) = quadratic, which Simpson integrates **exactly**.
- 2nd-order sliding window: `1e-3`. The integrand is quartic and Simpson is exact
  only through cubics, so a small quadrature residual is correct behaviour.
- algebraic estimators: `2e-2`. Backward-difference discretization, `O(Ts)`.

**Sign guards.** Two historical bugs are pinned explicitly, because both show up
as a sign flip rather than a magnitude error:

- 2nd-order prefactor `-60/Tw^5` returns `-F`
- 1st-order negated `u` kernel returns `F + (extra)*u`

**MIMO (§7-8).** Two more properties, pinned directly rather than only via the
closed-loop golden traces: the 2nd-order decoupled estimator's `den_raw`/`den_filt`
stay scalar while `num_raw`/`num_filt` are per-channel, a diagonal `alpha` does
not cross-couple channels (bit-identical to running the scalar estimator twice),
and `mfc_command_mimo`'s matrix solve reduces to `mfc_command` when `alpha` is
diagonal but is *not* equivalent to elementwise division once it isn't.

## `test_estimators.m` — the object layer (MATLAB)

Drives the real System objects through `step()`. Beyond tracking and
steady-state `F`, it checks three things only the object layer can get wrong:

- the undefined **coupled + sliding-window** combination is *rejected*, not
  silently mis-run
- the **stage blocks wired by hand reproduce `mfc_siso_core`** to `1e-12`,
  coupled and decoupled (the explicit PID stands in for the stock Simulink one)
- `reset()` genuinely restores the initial state (same inputs → same outputs)

> This file was stale for several commits — it referenced
> `mfc_siso_non_algebraic`, `mfc_siso_decoupled` and `mfc_siso_run`, all deleted
> in `228cf5d`. It has been rewritten against the current API. It still needs one
> run in MATLAB to confirm against your toolbox version.
