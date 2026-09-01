---
tags:
  - retro
created: 2026-08-31
subject: Make the algebraic estimators' smoother settings runtime-tunable
---
# Retro — tunable `est_filter_window` / `est_hold_time`

Moved both properties out of `properties (Nontunable)` in the algebraic
F-estimator blocks so a downstream GA can sweep them under Simulink fast restart
instead of recompiling per candidate.

## Toolbox inventory (verbatim, `detect_matlab_toolboxes`)

MATLAB Version: 26.1.0.3312084 (R2026a) Update 4, macOS 26.6.2.

```
MATLAB, Simulink, Aerospace Blockset, Aerospace Toolbox, Computer Vision Toolbox,
Control System Toolbox, Curve Fitting Toolbox, DSP System Toolbox, Deep Learning
Toolbox, Embedded Coder, Global Optimization Toolbox, Image Processing Toolbox,
MATLAB Coder, MATLAB Compiler, MATLAB Compiler SDK, MATLAB MCP Server Toolbox
(0.3.0), MATLAB Report Generator, MATLAB Test, Navigation Toolbox, Optimization
Toolbox, Reinforcement Learning Toolbox, Robust Control Toolbox, Sensor Fusion and
Tracking Toolbox, Signal Processing Toolbox, Simscape, Simscape Electrical,
Simulink Check, Simulink Coder, Simulink Compiler, Simulink Control Design,
Simulink Coverage, Simulink Design Optimization, Simulink Design Verifier,
Simulink Desktop Real-Time, Simulink Fault Analyzer, Simulink Test, Stateflow,
Symbolic Math Toolbox, System Composer, System Identification Toolbox, UAV
Toolbox, Vehicle Dynamics Blockset
```

All at 26.1 (R2026a) unless noted.

- **Global Optimization Toolbox: present.** `ga` is available.
- **Parallel Computing Toolbox: absent.** Confirmed by its absence from the list
  above, not inherited from the docs. No `parpool`, no `parsim`, no `UseParallel`
  — a GA campaign here is serial.

## Changed — the property move (11 blocks)

`mfc_fhat_alg{1,2}_{coupled,decoupled}_block`, their four `_mimo_` twins,
`mfc_fhat_riachy2_block`, `mfc_fhat_riachy2_mimo_block`,
`mfc_fhat_decoupled_dev_block`.

Confirmed in each before editing: neither property appears in `setupImpl`,
`resetImpl` or `getDiscreteStateSpecificationImpl`. The eight pure algebraic
blocks have no `setupImpl` at all and size their states from `n` alone. The
riachy2 and dev blocks size buffers from a separate `window_samples`, which was
left `Nontunable`.

## Deliberately left alone

- **`mfc_siso_core`, `mfc_mimo_core` — contradicts the brief.** Both use
  `est_filter_window` in `bufferLength`, which feeds `resetImpl` *and*
  `getDiscreteStateSpecificationImpl`; on the sliding-window path it is
  `n_buf = est_filter_window + 1`. It genuinely sizes a discrete state there, so
  it stays `Nontunable`. The brief's "if they expose them" would have broken
  them. A GA that wants to sweep the window on a core must compose the loop from
  stage blocks instead.
- `window_samples` on `mfc_fhat_window*`, and `n` / `Ts` everywhere — untouched,
  as instructed.

So the grid is not uniform, but the split is principled: tunable exactly where
the value reaches only `stepImpl`. `Knowledge/codegen-constraints.md` §6 was
rewritten to say this; it previously claimed the property was Nontunable
everywhere because it sizes buffers, which is only true of the cores and the
window blocks.

## `loadObjectImpl` — checked, and I could not reproduce the defect

The eight pure algebraic blocks need nothing: no private properties, no
`setupImpl`, so there is no un-serialized state to lose.

The audit did turn up two blocks with a private `kernel` built in `setupImpl` and
**no** `saveObjectImpl`/`loadObjectImpl` — `mfc_fhat_decoupled_dev_block` and
`mfc_fhat_window_mimo_block`, the latter being the direct twin of the block that
carries the fix and its measured-2026-08-21 comment. Added the same pair to both.

**But: the failure does not reproduce on R2026a.** Without the fix, the dev block
on its sliding-window setting ran three consecutive fast-restart sims and a
3-element array `sim()` with no error and identical output every run. So the
addition rests on structural analogy to a failure the author measured, not on one
observed this session. Cost is zero and the twin asymmetry was a trap; if a later
release makes it demonstrably dead, it can go.

## Fast restart — demonstrated

Trivial model, `mfc_fhat_alg2_decoupled_block`, `FastRestart` on:

```
run 1  W=10   0.70 s   status after run = compiled
run 2  W=60   0.14 s   status after run = compiled
run 3  W=10   0.03 s
max|F(W=10) - F(W=60)|      = 1.18075   the value bites
max|F(W=10) - F(W=10 again)| = 0        repeatable
```

Counterfactual on the pre-change tree: the same `set_param` was rejected with
`Simulink:Parameters:NoChangeWhenInFastRestart`. Post-change, `Ts` is still
rejected with that identifier, so fast restart is genuinely engaged and it is the
declaration that changed, not the enforcement.

## Tests

`test_estimators` (MATLAB, the block layer — the layer this change touches):
**86 passed, 0 failed.**

**Octave is not installed in this sandbox**, so the six Octave suites could not be
run as specified. Ran them under MATLAB instead (`test_golden.m` documents that it
also runs there). `octave_sanity` and `test_riachy` pass fully. `test_golden`,
`test_composed`, `test_golden_mimo` and `test_golden_riachy` report failures at
1e-13..1e-11 on the 2nd-order algebraic variants.

**These are pre-existing and not caused by this change.** Verified by stashing the
edits and re-running: the baseline output is identical digit for digit, same
variants, same magnitudes. The suites exercise `functions/` and the golden CSVs,
neither of which was touched. The goldens were captured under Octave and the
contract is bit-identical with no tolerance, so the residual is a
MATLAB-vs-Octave floating-point difference in the 2nd-order algebraic path.
Nothing was re-blessed. **The bit-identical acceptance gate has therefore not
been run on its own terms** — it needs an Octave install, and someone should
close that out.

## Downstream

No block interface changed: no port, name, default or mask entry moved. Setting
these from a block dialog or `set_param` behaves as before, so no downstream
model needs an edit. `library/mfc_lib.mdl` was not regenerated — no block was
added, removed or renamed.

Codegen is a real, if narrow, behaviour change: in the eleven blocks the two
values become runtime parameters rather than compile-time constants. Nothing in
this repo generates code from those blocks today, and the cores — the documented
codegen entry point in §"Verifying" — are unchanged. Not re-verified with
`codegen`; that check was not run this session.

## See also

- [[codegen-constraints]] — §6, rewritten by this change
- [[estimator-algebraic-2nd]], [[iir-smoother]] — the math that makes it safe
