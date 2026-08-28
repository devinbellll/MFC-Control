# MFC_SISO

Reusable Model-Free Control SISO blocks and estimators (MATLAB/Simulink). Run
`setup.m` to put the paths on the MATLAB path.

`blocks/` and `library/` are a shared surface consumed by downstream projects, so
a change here can ripple outward — keep the block interfaces and the golden tests
honest.

## Layout

- `functions/` — **all the math.** `mfc_siso` (the stage functions and `step`),
  the three estimator kernels, `mfc_iir_smoother`.
- `blocks/` — thin `matlab.System` wrappers: ports, state, masks. **No math
  lives here.** This is the shared, downstream-consumed surface.
- `library/` — `build_mfc_lib.m` generates `mfc_lib.mdl` from the classes in
  `blocks/`, so the classes stay the single source of truth. Run it once in
  MATLAB and commit the regenerated `.mdl`.
- `examples/`, `models/`, `scripts/` — MATLAB/Simulink.
- `tests/` — including `tests/golden/`; run before submitting a block change.
- `Knowledge/` — algorithm/design notes: `block-library-signal-flow.md`,
  `ultra-local-model.md`, `control-law-coupled-vs-decoupled.md`,
  `codegen-constraints.md`, the `estimator-*` notes, and `Plans/`.

## The shape of the library

`mfc_siso_core` is the all-in-one controller and the reference implementation;
`mfc_mimo_core` is its n-channel twin and the whole-loop development bench — the
whole non-Riachy grid (order x structure x estimator x channel count) on one
mask, so a variant sweep is a parameter sweep. At `n = 1` it takes the scalar
code path and reproduces `mfc_siso_core` bit for bit.
`mfc_fhat_decoupled_dev_block` is the same idea for the ESTIMATOR alone: order,
estimator and `n` on the mask, `F_hat` out, ports that never change — reaching
the decoupled non-Riachy grid only, which is precisely the set that shares one
wiring diagram. Everything else composes a loop by hand out of a smoother, an
estimator and the command block.

Every SISO block has an n-channel twin (`*_mimo_block`): same math on n-by-1
signals with a square n-by-n `alpha`, and square `Kp`/`Kd` where the estimator
folds them. The grid is complete except for coupled + sliding window, which is
undefined at any width (an FIR window has nowhere to fold poles), and
first-order Riachy, which is meaningless (an iP has no derivative to remove).

Two rules explain what is and is not a block:

1. **A structural choice is a different block, not a parameter.** Coupled vs
   decoupled and 1st vs 2nd order are four separate estimator classes. A coupled
   estimator takes `err` and carries `Kp` (and `Kd` at 2nd order) on its own
   mask; a decoupled one takes `y` and has no gain parameters. The first input
   port name tells you the structure. Do not reintroduce `a_fold`/`b_fold`-style
   knobs.
   The one exception is `mfc_fhat_riachy2_block` (and its vector twin
   `mfc_fhat_riachy2_mimo_block`), where the estimator kind
   (algebraic or sliding window) *is* a mask parameter: both choices have the
   same ports, the same wiring and the same $\mathfrak F$, so it is numerics,
   not structure.
2. **Anything Simulink already ships is not wrapped.** The explicit feedback is
   a stock Discrete PID into `mfc_command_block`'s `fb` input (or Ground, when
   the estimator is coupled); output filtering is a stock Discrete Filter.

Consequences worth knowing before "fixing" something:

- A composed loop has **no anti-windup**. `mfc_command_block`'s clamp is a plain
  limiter that owns no integrator. Only `mfc_siso_core` freezes its integral in
  the same sample the clamp bites.
- The algebraic estimators smooth numerator and denominator internally with one
  shared window. That is mathematically essential, not a convenience — it is why
  there is no divide block and why the raw ingredients are not exposed.
- The n-channel estimators share ONE integration window across the vector: the
  numerator is per-element, the denominator (`t`, `t^2`, or the FIR taps) is
  scalar. Do not give a channel its own window.
- The coupled n-channel estimators are the tuning-fragile corner: with unequal
  per-channel dynamics, gains that work per-channel can diverge, because nothing
  explicit is left outside `F_hat` to absorb the mismatch.
- Estimator blocks emit `F_hat` only; the old `valid` / `num_raw` / `den_raw` /
  `integral` debug ports are gone.

## Wiring invariants

- `err = y - sp_filt` — measurement minus **filtered** setpoint, and the command
  law **subtracts** `fb`. The coupled estimators fold $-K_p$, $-K_d$ against this
  convention; feeding `sp - y` inverts the folded poles and the loop diverges.
- `u_prev` needs a real Unit Delay, and must be the command that actually reached
  the plant.
- `alpha` must match between the estimator and `mfc_command_block`; `Ts` must
  match across every stateful block.
- `ff` is `ddot_sp` for a 2nd-order model, `dot_sp` for 1st — the only place
  model order enters the command law.

## Testing

```bash
octave --no-gui -q --path tests --path functions --eval octave_sanity
octave --no-gui -q --path tests --path functions --eval test_golden
octave --no-gui -q --path tests --path functions --eval test_composed
octave --no-gui -q --path tests --path functions --eval test_riachy
octave --no-gui -q --path tests --path functions --eval test_golden_mimo
octave --no-gui -q --path tests --path functions --eval test_golden_riachy
```

```matlab
>> setup
>> test_estimators     % MATLAB only: ports, masks, reset, block composition
```

`tests/test_riachy.m` covers Riachy's trick (the transform, the estimate, and
a composed loop against the ordinary iPD), §4 the NxN port; `tests/test_estimators.m`
§7 covers its block layer and §8 the NxN one, §9 the rest of the n-channel
estimator surface, §10 `mfc_mimo_core` and §11 `mfc_fhat_decoupled_dev_block`
(every setting against the specific block it stands in for). `tests/test_golden_riachy.m` is
Riachy's no-regression contract, SISO and NxN.

`tests/golden/*.csv` is the no-regression contract: six SISO variants
(`test_golden`), five n-channel ones (`test_golden_mimo`) and four Riachy ones
(`test_golden_riachy`). The capture scripts take variant names --
`golden_capture_mimo('name')` -- so ADDING a variant never rewrites the
existing files, which is how a regression gets blessed by accident. `test_composed` and `test_estimators` §5 both assert that a
hand-assembled loop reproduces `mfc_siso_core` exactly (`mfc_siso.feedback`
stands in for the stock PID so the comparison stays bit-exact). If you change a
stage, those are the tests that catch it.

## Rules

- Changing a block interface is a cross-project event — verify the golden tests
  and flag downstream impact rather than silently reshaping a block.
- Math goes in `functions/`, never in `blocks/`. If a block starts computing,
  it is in the wrong file.
- Width is not a fork in the math: `mfc_siso.step` runs both widths, and the
  kernels are vector-safe. If you find yourself writing a `*_mimo` copy of a
  function in `functions/`, stop.
- Regenerate `library/mfc_lib.mdl` with `build_mfc_lib` after any block add,
  remove or rename.
- Obsidian conventions: `[[wikilinks]]`, flat YAML frontmatter. Filename stem is
  the join key — never rename half a set.
- Prefer the minimal change with the load-bearing property.
