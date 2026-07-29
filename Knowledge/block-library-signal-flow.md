# The block library and signal flow

Eight blocks in `blocks/`. One is the assembled controller; seven are the pieces
you build a loop from. All of them are thin `matlab.System` wrappers — **no math
lives in `blocks/`**, only ports, state and masks. The math is in `functions/`.

Generate the Simulink library with `build_mfc_lib` (run once in MATLAB); it
builds `library/mfc_lib.mdl` from these classes, so the classes stay the single
source of truth.

## The design rule

**A structural choice is a different block, not a parameter.** Coupled vs
decoupled changes what you wire in, what gains exist, and where the feedback
lives — so it is two blocks with different port names, not a fold coefficient
you have to remember to set. Model order is likewise separate blocks.

**Anything Simulink already does well is not wrapped.** There is no feedback
block (stock Discrete PID), no command-filter block (stock Discrete Filter), and
no divide block.

## The blocks

| Block | In → Out |
|---|---|
| `mfc_siso_core` | `y_sp, y_m, t` (+`u_applied`,`alpha`) → `u, F_hat, sp_filt, err, u_raw` |
| `mfc_smoother_block` | `x` → `x_filt` (+`dot_x`,`ddot_x`) |
| `mfc_fhat_alg1_decoupled_block` | `y, u_prev, t` (+`alpha`) → `F_hat` |
| `mfc_fhat_alg2_decoupled_block` | `y, u_prev, t` (+`alpha`) → `F_hat` |
| `mfc_fhat_window_block` | `y, u_prev, t` (+`alpha`) → `F_hat` |
| `mfc_fhat_alg1_coupled_block` | `err, u_prev, t` (+`alpha`) → `F_hat` |
| `mfc_fhat_alg2_coupled_block` | `err, u_prev, t` (+`alpha`) → `F_hat` |
| `mfc_command_block` | `F_hat, ff, fb` (+`alpha`) → `u` |

The estimator's **first input port name tells you the structure**: `y` means
decoupled (true-plant $\hat F$, you supply the PID), `err` means coupled (the
poles are folded in, `fb` is Ground). There is no coupled sliding-window block
and there cannot be one — see [[estimator-sliding-window]].

Gains appear only on the block that actually uses them: `Kp`/`Kd` on the coupled
2nd-order estimator, `Kp` only on the coupled 1st-order one (no derivative room
at first order), neither on any decoupled block.

`mfc_command_block` is **stateless** and can run at any rate. Every other block
fixes its own discrete rate from `Ts` and must not inherit — the recursions
assume exactly one advance per sample.

## The assembled loop

Decoupled (the PID is yours):

```
 setpoint ──► [smoother] ──┬── sp_filt ──►(−)──► err ──► [Discrete PID] ──► fb ──┐
                           └── ddot_x ──────────── ff ─────────────────┐         │
                                                                       ▼         ▼
        y ──────────────┬──────────► [F-hat estimator: y] ──► F_hat ──► [command] ──► u
                        │                     ▲                                    │
                        │                  u_prev                                  │
                        │                     │   ┌────────────┐                   │
                        └── plant ◄───────────┴───┤ unit delay │◄──────────────────┘
                                                  └────────────┘
```

Coupled (`Kp`, `Kd` live inside `F_hat`; `fb` is Ground unless you want `Ki`):

```
 setpoint ──► [smoother] ──┬── sp_filt ──►(−)──► err ──┐   [Ground] ──► fb ──┐
                           └── ddot_x ────── ff ───┐   │                     │
                                                   ▼   ▼                     ▼
        y ──────────────┬──────────► [F-hat estimator: err] ──► F_hat ──► [command] ──► u
                        │                     ▲                                     │
                        └── plant ◄───────────┴──────── [unit delay] ◄──────────────┘
```

## The two wiring rules

### 1. `u_prev` must be routed back explicitly

The estimator needs the command that was actually applied over the **last**
sample — the current one does not exist yet when it runs.

`mfc_siso_core` keeps this internally as `state.u_km1`. An assembled loop needs a
real Unit Delay from `u` back to the estimator's `u_prev`. Same signal; it is
simply visible now.

If a real actuator or an external limiter modifies the command, feed *that*
signal to `u_prev` instead — it is what the plant was actually driven with.

### 2. The error sign is `y - sp_filt`

Measurement minus **filtered** setpoint, matching `mfc_siso.step`, and the
command law **subtracts** `fb`. This matters twice:

- the coupled estimators fold $-K_p$ (and $-K_d$) against this convention;
  feeding `sp - y` inverts the folded poles and the loop diverges;
- your PID must be driven by the same error, not its negation.

Use the raw `y_sp` rather than `sp_filt` and you lose the feedforward/feedback
consistency that makes the step response clean.

## No anti-windup in a composed loop

`mfc_command_block`'s clamp is a plain limiter. It owns no integrator, so it has
nothing to freeze, and there is no `sat`/`freeze` handshake anywhere in the
library any more.

With `Ki = 0` — the usual MFC case — this costs nothing. If you need integral
action **and** a real actuator limit, either use a PID block with its own
back-calculation anti-windup, or use `mfc_siso_core`, which freezes its integral
in the same sample the clamp bites (`mfc_siso.step` applies `frozen` to the
integral it just computed).

## Matching parameters across blocks

Split into separate blocks, consistency is yours to maintain. What still has to
agree:

| must match | why |
|---|---|
| `alpha` on the estimator ↔ on `mfc_command_block` | otherwise $\hat F$ and the inversion assume different models |
| `Ts` on every stateful block | the stages must step in lockstep |
| `ff` ↔ the estimator's order | `ddot_x` for a 2nd-order block, `dot_x` for 1st |

The old "does `coupled` match your estimator?" trap is gone: with a coupled
estimator there is no P or D term to double, because you wire `fb` to Ground.

## Smoothing $\hat F$

The algebraic estimators smooth numerator and denominator **before** dividing,
with one shared window — mathematically essential (see
[[estimator-algebraic-2nd]]), and therefore internal. It is not something you
wire up, and the raw ingredients are not exposed.

A smoother placed **after** an estimator is a genuinely different operation and
is the right tool for [[estimator-sliding-window]], which has no internal
smoothing at all:

```
F_hat ─► [mfc_smoother_block] ─► F_hat_filt
```

## Is the decomposition honest?

Yes, and it is tested three ways:

| test | runs in | proves |
|---|---|---|
| `tests/test_composed.m` | Octave | the stage **maths** compose back into `mfc_siso.step`, all 6 variants, bit-identical |
| `tests/test_estimators.m` §5 | MATLAB | the stage **blocks** reproduce `mfc_siso_core` to `1e-12`, coupled and decoupled |
| `examples/val_mfc_composed.m` | MATLAB | the same, as a runnable overlay plot |

`mfc_siso.feedback` and `mfc_siso.limit` still exist in `functions/` — no block
wraps them, but `mfc_siso.step` (and so `mfc_siso_core`) does, and the tests use
`mfc_siso.feedback` as the stand-in for the stock PID so the comparison stays
exact.

## See also

- [[ultra-local-model]] — what the pipeline is computing
- [[control-law-coupled-vs-decoupled]] — what the two estimator families mean
- [[iir-smoother]] — the smoother's roles
- [[codegen-constraints]] — why the blocks are written the way they are
