# The block library and signal flow

Nine blocks in `blocks/`. One is the assembled controller; eight are the stages
it is made of. All of them are thin `matlab.System` wrappers — **no math lives
in `blocks/`**, only ports, state and masks. The math is in `functions/`.

Generate the Simulink library with `build_mfc_lib` (run once in MATLAB); it
builds `library/mfc_lib.mdl` from these classes, so the classes stay the single
source of truth.

## The blocks

| Block | Stage | In → Out |
|---|---|---|
| `mfc_siso_core` | *all of them* | `y_sp, y_m, t` (+`u_applied`,`alpha`) → `u, F_hat, sp_filt, err, u_raw, F_valid` |
| `mfc_smoother_block` | 1, 3, post | `x` → `x_filt` (+`dot_x`,`ddot_x`) |
| `mfc_fhat_alg1_block` | 2 | `z, u_prev, t` (+`alpha`) → `F_hat, valid` (+`num_raw`,`den_raw`) |
| `mfc_fhat_alg2_block` | 2 | same, plus `a_fold` |
| `mfc_fhat_window_block` | 2 | `y, u_prev, t` (+`alpha`) → `F_hat, valid` (+`integral`) |
| `mfc_fhat_divide_block` | 3 | `num, den, t` → `F_hat, valid` |
| `mfc_feedback_block` | 3 | `err` (+`freeze`) → `fb, int_err, dot_err` |
| `mfc_command_block` | 4 | `F_hat, ff, fb` (+`alpha`) → `u_raw` |
| `mfc_command_filter_block` | 5 | `u_raw` (+`u_prev`) → `u, sat` |

`mfc_fhat_divide_block` and `mfc_command_block` are **stateless** and can run at
any rate. Every other block fixes its own discrete rate from `Ts` and must not
inherit — the recursions assume exactly one advance per sample.

## The assembled loop

```
              ┌──────────────────────────────────────────────┐
              │                                              │
 setpoint ──► [smoother] ──┬── sp_filt ──►(−)──► err ──► [feedback] ──► fb ──┐
   (role 1)                └── ddot_x ──────────── ff ──────────────┐        │
                                                                    ▼        ▼
        y ──────────────┬──────────────► [F-hat estimator] ──► F_hat ──► [command]
                        │                        ▲                            │
                        │                        │                            ▼
                        │                     u_prev                   [command filter]
                        │                        │                       │        │
                        │                        │                       u       sat
                        │                   ┌────┴───────┐               │        │
                        └── plant ◄─────────┤ unit delay │◄──────────────┘        │
                                            └────────────┘                        │
                                                  ▲     ┌────────────┐            │
                                                  └─────┤ unit delay │◄───────────┘
                                        to feedback.freeze └────────────┘
```

## The two wiring rules

Everything `mfc_siso_core` hides from you comes down to these.

### 1. `u_prev` must be routed back explicitly

The estimator needs the command that was actually applied over the **last**
sample — the current one does not exist yet when it runs.

`mfc_siso_core` keeps this internally as `state.u_km1`. An assembled loop needs
a real Unit Delay from the command filter's `u` back to the estimator's
`u_prev`. Same signal; it is simply visible now.

If a real actuator or an external limiter modifies the command, feed *that*
signal instead — via `mfc_command_filter_block`'s `u_prev` input and the
estimator's `u_prev` port. **Give both blocks the same signal**, or they
disagree about what the plant was actually driven with.

### 2. Anti-windup is a delayed handshake

`mfc_command_filter_block` emits `sat` when the clamp bit this sample.
`mfc_feedback_block` takes an optional `freeze` input that discards this
sample's integration.

That path **also needs a Unit Delay** — without one, `feedback → command →
command filter → feedback` is an algebraic loop and Simulink will refuse to
solve it.

> **Consequence, by design:** a composed loop freezes the integrator **one
> sample later** than `mfc_siso_core`, which does both in the same sample
> (`mfc_siso.step` applies `frozen` to the integral it just computed). The
> difference is one sample of extra windup at the moment of saturation.
> `mfc_siso_core` remains the reference implementation.

With saturation off there is no difference at all — which is why
`tests/test_composed.m` matches bit-for-bit.

## Matching flags across blocks

Split into separate blocks, consistency is now yours to maintain. Three pairs
must agree:

| must match | why |
|---|---|
| `mfc_feedback_block.coupled` ↔ your estimator's folding | coupled folds $K_p$ into $\hat F$; leaving `coupled = false` applies it **twice** |
| `alpha` on the estimator ↔ on `mfc_command_block` | otherwise $\hat F$ and the inversion assume different models |
| `Ts` on every stateful block | the stages must step in lockstep |

And the feedforward: `ff` is `ddot_x` for a second-order model, `dot_x` for
first order. That is the only place model order enters the command path.

## Dissecting the F-hat filtering

The algebraic estimators smooth numerator and denominator **before** dividing,
which is mathematically essential (see [[estimator-algebraic-2nd]]). So a
standalone "F-hat filter" block after the divide is a *different animal*. The
library gives you both:

```
default (self-contained):
    [alg2] ─────────────────────────────────────► F_hat

dissected (expose_raw + internal_filter off):
    [alg2] ─ num_raw ─► [smoother W] ─┐
           ─ den_raw ─► [smoother W] ─┴─► [divide + hold] ─► F_hat

post-filter (a genuinely different operation):
    F_hat ─► [smoother] ─► F_hat_filt
```

Both smoothers in the dissected chain must use the **same** $W$. The
reconstruction is verified bit-identical to the internal path, for both orders,
in `tests/`. `build_mfc_lib` ships this arrangement as a ready-made
"Dissected Estimator" subsystem.

The post-filter is the right tool for [[estimator-sliding-window]], which has no
internal smoothing to dissect.

## Is the decomposition honest?

Yes, and it is tested three ways:

| test | runs in | proves |
|---|---|---|
| `tests/test_composed.m` | Octave | the stage **maths** compose back into `mfc_siso.step`, all 6 variants, bit-identical |
| `tests/test_estimators.m` §5 | MATLAB | the stage **blocks** reproduce `mfc_siso_core` to `1e-12` |
| `examples/val_mfc_composed.m` | MATLAB | the same, as a runnable overlay plot |

## See also

- [[ultra-local-model]] — what the pipeline is computing
- [[control-law-coupled-vs-decoupled]] — which flags to match
- [[iir-smoother]] — the smoother's four roles
- [[codegen-constraints]] — why the blocks are written the way they are
