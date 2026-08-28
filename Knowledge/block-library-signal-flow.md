# The block library and signal flow

Nineteen blocks in `blocks/`. Two are assembled controllers — `mfc_siso_core`
and `mfc_mimo_core` — one is an estimator bench (`mfc_fhat_decoupled_dev_block`),
and the other sixteen are the pieces you build a loop from: eight SISO (scalar
`alpha`) and eight n-channel (square-matrix `alpha`).
All of them are thin `matlab.System` wrappers — **no math lives in `blocks/`**,
only ports, state and masks. The math is in `functions/`.

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

**The benches are the deliberate exception.** `mfc_siso_core` and
`mfc_mimo_core` put the whole variant grid on one mask; `mfc_fhat_decoupled_dev_block`
does the same one stage down, for the estimator alone. The cores are the
assembled reference implementation, and both are where you find out which
variant you want — a sweep becomes a parameter change instead of a rewiring
job. Ship the specific blocks, whose port names state the structure.

The estimator bench can only reach the **decoupled, non-Riachy** grid, and that
is what makes it possible rather than what limits it: those four variants (times
any `n`) share one wiring diagram, so they are interchangeable behind fixed
ports. Coupled estimators take `err` and carry the gains they fold; Riachy's
$\hat F$ means something else and needs its own `ff`. Neither can hide behind an
unchanged port list.

## The blocks

| Block | In → Out |
|---|---|
| `mfc_siso_core` | `y_sp, y_m, t` (+`u_applied`,`alpha`) → `u, F_hat, sp_filt, err, u_raw` |
| `mfc_smoother_block` | `x` → `x_filt` (+`dot_x`,`ddot_x`) |
| `mfc_fhat_alg1_decoupled_block` | `y, u_prev, t` (+`alpha`) → `F_hat` |
| `mfc_fhat_alg2_decoupled_block` | `y, u_prev, t` (+`alpha`) → `F_hat` |
| `mfc_fhat_window_block` | `y, u_prev, t` (+`alpha`) → `F_hat` |
| `mfc_fhat_riachy2_block` | `y, u_prev, t` (+`alpha`) → `F_hat` (+`Y`) |
| `mfc_fhat_alg1_coupled_block` | `err, u_prev, t` (+`alpha`) → `F_hat` |
| `mfc_fhat_alg2_coupled_block` | `err, u_prev, t` (+`alpha`) → `F_hat` |
| `mfc_command_block` | `F_hat, ff, fb` (+`alpha`) → `u` |
| `mfc_mimo_core` | `y_sp, y_m, t` (n-by-1) (+`u_applied`,`alpha`) → `u, F_hat, sp_filt, err, u_raw` (n-by-1) |
| `mfc_fhat_decoupled_dev_block` | `y, u_prev, t` (n-by-1) (+`alpha`) → `F_hat` (n-by-1) — order/estimator/`n` on the mask |
| `mfc_smoother_mimo_block` | `x` (n-by-1) → `x_filt` (+`dot_x`,`ddot_x`), all n-by-1 |
| `mfc_fhat_alg1_decoupled_mimo_block` | `y, u_prev, t` (n-by-1) (+`alpha`, n-by-n) → `F_hat` (n-by-1) |
| `mfc_fhat_alg2_decoupled_mimo_block` | `y, u_prev, t` (n-by-1) (+`alpha`, n-by-n) → `F_hat` (n-by-1) |
| `mfc_fhat_window_mimo_block` | `y, u_prev, t` (n-by-1) (+`alpha`, n-by-n) → `F_hat` (n-by-1) |
| `mfc_fhat_riachy2_mimo_block` | `y, u_prev, t` (n-by-1) (+`alpha`) → `F_hat` (+`Y`), n-by-1 |
| `mfc_fhat_alg1_coupled_mimo_block` | `err, u_prev, t` (n-by-1) (+`alpha`) → `F_hat` (n-by-1) |
| `mfc_fhat_alg2_coupled_mimo_block` | `err, u_prev, t` (n-by-1) (+`alpha`) → `F_hat` (n-by-1) |
| `mfc_command_mimo_block` | `F_hat, ff, fb` (n-by-1) (+`alpha`, n-by-n) → `u` (n-by-1) |

The estimator's **first input port name tells you the structure**: `y` means
decoupled (true-plant $\hat F$, you supply the PID), `err` means coupled (the
poles are folded in, `fb` is Ground). There is no coupled sliding-window block
and there cannot be one — see [[estimator-sliding-window]].

`mfc_fhat_riachy2_block` is decoupled too, but its `F_hat` is
$\mathfrak F = F + K_D\dot y$, not $F$: it estimates from $Y = y + K_D\int y$,
so the derivative feedback arrives inside the estimate. Its PID therefore needs
**D = 0** and its `ff` needs $\ddot{sp} + K_D\,\dot{sp}$ — see
[[riachy-trick]].

The bottom block is the **n-channel surface**: the same math on n-by-1 vector
signals with a square n-by-n `alpha` (and, for the coupled estimators, square
`Kp`/`Kd`). It mirrors the SISO surface block for block, with one absence that
is structural rather than an omission — see
[The n-channel (matrix-alpha) surface](#the-n-channel-matrix-alpha-surface)
below.

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

## The n-channel (matrix-alpha) surface

Every SISO block has a vector-valued twin: same math, same wiring diagrams as
above, but `y`, `u_prev`, `F_hat`, `ff`, `fb` and `u` are all n-by-1, and
`alpha` is a square, invertible n-by-n matrix instead of a scalar:

```
       F_hat, ff, fb, u : n-by-1
       alpha            : n-by-n, invertible

  estimator:  ddot_y = F + alpha*u              (vector y, u, F)
  command:    alpha * u = -F_hat + ff - fb   ->  u = alpha \ (-F_hat + ff - fb)
```

The channels are cross-coupled **through `alpha` alone** in the decoupled
blocks — and additionally through matrix `Kp`/`Kd` in the coupled ones, and
through matrix `Kd` in `mfc_fhat_riachy2_mimo_block`. `t` stays a single scalar clock shared by the whole
vector (not one per channel): the estimator's `den_raw = t^2` is scalar, and
only the numerator (`num_raw`, built from `z` and `alpha*u_prev`) is
per-element. That is what lets one growing window serve every channel, and it
is the property [`tests/test_golden_mimo.m`](../tests/test_golden_mimo.m) and
`tests/octave_sanity.m` §7 pin directly.

**Closing a MIMO loop** follows the same two wiring rules as the SISO
decoupled loop below, plus one more:

- `err = y - sp_filt` and `u_prev` through a real Unit Delay, exactly as SISO
  — `mfc_smoother_block` and a per-channel PID (or a vector Discrete PID) can
  stay as-is, since they don't need to know about `alpha`.
- The grid is **complete except where it is undefined**: 1st and 2nd order,
  coupled and decoupled, algebraic and sliding window, plus Riachy's trick.
  The one gap is coupled + sliding window, which does not exist at any width
  — an FIR window has nowhere to fold poles, see [[estimator-sliding-window]].
  There is no `mfc_fhat_riachy1_*` either, and there never will be: a
  first-order iP has no derivative feedback to remove, so the trick has
  nothing to do.
- The **coupled** estimators fold MATRIX `Kp` (and `Kd` at 2nd order), so the
  closed loop is a matrix polynomial and off-diagonal entries fold
  *cross-channel* P/D action into $\hat F$. That is the one thing they do
  that n SISO coupled blocks side by side cannot. They are also the
  tuning-fragile corner: with the channels' own dynamics unequal the same
  gains that work per-channel can diverge, because nothing explicit is left
  outside $\hat F$ to absorb the mismatch. The decoupled variants tolerate
  it. `tests/golden/mimo_2nd_coupled_alg.csv` is captured on a plant with
  identical per-channel dynamics for exactly this reason.
- `mfc_fhat_window_mimo_block` does **no internal smoothing** (nothing FIR
  does) — follow it with `mfc_smoother_mimo_block` if $\hat F$ is noisy.
- **The same `alpha` matrix must be fed to both blocks.** This matters more
  than in the scalar case: an estimator/command mismatch does not just bias
  the model gain, it can silently swap which physical channels are coupled.
- `alpha` must be well-conditioned. It is inverted every command-block
  sample (`alpha \ (...)`, a genuine linear solve, not an elementwise
  divide); a near-singular `alpha` blows up `u` the same way a near-zero
  scalar `alpha` would in the SISO case.
- **Sanity-check a new `alpha` against the SISO block first.** When `alpha`
  is diagonal, both blocks reduce exactly to running independent SISO loops
  with each diagonal entry as that channel's scalar `alpha` — that is how
  `tests/test_golden_mimo.m` cross-checks the MIMO pair against
  `tests/golden/2nd_decoupled_alg.csv` rather than only against its own
  capture. If a new plant's diagonal-`alpha` case doesn't match its SISO
  equivalent, look at the MIMO wiring before the plant.

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
| `tests/test_composed.m` | Octave | the stage **maths** compose back into `mfc_siso.step`, all 6 SISO variants, bit-identical |
| `tests/test_estimators.m` §5 | MATLAB | the stage **blocks** reproduce `mfc_siso_core` to `1e-12`, coupled and decoupled |
| `examples/val_mfc_composed.m` | MATLAB | the same, as a runnable overlay plot |

The n-channel surface has its own no-regression contract,
`tests/test_golden_mimo.m` (traces `mimo_diag`, `mimo_cross`,
`mimo_1st_decoupled_alg`, `mimo_2nd_decoupled_win`, `mimo_2nd_coupled_alg`),
and Riachy's trick has `tests/test_golden_riachy.m` (`riachy_win`,
`riachy_alg`, `riachy_mimo_win`, `riachy_mimo_alg`, the last two with a
non-diagonal $K_D$). `tests/test_estimators.m` §10 pins the other half of the
claim: `mfc_mimo_core` at n = 1 reproduces `mfc_siso_core` bit for bit, and at
n = 2 reproduces a loop hand-composed from the individual blocks.

`mfc_siso.feedback` and `mfc_siso.limit` still exist in `functions/` — no block
wraps them, but `mfc_siso.step` (and so `mfc_siso_core`) does, and the tests use
`mfc_siso.feedback` as the stand-in for the stock PID so the comparison stays
exact.

## See also

- [[ultra-local-model]] — what the pipeline is computing
- [[control-law-coupled-vs-decoupled]] — what the two estimator families mean
- [[iir-smoother]] — the smoother's roles
- [[codegen-constraints]] — why the blocks are written the way they are
