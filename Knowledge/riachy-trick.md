# Riachy's trick — an iPD with no derivative

`functions/mfc_riachy_transform.m` · block `mfc_fhat_riachy2_block` ·
tests `tests/test_riachy.m`

The decoupled second-order iPD of [[control-law-coupled-vs-decoupled]] needs
$\dot e$, and the only honest source of $\dot e$ is a differentiated
measurement. Riachy's trick removes it: the derivative feedback is turned into
an **integral of the measurement**, which costs one state and no noise gain.

## The rewrite

Start from the ultra-local model [[ultra-local-model]] and add $K_D\dot y$ to
both sides:

$$
\ddot y + K_D \dot y = F + K_D \dot y + \alpha u
$$

Define the auxiliary output, with $0 \le c < t$ free:

$$
Y(t) = y(t) + K_D \int_c^t y(\sigma)\,d\sigma
\qquad\Longrightarrow\qquad
\ddot Y = \ddot y + K_D \dot y
$$

and lump the derivative into the disturbance,
$\mathfrak F = F + K_D\dot y$. What is left is the **same second-order
ultra-local model, in $Y$**:

$$
\ddot Y = \mathfrak F + \alpha u
$$

So *any* second-order estimator — [[estimator-algebraic-2nd]] or
[[estimator-sliding-window]] — returns $\mathfrak F$ when you feed it $Y$
instead of $y$. Nothing about the estimator changes; only its input does.

## The control law

$$
u = -\,\frac{\mathfrak F_{\text{est}} - \ddot y^* - K_D \dot y^* + K_P e + K_I\!\int_c^t e}{\alpha}
$$

$K_I = 0$ gives the iPD. In this repo's sign convention ($e = y - sp_{filt}$,
and the command law *subtracts* `fb`, see [[block-library-signal-flow]]) that
maps onto the stock blocks as:

| term | where it goes |
|---|---|
| $\mathfrak F_{\text{est}}$ | `mfc_fhat_riachy2_block` → `F_hat` |
| $\ddot y^* + K_D\dot y^*$ | `ff`, built from the smoother's two derivative outputs with a Gain and a Sum |
| $K_P e + K_I\int e$ | `fb`, a stock Discrete PID with **D = 0** |

Substituting a perfect estimate gives

$$
\ddot e + K_D \dot e + K_P e + K_I\!\int e = 0
$$

— the poles of the ordinary iPD, obtained without differentiating anything.

**The two ways to get $K_D$ wrong**, both silent:

- leaving D on the PID applies $K_D$ **twice** (once inside $Y$, once outside);
- dropping $K_D\dot y^*$ from `ff` leaves the feedforward inconsistent with $Y$.
  This is invisible on a settled step and shows up only on a moving setpoint —
  `tests/test_riachy.m` pins it on a ramp ($2.7\times10^{-3}$ tracking error
  with the term, $2.8\times10^{-2}$ without).

## Discretization

Trapezoidal, not forward Euler:

$$
I[k] = I[k-1] + \tfrac{T_s}{2}\big(y[k] + y[k-1]\big), \qquad Y[k] = y[k] + K_D I[k]
$$

The estimators downstream weight their window by powers of $t$, so a
half-sample bias in $Y$ is not free. With $I[-1] = y[-1] = 0$ the first sample
counts half, so $I$ equals the continuous integral shifted by $T_s/2$ — harmless,
because $c$ is free and a constant offset in $I$ does not reach $\ddot Y$.

## Which estimator, and why it is a parameter

`mfc_fhat_riachy2_block` selects algebraic vs sliding window on the **mask**,
which is a deliberate exception to the repo's "a structural choice is a
different block" rule: both choices have the same ports, the same wiring and
estimate the same $\mathfrak F$. Only the numerics differ.

The numerics differ in a way that matters here more than elsewhere, though,
because $\mathfrak F$ contains $K_D \dot y$ and therefore **moves fast**:

- **Sliding window.** A fixed window returns the average over $[t-T_w, t]$, so
  for a ramping $\mathfrak F$ it returns exactly $\mathfrak F(t - T_w/2)$ — a
  clean half-window lag, verified to $10^{-3}$ in `tests/test_riachy.m`. Finite
  memory also means the unbounded $\int y$ never accumulates weight.
- **Algebraic.** The window grows from $t=0$ and weights by $t^2$, so it lags
  more, by an amount that depends on $t$ — and $\int y$ is unbounded, so with a
  non-zero steady state $Y$ ramps forever underneath those growing weights.

Prefer the sliding window for this trick unless you have a reason not to.

## NxN

The rewrite is linear, so it survives being made vector-valued verbatim. With
$y, u, F \in \mathbb R^n$ and **square** $K_D, \alpha \in \mathbb R^{n\times n}$:

$$
Y = y + K_D \int_c^t y \qquad\Longrightarrow\qquad \ddot Y = \mathfrak F + \alpha u,
\qquad \mathfrak F = F + K_D \dot y
$$

and the command law becomes the linear solve of [[block-library-signal-flow]]'s
matrix pair, $u = \alpha^{-1}(-\mathfrak F_{\text{est}} + \ddot y^* + K_D \dot y^*
- K_P e - K_I\!\int e)$, closing $\ddot e + K_D \dot e + K_P e = 0$ as a matrix
polynomial. $K_D$ need not be diagonal: a full $K_D$ folds *cross-channel*
derivative feedback into $Y$, which is the one thing this port buys over running
$n$ SISO blocks side by side.

Both estimators act channel-wise through **one shared window** — the algebraic
denominator $t^2$ is scalar, and the sliding window's taps are scalar quadrature
weights. That is the same asymmetry (per-channel numerator, shared denominator)
that `tests/test_golden_mimo.m` pins for the plain algebraic estimator.

Block `mfc_fhat_riachy2_mimo_block`, feeding `mfc_command_mimo_block`; tests
`tests/test_riachy.m` §4 (numerics, and a 2×2 loop on a cross-coupled plant)
and `tests/test_estimators.m` §8 (ports, masks, reduction to the SISO block).
`mfc_fhat_sliding_window` grew a column per channel in its buffers to make this
work; the scalar path is byte-for-byte unchanged.

## There is no first-order Riachy

Not missing — **meaningless**. The trick removes $K_D\dot e$ from an iPD by
moving it inside $\mathfrak F$. A first-order ultra-local model
($\dot y = F + \alpha u$) closes with an iP: a single pole at $-K_P$, no
derivative feedback anywhere, nothing to remove. There is no $K_D$ to fold into
$Y$, and folding $K_P$ instead would give $\dot Y = \dot y + K_P y$, i.e. the
*coupled* first-order estimator ([[control-law-coupled-vs-decoupled]]) by a
longer route. So the grid stops at second order, at every width, on purpose.

## The no-regression contract

`tests/test_riachy.m` checks that the trick *behaves*; `tests/test_golden_riachy.m`
pins the actual numbers, over four traces: `riachy_win`, `riachy_alg`,
`riachy_mimo_win`, `riachy_mimo_alg` (the last two with a **non-diagonal**
$K_D$, so cross-channel derivative folding is in the contract). It also
re-derives $Y$ from the logged measurement at every sample, which is what would
catch a forward Euler slipping in where the trapezoid belongs.

Two tuning facts fell out of capturing those traces, and both are the same
fact: $\mathfrak F$ carries $K_D\dot y$ and therefore **moves fast**, so any
estimator lag is expensive. At $T_s = 10$ ms with a 4 rad/s target, $T_w = 0.4$ s
diverges and $T_w = 0.1$ s tracks to $10^{-6}$; the algebraic variant needs its
num/den smoother halved (10 → 5 samples) for the same reason. Neither is a
property of the plant — it is the price named in *What it costs* below, showing
up as a tuning constraint.

## The noise result, which is arguably the bigger win

Removing $K_D \dot e$ removes the loop's dominant noise path. In a
representative decoupled loop ($\sigma = 10^{-4}$, $T_s = 1$ ms, $K_d = 16$) the
derivative term carried ten times more command noise than the estimator did;
Riachy's trick replaces it with an *integral* of the measurement, which
attenuates noise instead of amplifying it. Same plant, same closed-loop poles,
no filter and no added lag: $\operatorname{std}(u)$ fell from 1.152 to 0.120.

That is a better argument for the trick than the one it is usually sold on. See
[[noise-and-filter-placement]] for the diagnosis that decides when it applies,
and for why filtering the derivative term instead is the expensive option.

## What it costs

Against the ordinary decoupled iPD on the same plant and the same gains
($p = 8$, $K_P = p^2$, $K_D = 2p$, $T_s = 1$ ms), the Riachy loop tracks the
same step with a peak deviation of $5.6\times10^{-2}$ — the estimator's lag on
the fast-moving $K_D\dot y$ term. That is the price: the derivative is not
free, it is moved inside $\mathfrak F$ where the estimator's bandwidth applies
to it. What you buy is that no differentiator ever touches the measurement.

## See also

- [[ultra-local-model]] — the model being rewritten
- [[control-law-coupled-vs-decoupled]] — the iPD this replaces, and the third
  option (folding $K_P$, $K_D$ into a coupled estimator instead)
- [[estimator-sliding-window]] — the recommended estimator for $Y$
- [[estimator-algebraic-2nd]] — the alternative, and its growing window
- [[estimator-free-derivatives]] — the other way to get derivatives without
  differentiating: read them out of the estimate itself
- [[block-library-signal-flow]] — ports, sign conventions, the unit delay
- [[noise-and-filter-placement]] — where to filter once the D term is gone
