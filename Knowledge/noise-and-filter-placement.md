# Noise, and where to put a filter

`functions/mfc_iir_smoother.m` · `mfc_siso.window_kernel` · blocks
`mfc_smoother_block`, `mfc_fhat_riachy2_block`

A composed loop offers four places to put a low-pass: on $y$ before the
estimator, inside the estimator, on $\hat F$ after it, and on $u$. They are
**not interchangeable**, and two of them are actively wrong. This note is the
reasoning that decides between them.

Two questions settle almost every case:

1. **Which path is actually carrying the noise?** Filtering the other one is
   pure lag.
2. **What does the lag land on?** Lag on the estimate is survivable. Lag on the
   stabilising gains is not.

## The two filter families in play

| | sliding-window taps | `mfc_iir_smoother` |
|---|---|---|
| type | FIR | IIR |
| poles | none | two, repeated real at $W/(W+1)$ |
| memory | exact: gone after $T_w$ | infinite tail |
| stability | cannot fail | a design question |
| cost | $n{+}1$ multiplies, $n{+}1$ stored | 3 multiplies, 2 stored |

Everything else in the discussion — Butterworth, Bessel, exponential
forgetting, an observer — is IIR. Poles mean recursion.

**On Butterworth specifically**, since it comes up: it is the *same family* as
our smoother, one parameter apart. Both are second-order all-pole low-pass with
unity DC gain; ours is critically damped ($\zeta = 1$, repeated real pole),
Butterworth is $\zeta = 0.707$ (complex pair), plus two zeros at Nyquist from
the bilinear transform. Matched at the same $-3$ dB point, Butterworth rolls off
harder ($-40.6$ vs $-32.3$ dB a decade out) and averages slightly more noise per
unit cutoff — and overshoots a step by 4.3 % where ours overshoots by exactly
zero, with 15 % more delay.

The reason this repo uses the critically damped one is **not** filter quality:

- a reference trajectory that overshoots is an instruction to the plant to
  overshoot, and `sp_filt` sets `err`;
- the impulse response is non-negative, so the output can never leave the range
  of its recent inputs — no ringing, ever;
- $W = 0$ is an exact pass-through, which is why `mfc_smoother_block` needs no
  enable flag;
- $W$ is *memory in samples*, the same unit the estimator windows use, rather
  than a cutoff in Hz that has to be re-derived whenever $T_s$ changes.

For the $\hat F$ post-filter role, where there is no setpoint and overshoot
costs nothing, a Butterworth would genuinely be the better choice on the
numbers. See [[iir-smoother]].

## Rule 1 — filter ORDER, not filter length

The single biggest lever, and the least obvious. Estimating $\ddot y$ amplifies
noise like $\omega^2$. Put a two-pole low-pass against that and the noise power
integral still has the **entire band above cutoff** contributing — two poles is
marginally insufficient for a second derivative. Add a second stage and it
finally collapses.

Consequence, measured on the algebraic estimator at matched delay (each 2-pole
stage costs $2W$ samples of lag):

| configuration | delay [samples] | $\hat F$ noise |
|---|---|---|
| internal $W = 160$ | 320 | 0.071 |
| internal $W = 10$ + post-filter $W = 40$ | 100 | 0.077 |
| internal $W = 640$ | 1280 | 0.020 |
| internal $W = 40$ + post-filter $W = 40$ | 160 | 0.024 |

**Same noise for 3–8x less lag**, purely from splitting the same filtering into
two stages. So: when a single window is not quiet enough, cascade an
`mfc_smoother_block` rather than lengthening the one you have.

Rule of thumb: **filter order must exceed the derivative order being
estimated.** For $\ddot y$ that means 4 poles — two of our smoothers.

## Rule 2 — never pre-filter $y$

Filtering the measurement *before* the estimator looks equivalent to filtering
$\hat F$ after it (for noise, it is — they commute). It is not, and the damage
does not need any noise to appear. On a **noiseless** run with a moving input
$u = \cos(2\pi \cdot 3t)$ and true $F = 3$:

| configuration | mean $\hat F$ | ripple |
|---|---|---|
| no extra filter | 2.985 | 0.018 |
| pre-filter $y$, $W = 40$ | 2.820 | **0.957** |
| post-filter $\hat F$, $W = 40$ | 2.980 | 0.011 |

Every estimator here computes something of the form $(\text{terms in } y) -
\alpha u$. Delay $y$ and you have delayed one side of that subtraction and not
the other, so the input cancellation no longer lines up. The algebraic variants
take a second injury: their terms are weighted by explicit powers of $t$ matched
to the actual $y(t)$, and a delayed $y$ no longer matches its own weights.

Post-filtering $\hat F$ filters both halves together, after the subtraction. Same
noise benefit, none of the damage. **If you need to remove something specific
from $y$** — mains hum, a rotor tone — use a notch, which is near-zero phase away
from its notch frequency and does not desynchronise $y$ from $u$.

## Rule 3 — never put the stabilising gains behind a filter

One principle, three faces:

- **Coupled estimators.** $K_p$ and $K_d$ live *inside* $\hat F$, so an $\hat F$
  post-filter delays your poles. Measured on identical plants and gains: a
  decoupled loop was untouched by a post-filter out to 320 samples of lag; the
  coupled loop showed 10 % overshoot at 40 samples and went **unstable at 160**.
  Coupled loops have almost no filtering headroom — take their noise reduction
  from `est_filter_window`, which sits inside the numerator/denominator symmetry
  and does not delay the folded gains relative to the estimate.
- **Filtering $u$.** Since $u = (-\hat F + ff - fb)/\alpha$, a filter there
  delays *all three* terms — including $fb$, the stabilising PD, and $ff$, which
  was not in the loop at all. It is the most expensive place to spend lag.
- **"Matching" phase between the two paths is not a goal.** $\hat F$ and $fb$ are
  not meant to cancel each other — one cancels the plant, the other places the
  poles. Adding lag to a healthy path so it matches a laggy one is strictly
  worse. Minimise lag in each path independently.

## Rule 4 — `u_prev` must be the command that actually reached the plant

This is the invariant in [[block-library-signal-flow]], and MFC is unusually
brutal about it. The window estimator's input taps sum to exactly $-1$, so the
input term contributes $\approx -\alpha u$; and the command law gives
$-\alpha u = \hat F - ff + fb$. **$\hat F$ therefore appears on the right-hand
side of its own estimate with coefficient $+1$** — an algebraic loop with unity
gain.

It is harmless only because the measurement path cancels it exactly: the plant
answers $\ddot y = F + \alpha u$, the $y$ taps return that, and the $+\alpha u$
and $-\alpha u$ annihilate, leaving $F$.

Tell the estimator about a $u$ the plant never received and that cancellation
breaks by exactly the discrepancy — leaving unity-gain positive feedback with
nothing against it. So if you filter $u$, the Unit Delay must tap the wire
**after** the filter. Measured: with $u_{prev}$ taken before the filter the loop
went unstable at half the filter strength that was still comfortable with it
taken after.

## Diagnose before you filter

Compare the two paths into $u$:

$$\text{D-term noise} \approx \frac{K_d\,\sigma}{T_s}
\qquad\text{vs}\qquad
\text{estimator noise} \approx \frac{\sigma_{\hat F}}{\alpha}$$

In a representative decoupled loop ($\sigma = 10^{-4}$, $T_s = 1$ ms,
$K_d = 16$) these came out **1.16 and 0.12** — the derivative term carried ten
times more noise than the estimator did. Post-filtering $\hat F$ in that loop
changed $\operatorname{std}(u)$ from 1.152 to 1.160: nothing, at the cost of real
lag.

Filtering $u$ is the only option that catches both paths, which is exactly why
it is tempting and why it costs the most. A short one is genuinely good value —
$W = 5$–$10$ bought a 16–29x reduction in command noise with no measurable
change in the step response — but it turns fast: $W = 20$ was already showing
overshoot and $W = 80$ was gone.

## The structural answer beats the filter

If the derivative term is the dominant noise path, the best move is not to
filter it but to **not have one**. That is [[riachy-trick]]: $K_d$ moves inside
$Y = y + K_d \int y$, and an *integral* of the measurement attenuates noise
instead of amplifying it. Same plant, same closed-loop poles, no filter and no
added lag:

| | $\operatorname{std}(u)$ | from $\hat F$ | from the D term |
|---|---|---|---|
| classic iPD | 1.152 | 0.117 | 1.157 |
| Riachy (no D) | **0.120** | 0.119 | — |

A 10x quieter command for free. The residue is then the $\hat F$ path, which
*is* worth post-filtering, because now it is the term that dominates.

This is arguably a stronger argument for Riachy's trick than the one it is
usually sold on ("an iPD without estimating derivatives").

## What lag costs even when it is safe

Estimator lag never destabilises a decoupled loop on a plant its PD could hold
alone, but it is not free:

- **Disturbance rejection degrades roughly linearly with lag**, always. $\hat F$
  *is* the disturbance rejection; a delayed estimate is delayed rejection.
- **Step response degrades in proportion to how much of the plant's own dynamics
  $\hat F$ is cancelling.** On an open-loop stable plant, nothing measurable out
  to 640 samples of lag. On an open-loop unstable one, overshoot went from
  1.0004 to 1.78 across the same sweep — there $\hat F$ is doing the
  stabilising, so its lag is loop lag.

## Decision order

1. **Diagnose.** $K_d\sigma/T_s$ vs $\sigma_{\hat F}/\alpha$.
2. **If the D term dominates**, remove it — [[riachy-trick]] — before reaching
   for any filter.
3. **If $\hat F$ dominates**, lengthen the estimator window first: the
   sliding-window taps are a better-shaped average than a bolted-on IIR, because
   they *are* the estimator ([[estimator-sliding-window]]).
4. **Then cascade** an `mfc_smoother_block` on the $\hat F$ wire. Never on $y$.
5. **A short filter on $u$** ($W \approx 5$–$10$) last, and only with `u_prev`
   tapped after it.
6. **On a coupled loop**, skip steps 3–5 and use `est_filter_window`.

Heavy post-filtering also extends the effective startup: the estimate has to
climb out of its initial zero through the filter, so raise `est_hold_time` to
match or expect a longer settle.

## See also

- [[riachy-trick]] — the derivation, and the D-term noise path it removes
- [[estimator-sliding-window]] — the FIR taps, and why they are the better average
- [[estimator-algebraic-2nd]] — the recursion whose noise this is mostly about
- [[iir-smoother]] — the two-pole filter itself, and its two internal roles
- [[control-law-coupled-vs-decoupled]] — why coupled has no filtering headroom
- [[block-library-signal-flow]] — the `u_prev` invariant, stated as wiring
