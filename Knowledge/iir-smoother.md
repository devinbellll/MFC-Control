# The IIR smoother

`functions/mfc_iir_smoother.m` · block `mfc_smoother_block`

One filter, used in four places. Eleven lines of code, and the only filter in
the whole controller.

## Definition

Unity-DC-gain, critically damped, second-order IIR low-pass:

$$
x_f[k] = \frac{x[k] + (2W^2 + 2W)\,x_f[k-1] - W^2\,x_f[k-2]}{W^2 + 2W + 1}
$$

```matlab
x_filt = (x_raw + (2*window^2 + 2*window)*x_filt_km1 - window^2*x_filt_km2) ...
         / (window^2 + 2*window + 1);
```

$W$ is **dimensionless** — memory length in samples, not seconds. It does not
need to be an integer.

## Properties

Factoring the denominator, the transfer function is

$$
H(z) = \frac{1}{(W+1)^2}\cdot\frac{1}{\left(1 - \frac{W}{W+1}z^{-1}\right)^{2}}
$$

- **repeated real pole** at $z = W/(W+1)$ → critically damped: no overshoot, no
  ringing. This matters for a reference trajectory, whose second derivative is
  fed straight through $1/\alpha$ to the actuator.
- **unity DC gain**: $H(1) = 1$. A constant passes through unchanged, so the
  filter never introduces steady-state error.
- **time constant** ≈ $W$ samples.
- $W = 0$ is an **exact** pass-through, not an approximate one — the expression
  collapses to $x_f = x$. This is used deliberately: `mfc_fhat_alg*_block`
  so "no filtering" needs no separate code path — and `mfc_smoother_block`
  therefore has no enable flag, only $W$.

All four properties are asserted in `tests/octave_sanity.m`.

## The four roles

| role | where | exposed as |
|---|---|---|
| 1. reference trajectory filter | `mfc_siso.ref_traj` | `mfc_smoother_block` |
| 2. estimator **numerator** filter | inside the algebraic estimators | — (internal) |
| 3. estimator **denominator** filter | inside the algebraic estimators | — (internal) |
| 4. $\hat F$ post-filter | optional, after any estimator | `mfc_smoother_block` |

One class in all four, so changing the smoother math is a single edit that
propagates everywhere. Role 1 also supplies $\dot y^*, \ddot y^*$ — enable the
derivative outputs.

**Roles 2 and 3 must use the same $W$**, which is exactly why they are *not*
blocks. They feed a division; identical filtering on both is what keeps the ratio
unbiased and largely cancels the filter's own lag, and a wiring you can get wrong
is a wiring that will be got wrong. The estimator blocks take one
`est_filter_window` and apply it to both. See [[estimator-algebraic-2nd]].

**Role 4 is a different operation** from roles 2–3, even though it is the same
filter. Post-filtering $\hat F$ lags the estimate; filtering num and den before
dividing does not, to first order. Prefer the internal path where it exists;
role 4 is for [[estimator-sliding-window]], which has no internal smoothing.

## Role 1: why a reference *trajectory*

The command law feeds $\ddot y^*$ straight through $1/\alpha$. Differentiate a
raw step twice and you get an impulse pair — an actuator command spike of
magnitude $\sim 1/T_s^2$.

So the setpoint is filtered into something twice-differentiable first, and the
derivatives are backward differences of the **filtered** history:

```matlab
sp_filt = mfc_iir_smoother(setpoint, sp_km1, sp_km2, window);
dot_sp  = (sp_filt - sp_km1) / Ts;
ddot_sp = (sp_filt - 2*sp_km1 + sp_km2) / Ts^2;
```

With $W = 0$ the smoother is an exact pass-through and the derivatives become
finite differences of the **raw** setpoint — impulsive for a step. That mode
exists for feeding an already-smooth externally generated trajectory, not for
step commands. `mfc_smoother_block` therefore has no separate enable flag:
bypassing it and setting $W = 0$ are the same operation.

$W$ here trades tracking aggressiveness against command smoothness. It is the
second knob to reach for, after $\alpha$.

## Choosing $W$

| | small $W$ | large $W$ |
|---|---|---|
| role 1 | snappy tracking, harsh command | gentle command, laggy tracking |
| roles 2–3 | responsive $\hat F$, noisy | smooth $\hat F$, slow to react |
| role 4 | little effect | visible lag in the loop |

Defaults are $W = 10$ for both `ref_filter_window` and `est_filter_window`, at
$T_s = 0.01$ — roughly a 0.1 s memory in both.

## See also

- [[estimator-algebraic-2nd]] — why num and den are filtered separately
- [[estimator-sliding-window]] — the estimator with no internal smoothing
- [[block-library-signal-flow]] — wiring the smoother into each role
