# Algebraic estimator, second order

`functions/mfc_fhat_algebraic_second_order.m` · blocks
`mfc_fhat_alg2_decoupled_block`, `mfc_fhat_alg2_coupled_block`

The default estimator. Recovers $F$ from a **growing** window starting at
$t = 0$, using operational calculus rather than quadrature.

## The model

$$
\ddot z = F + \alpha u + a_{\text{fold}}\,\dot z + b_{\text{fold}}\, z
$$

$z$ is the drive signal: the measurement (decoupled, $a=b=0$) or the tracking
error (coupled, $a = -K_d$, $b = -K_p$). See
[[control-law-coupled-vs-decoupled]].

## Derivation

**1. To the Laplace domain**, keeping the initial conditions explicitly — they
are unknown, and getting rid of them is the entire trick. With $F$ treated as
locally constant ($F/s$):

$$
s^2 Z(s) - s z(0) - \dot z(0) = \frac{F}{s} + \alpha U(s) + a\big(sZ(s) - z(0)\big) + b\,Z(s)
$$

**2. Annihilate the initial conditions.** $z(0)$ and $\dot z(0)$ appear
multiplied by $s^1$ and $s^0$. Differentiating twice with respect to $s$ kills
any term polynomial in $s$ of degree $< 2$ — both vanish at once, and no
knowledge of the plant's state is needed.

$$
\frac{d^2}{ds^2}\Big[\,\cdot\,\Big]
$$

**3. Back to the time domain.** Under the Laplace transform,

$$
\frac{d^n}{ds^n} \;\longleftrightarrow\; (-t)^n
$$

so differentiating in $s$ becomes **multiplying by powers of $t$** in the time
domain. Every $s$-derivative becomes a time-weighted integral of signals you
already have. No differentiation of measurements, and no plant model.

The result is a ratio $\hat F = \text{num}/\text{den}$ where $\text{den} = t^2$.

## Discretization

Backward differences on the time-weighted terms (`t` is the growing window
length, so these weights *grow* with the run):

```matlab
s_dz   = -(t*z - (t - Ts)*state.z_km1) / Ts;                                       % d/ds (sZ)
s2_d2z =  (t^2*z - 2*(t - Ts)^2*state.z_km1 + (t - 2*Ts)^2*state.z_km2) / Ts^2;
s_d2z  =  (t^2*z - (t - Ts)^2*state.z_km1) / Ts;
dz     = -t*z;
d2z    =  t^2*z;
d2u    =  t^2*u_prev;

num_raw = 2*z + 4*s_dz + s2_d2z - a_fold*(2*dz + s_d2z) - b_fold*d2z - alpha*d2u;
den_raw = t^2;
```

Note the folded terms enter `num_raw` on exactly the same footing as
$\alpha u$ — all are known-coefficient terms, all annihilated by the same $t^2$
transform. Folding is free.

## Smoothing, and why it happens *before* the divide

$$
\hat F = \frac{\mathcal{S}[\text{num}]}{\mathcal{S}[\text{den}]}
$$

```matlab
num_filt = mfc_iir_smoother(num_raw, state.num_filt_km1, state.num_filt_km2, filter_window);
den_filt = mfc_iir_smoother(den_raw, state.den_filt_km1, state.den_filt_km2, filter_window);
```

Numerator and denominator are smoothed **separately and identically**, then
divided. This is not a stylistic choice:

- smoothing the ratio instead would bias it — $\mathbb{E}[a/b] \neq \mathbb{E}[a]/\mathbb{E}[b]$
- **identical** filtering on both means the smoother's transfer function largely
  cancels in the ratio, so the estimate is smoothed without being lagged as much
  as a naive post-filter would lag it

This is why "put a filter on $\hat F$" is a *different operation* from what the
estimator does internally, and why the pipeline exposes both — see
[[block-library-signal-flow]].

## The startup hold

$\text{den} = t^2$ is zero at $t=0$ and tiny just after. Dividing by it during
the startup transient produces an enormous $\hat F$, which the controller
faithfully converts into an enormous command.

```matlab
valid = (den_filt ~= 0) && (t > hold_time);
if valid, F_hat = num_filt / den_filt; else, F_hat = 0; end
```

`est_hold_time` (default `0.1 s`) refuses to divide until the denominator has
grown, and $\hat F$ is exactly `0` until then. The functions still return a
`valid` flag in their debug struct, but no block exposes it — during the hold
$\hat F = 0$ says the same thing, and the flag was noise on the interface.

## The time origin matters

The window **grows from $t = 0$** and the weights are powers of $t$. Consequences
you have to design around:

- feed it a clock that **starts with the run**. A wall clock, or a clock that
  survives a restart, gives meaningless weights.
- as $t$ grows the window covers ever more history, so the estimator gets
  progressively more sluggish. On a long run this matters.
- $t^2$ grows without bound — at very long runtimes, numerator and denominator
  both become large and precision suffers.

If any of that bites, use [[estimator-sliding-window]], which has finite memory
and no time origin at all.

## Why the num/den filtering is not exposed

The smoothing sits **before** the division and both sides must use the same $W$;
that is what keeps the ratio unbiased. Splitting it into wireable blocks would
make an essential invariant into something you can mis-wire, for no gain — so the
estimator blocks take one `est_filter_window` and apply it to both sides
internally. The raw numerator and denominator are not output at all.

If you want a *smoothed $\hat F$*, that is a genuinely different operation:
put an `mfc_smoother_block` after the estimator. It lags the estimate, which
filtering num and den does not, to first order. See [[iir-smoother]].

## Parameters

| Parameter | Meaning |
|---|---|
| `Ts` | sample time; the backward differences assume one advance per `Ts` |
| `est_filter_window` | $W$ of the num/den smoother, in samples (both sides) |
| `est_hold_time` | hold $\hat F$ at 0 while $t \le$ this |
| `Kp`, `Kd` | **coupled block only**; folded in as $b_{\text{fold}} = -K_p$, $a_{\text{fold}} = -K_d$. The decoupled block has neither. |
| `alpha` | ultra-local model gain; must match `mfc_command_block` |

## See also

- [[estimator-algebraic-1st]] — same method, one derivative, no `a_fold`
- [[estimator-sliding-window]] — finite memory, time-origin independent
- [[iir-smoother]] — the filter used on num and den
- [[control-law-coupled-vs-decoupled]] — what to fold and when
- [[estimator-algebraic-2nd-stock-blocks]] — the same estimator drawn with core
  Simulink blocks, discrete and continuous
- [[estimator-free-derivatives]] — what the $n=0$ and $n=1$ equations give you
  once you stop throwing the initial conditions away
