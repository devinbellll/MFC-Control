# Algebraic estimator, second order — built from stock blocks

`functions/mfc_fhat_algebraic_second_order.m` · block
`mfc_fhat_alg2_decoupled_block`

What [[estimator-algebraic-2nd]] computes, drawn entirely with core Simulink
blocks — no `matlab.System`, no MATLAB Function block. This is a **reading
aid and a porting reference**, not a replacement: the shipped block stays the
implementation, because the invariants at the bottom of this note are trivial
to break by hand.

Decoupled only ($a_{\text{fold}} = b_{\text{fold}} = 0$, drive signal $z = y$).
The coupled folding adds two more time-weighted terms on the same footing;
see [[control-law-coupled-vs-decoupled]].

## The restructuring that makes it drawable

The three operational-calculus terms look like three unrelated expressions, but
they are plain backward differences of **two** time-weighted products. Write

$$
p[k] = t_k\, z[k], \qquad q[k] = t_k^2\, z[k]
$$

Because $t_{k-1} = t - T_s$ and $t_{k-2} = t - 2T_s$, the shifted samples the
code writes out longhand *are* the delayed products:    

$$
p[k-1] = (t - T_s)\,z[k-1], \qquad q[k-2] = (t - 2T_s)^2\,z[k-2]
$$

so the discretization in [[estimator-algebraic-2nd]] collapses to

$$
\texttt{s\_dz} = -\frac{\Delta p}{T_s}, \qquad
\texttt{s2\_d2z} = \frac{\Delta^2 q}{T_s^2}, \qquad
\Delta = 1 - z^{-1}
$$

and the whole estimator is

$$
\text{num} = 2z - \frac{4}{T_s}\Delta\!\left[t z\right]
             + \frac{1}{T_s^2}\Delta^2\!\left[t^2 z\right]
             - \alpha\, t^2 u_{\text{prev}},
\qquad \text{den} = t^2
$$

$$
\hat F = \frac{\mathcal{S}[\text{num}]}{\mathcal{S}[\text{den}]}
\quad\text{(hold at 0 while } t \le \texttt{hold\_time)}
$$

**No history of $z$ is carried by hand.** `z_km1`/`z_km2` disappear into two
Discrete Filter blocks — the state moves from the signal to the differencer.

## Discrete

```
                        ┌──────────────┐
  t ──┬─────────────────►    t·z       │      ┌────────────┐
      │              ┌──►  (Product)   ├─────►│ [1 −1]     │──►(×) −4/Ts ─┐
      │              │  └──────────────┘      │ [1  0]     │              │
      │              │                        └────────────┘              │
      │  ┌────┐      │                          Δ, Ts                     │
      ├──►(u²)├──┬───┼──► t²                                              │
      │  └────┘  │   │   ┌──────────────┐    ┌────────────┐               ▼
      │          │   │   │   t²·z       │    │ [1 −2  1]  │            ┌──────┐
  z ──┴──────────┼───┴──►│  (Product)   ├───►│ [1  0  0]  │──(×)1/Ts² ─►      │
  (=y)           │       └──────────────┘    └────────────┘            │ Σ    │
      └──────────┼──────────────────► (×) 2 ───────────────────────────►      │
                 │                                                     │ ++−− │
  u_prev ────────┴────► t²·u_prev (Product) ──► (×) α ─────────────────►└──┬───┘
                                                                          │ num
                 t² ──────────────────────────────────► den               │
                          │                                               │
                    ┌─────▼───────────────┐                    ┌──────────▼──────────┐
                    │ Discrete Filter  S  │                    │ Discrete Filter  S  │
                    │ num [1]             │                    │ num [1]             │
                    │ den [a0 −a1 a2]     │                    │ den [a0 −a1 a2]     │
                    └─────┬───────────────┘                    └──────────┬──────────┘
                          │ den_filt                                      │ num_filt
                          └────────────────►┌─────────┐◄──────────────────┘
                                            │ Divide  │
                                            └────┬────┘
                                                 │
   t ──►[t > hold_time]──┐                       │
   den_filt ──►[≠ 0]─────┴──►(AND)──► ┌──────────▼─────────┐
                                      │ Switch  (else 0)   │──► F_hat
                                      └────────────────────┘
```

| Stock block | Setting |
|---|---|
| Clock | $t$, **starting with the run** — the growing window depends on the origin |
| Product / Math `u^2` | $t z$, $t^2$, $t^2 z$, $t^2 u_{\text{prev}}$ |
| Discrete Filter ($\Delta$) | num `[1 -1]`, den `[1 0]`, sample time `Ts`; then gain $-4/T_s$ |
| Discrete Filter ($\Delta^2$) | num `[1 -2 1]`, den `[1 0 0]`, sample time `Ts`; then gain $1/T_s^2$ |
| Discrete Filter ($\mathcal S$) ×2 | num `[1]`, den `[a0 −a1 a2]`, sample time `Ts` |
| Gain | $\alpha$ on the $t^2 u_{\text{prev}}$ branch |
| Sum | `+ + + −` |
| Divide | `num_filt / den_filt` |
| Compare To Constant ×2 + AND + Switch | $t > \texttt{hold\_time}$, $\texttt{den\_filt} \neq 0$, else `0` |

The smoother coefficients come straight from [[iir-smoother]]:

$$
a_0 = W^2 + 2W + 1, \qquad a_1 = 2W^2 + 2W, \qquad a_2 = W^2
$$

Every stateful block runs at `Ts` and **must not inherit** — the difference
operators assume exactly one advance per sample.

Verified against `mfc_fhat_algebraic_second_order` over a 400-step run
(random-perturbed sinusoidal $z$, $\alpha = 0.7$, $W = 10$): max deviation
$3.3\times10^{-12}$, i.e. float re-association only.

## Continuous — pure integrators, no differentiators

The continuous form is **not** the discrete graph with $\Delta/T_s$ replaced by
`du/dt`. Written that way it would need a second derivative of a measured
signal, which is exactly what the algebraic method exists to avoid. The honest
continuous estimator is built from **Integrator blocks only**, and it is the
form to port to.

### Derivation

Start from the ultra-local model in the Laplace domain, initial conditions kept
explicit, $F$ locally constant:

$$
s^2 Y - s\,y(0) - \dot y(0) = \frac{F}{s} + \alpha U
$$

Differentiate twice w.r.t. $s$. Both initial conditions are polynomial in $s$ of
degree $< 2$, so both vanish:

$$
2Y + 4sY' + s^2 Y'' - \alpha U'' = \frac{2F}{s^3}
$$

This still contains positive powers of $s$ — that is where differentiators would
come from. Remove them by multiplying through by $s^{-2}$:

$$
2\,s^{-2}Y + 4\,s^{-1}Y' + Y'' - \alpha\,s^{-2}U'' = \frac{2F}{s^5}
$$

Every operator is now an *integration*. Map back with $Y^{(n)} \leftrightarrow
(-t)^n y$ and $s^{-n} \leftrightarrow$ $n$-fold integration:

$$
\boxed{\;
\underbrace{2\!\iint\! y \;-\; 4\!\int_0^t \!\tau\,y(\tau)\,d\tau \;+\; t^2 y
\;-\; \alpha \!\iint\! \tau^2 u}_{\textstyle \text{num}(t)}
\;=\; F\,\underbrace{\frac{t^4}{12}}_{\textstyle \text{den}(t)} \;}
$$

$\hat F = \text{num}/\text{den}$. No derivative of $y$ or $u$ appears anywhere —
only products with powers of $t$, and integrals.

### The diagram

```
                                ┌─────┐   ┌─────┐
  y ──┬──────────────► ×2 ─────►│ 1/s ├──►│ 1/s ├────────────────────►┐
      │                         └─────┘   └─────┘                     │
      │        ┌──────┐         ┌─────┐                               │
      ├──►(×)─►│ t·y  ├── ×(−4)►│ 1/s ├──────────────────────────────►│
      │   ▲    └──────┘         └─────┘                             ┌─┴──┐
      │   │                                                         │ Σ  │ num
      └──►(×)──────────────────────────────────────► t²·y ─────────►│    ├──┐
          ▲ t²                                                      └─┬──┘  │
          │              ┌─────┐   ┌─────┐                            │     │
  u ─────►(×)─► t²·u ───►│ 1/s ├──►│ 1/s ├──► ×(−α) ──────────────────►┘    │
                         └─────┘   └─────┘                                  │
                                                                       ┌────▼───┐
  t ──┬──►(u²)──┬──► t²  (to the products above)                        │   ÷    │
      │         │                                                      └────▲───┘
      └─────────┴──►(u²)─► t⁴ ──► ×(1/12) ──► den ───────────────────────────┘
                                                     │
                                                     └──► F̂, gated to 0 while t ≤ hold_time
```

| Stock block | Setting |
|---|---|
| Clock | $t$, starting with the run |
| Integrator ×5 | initial condition `0` — two on $2y$, one on $t y$, two on $t^2 u$ |
| Product / Math `u^2` | $t^2$, $t^4$, $t y$, $t^2 y$, $t^2 u$ |
| Gain | $2$, $-4$, $-\alpha$, $1/12$ |
| Sum | `+ + + +` (the signs are already in the gains) |
| Divide, Compare To Constant, Switch | as in the discrete version |

Five integrators, no `du/dt` block, no Discrete Filter. Note the denominator is
now $t^4/12$, **not** $t^2$ — the extra $s^{-2}$ that bought the integral form
shows up as two more powers of $t$.

### Verified

Numerically, against a plant with known $F$ and **nonzero** initial conditions
($y(0) = 1.7$, $\dot y(0) = -0.9$, $\alpha = 0.7$, $u = \cos 2t$, trapezoid
quadrature at $dt = 10^{-4}$):

| $t$ | $\hat F$ | true |
|---|---|---|
| 1.0 | 2.99999993 | 3 |
| 2.0 | 2.99999998 | 3 |
| 6.0 | 3.00000000 | 3 |

Residual is quadrature error, shrinking with $t$. The nonzero ICs are the point:
they never enter, which is the annihilation working.

### The fully-integral variant

The bare $t^2 y$ term above carries no integrator, so it passes measurement
noise straight through. Multiplying by $s^{-3}$ instead of $s^{-2}$ integrates
that term too:

$$
2\!\iiint\! y \;-\; 4\!\iint\! \tau y \;+\; \int\! \tau^2 y
\;-\; \alpha\!\iiint\! \tau^2 u \;=\; F\,\frac{t^5}{60}
$$

Seven integrators, every path filtered, denominator $t^5/60$. Both forms recover
$F = 3$ exactly in the test above; prefer this one if $y$ is noisy, the shorter
one otherwise.

### What the continuous form does *not* fix

- **The growing window still dominates.** Both forms weight the whole history
  from $t=0$, so a *time-varying* $F$ is tracked with increasing lag: against
  $F(t) = 3 + 1.5\sin(0.5t)$ the estimate was already off by $\sim1$ at
  $t = 6$ in both variants. Continuous integration removes the differentiators,
  not the sluggishness — that needs [[estimator-sliding-window]].
- **The startup hold is still required.** $\text{den} \to 0$ as $t \to 0$
  regardless of form.
- **$u$ vs $u_{\text{prev}}$.** The discrete version needs the explicit Unit
  Delay — the command actually applied over the *last* sample. In continuous
  time that delay is not a state: feed $u$ directly, or a Transport Delay of
  $T_s$ to reproduce the discrete behaviour.

## What must not be "simplified"

In the **discrete** version, both of these are load-bearing and one careless
edit away in a hand-wired copy — which is the argument for keeping the shipped
block:

- $\mathcal{S}$ is the **same** filter on numerator and denominator. Identical
  filtering is what keeps the ratio unbiased and largely un-lagged; two
  instances with different windows silently ruins the estimate. See
  [[estimator-algebraic-2nd]].
- `den` is the **smoothed** $t^2$, not raw $t^2$. Cancelling it to
  $\text{num\_filt}/t^2$ changes the estimate — the cancellation only holds in
  the DC limit.

The continuous form needs neither: the integrators *are* the smoothing, and its
$\text{den} = t^4/12$ is exact. If you add filtering anyway, the same rule
applies — one filter, both sides.

## See also

- [[estimator-algebraic-2nd]] — the derivation and the discretization
- [[iir-smoother]] — where $a_0, a_1, a_2$ come from
- [[block-library-signal-flow]] — the loop this estimator sits in
- [[codegen-constraints]] — why the shipped version is a `matlab.System` class
