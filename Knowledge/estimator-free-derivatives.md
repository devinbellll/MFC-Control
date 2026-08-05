# Free derivative and initial-condition estimates

`functions/mfc_fhat_algebraic_second_order.m` · models `library/mfc_F_est_*.mdl`

The algebraic estimator throws away more information than it uses. Recovering
it costs a handful of blocks and no new measurements: $\ddot y$ comes out
**exactly free**, and $y(0)$, $\dot y(0)$ come from two equations the derivation
discards. Current-time $\dot y$ and $y$ follow by integrating the model
forward — cheap, but with a drift caveat that differs between the coupled and
decoupled structures.

Nothing here is implemented in `blocks/`. This note is the design, the formulas
and the numbers; see [[estimator-algebraic-2nd-stock-blocks]] for the graph the
extra blocks attach to.

## Why anything is free

The ultra-local model in the Laplace domain has **three** unknowns:

$$
s^2 Y - s\,y(0) - \dot y(0) = \frac{F}{s} + \alpha U
$$

$F$, $y(0)$, $\dot y(0)$. Differentiating $n$ times w.r.t. $s$ gives one
equation per $n$. [[estimator-algebraic-2nd]] uses $n = 2$ *because* it
annihilates $y(0)$ and $\dot y(0)$ — but $n = 1$ and $n = 0$ are equally valid
equations, and they contain the initial conditions rather than destroying them.

Three equations, three unknowns, and they are **triangular**: $n=2$ gives $F$
with no ICs in it, $n=1$ then gives $y(0)$, $n=0$ then gives $\dot y(0)$. You
already solve the hard one.

## The two extra equations

Multiplied by $s^{-2}$ so every operator is an integration — the same
convention as the shipped continuous model, which is what makes the terms
reusable:

$$
n=2:\quad 2\!\iint\! y - 4\!\int\! \tau y + t^2 y - \alpha\!\iint\! \tau^2 u
      \;=\; F\,\frac{t^4}{12}
$$

$$
n=1:\quad y(0) \;=\; \frac{1}{t}\left( 2\!\int\! y \;-\; t\,y
      \;+\; \hat F\,\frac{t^3}{6} \;+\; \alpha\!\iint\! \tau u \right)
$$

$$
n=0:\quad \dot y(0) \;=\; \frac{1}{t}\left( y \;-\; \hat y(0)
      \;-\; \hat F\,\frac{t^2}{2} \;-\; \alpha\!\iint\! u \right)
$$

Solve in that order. Both divide by $t$ (not $t^2$), so they need the same
startup hold as $\hat F$ and no more.

## What each one costs

Against the block graph of `library/mfc_F_est_cont.mdl`:

| Quantity | New blocks | Reused |
|---|---|---|
| $\ddot y$ | 1 Product ($\alpha u$), 1 Sum | $\hat F$, $u$, $\alpha$ |
| $y(0)$ | 1 Product ($tu$), 2 Integrators, 1 Sum, 1 Divide | $2\!\int\! y$, $t y$, $\hat F$, $t$ |
| $\dot y(0)$ | 2 Integrators, 1 Sum, 1 Divide | $y$, $\hat y(0)$, $\hat F$, $t$ |
| $\dot y(t)$, $y(t)$ | 2 Integrators, 2 Sums | everything above |

Two reuses are worth pointing out, because they are why this is cheap:

- **$2\!\int\! y$ is already a wire.** The $2\!\iint\! y$ term is built as
  `Gain(2) → Integrator → Integrator`; the intermediate signal is exactly the
  $2\!\int\! y$ that the $n=1$ equation wants. Tap it, do not rebuild it.
- **$t\,y$ is already a wire** — `Product_tz`, built for the $-4\!\int\!\tau y$
  term.

## $\ddot y$ is free and exact

For a **decoupled** estimator, $\hat F$ *is* the plant lumped dynamics
$\ddot y - \alpha u$ by definition. So

$$
\ddot{\hat y} = \hat F + \alpha u
$$

is not an approximation, it is the ultra-local model rearranged — one Sum and
one Product, no state, no lag beyond whatever lag $\hat F$ already has. If you
are logging $\hat F$ you are one block away from logging $\ddot y$.

## Current-time $\dot y$ and $y$

The ICs are the state at $t = 0$, not now. Integrate the model forward from
them:

$$
\dot{\hat y}(t) = \dot{\hat y}(0) + \int_0^t \left(\hat F + \alpha u\right), \qquad
\hat y(t) = \hat y(0) + \int_0^t \dot{\hat y}
$$

Verified against a plant with known $F=3$ and nonzero ICs
($y(0)=1.7$, $\dot y(0)=-0.9$, $\alpha=0.7$, $u=\cos 2t$, $dt=10^{-4}$):

| quantity | recovered | true | current-time error, $t>0.5$ |
|---|---|---|---|
| $\hat F$ | 3.000000 | 3 | — |
| $\hat y(0)$ | 1.700000 | 1.7 | — |
| $\dot{\hat y}(0)$ | −0.900000 | −0.9 | — |
| $\ddot{\hat y}$ | | | $1.8\times10^{-7}$ |
| $\dot{\hat y}$ | | | $2.5\times10^{-6}$ |
| $\hat y$ | | | $1.9\times10^{-5}$ |

The ICs are exact from $t \approx 0.2$ onward — the same point $\hat F$ becomes
usable.

## The drift caveat, and why it is structural

Forward integration is open loop, so any bias in $\hat F$ accumulates. **How
badly depends on which estimator you used**, and the difference is large enough
to drive the choice:

- **Decoupled.** The error obeys $\ddot e = \Delta F$ — a double integrator.
  Unbounded: a 1 % bias on $F=3$ gives $|e| = \Delta F\,t^2/2 = 6.0$ at
  $t = 20$ s (measured: 6.0000).
- **Coupled.** $\hat F$ has the loop poles folded in, so forward integration
  runs $\ddot e + K_d \dot e + K_p e = \Delta F$ — a *stable* second-order
  system. The error is **bounded** at $\Delta F / K_p$: the same 1 % bias gives
  0.0075 at $t=20$ with $K_p=4$, $K_d=2.5$ (measured: 0.0075, i.e. exactly
  $\Delta F/K_p$).

So the coupled structure gives drift-free integrated states, at the cost that
what you integrate is the **tracking error** $z = \text{err}$ and its
derivative, not $y$ and $\dot y$. Recover the measurement-frame values with the
setpoint derivatives the smoother already outputs, per
[[block-library-signal-flow]]:
$\dot y = \dot z + \dot{sp}_{\text{filt}}$.

Note also that with a coupled estimator $\ddot z \neq \hat F + \alpha u$ — the
folded terms mean $\ddot z = \hat F + \alpha u - K_d \dot z - K_p z$, which
needs $\dot z$ and so is not free. The "$\ddot y$ for one Sum" trick is
**decoupled only**.

## What this does not fix

Everything downstream of the growing window still applies. With a
time-varying $F$ the estimate lags, and the recovered derivatives inherit that
lag and then integrate it: against $F(t) = 3 + 1.5\sin(0.5t)$ the same setup
gave $\ddot y$ error 2.2 and $y$ error 5.9. These are free *derivatives of the
model*, not an independent observer — they are exactly as good as $\hat F$ is,
and no better. If that is the binding constraint, the window is the thing to
change ([[estimator-sliding-window]]), not this.

## See also

- [[estimator-algebraic-2nd]] — the $n=2$ equation and where the ICs went
- [[estimator-algebraic-2nd-stock-blocks]] — the block graph these attach to
- [[ultra-local-model]] — why $\hat F + \alpha u$ is the whole model
- [[control-law-coupled-vs-decoupled]] — the structure that decides the drift
