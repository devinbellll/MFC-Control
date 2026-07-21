# Algebraic estimator, first order

`functions/mfc_fhat_algebraic_first_order.m` · block `mfc_fhat_alg1_block`

The same operational-calculus method as [[estimator-algebraic-2nd]], applied to
the first-order ultra-local model. One derivative, one initial condition, one
$s$-differentiation.

## The model

$$
\dot z = F + \alpha u + b_{\text{fold}}\, z
$$

There is **no $a_{\text{fold}}$**. That absence is the most important fact on
this page — see below.

## Derivation

In the Laplace domain, with $F$ locally constant and $z(0)$ unknown:

$$
sZ(s) - z(0) = \frac{F}{s} + \alpha U(s) + b\,Z(s)
$$

$z(0)$ appears as a constant in $s$ — degree $0$. Differentiating **once** with
respect to $s$ annihilates it:

$$
\frac{d}{ds}\Big[\,\cdot\,\Big] \quad\Longrightarrow\quad z(0) \text{ vanishes}
$$

Mapping back with $\frac{d}{ds} \leftrightarrow -t$ gives a ratio with
$\text{den} = t$ — one power of $t$, where the second-order case had $t^2$.

$b_{\text{fold}}\,z$ is an order-0, known-coefficient term, exactly like
$\alpha u$, so the same $t^1$ transform annihilates it. Folding costs nothing.

## Discretization

```matlab
num_raw = -z + (t*z - (t - Ts)*state.z_km1)/Ts - t*alpha*u_prev - t*b_fold*z;
den_raw = t;
```

Then, identically to second order, numerator and denominator are smoothed
**separately with the same filter** before dividing, and the result is held at
zero until `t > hold_time`. The reasoning for both is in
[[estimator-algebraic-2nd]] and applies unchanged.

## No derivative room — the structural limitation

The first-order model has no $\dot z$ term. There is therefore nowhere to fold a
derivative gain, and `mfc_fhat_alg1_block` correctly has no `a_fold` parameter.

**Coupled at first order:** $K_p$ folds via `b_fold = -Kp`. $K_d$ is silently
**unused** — `mfc_siso.feedback` applies only $K_i$ in the coupled branch, and
the estimator has no slot for $K_d$. So a coupled first-order loop supplies *no
damping whatsoever*, at any gain.

On a plant that needs derivative action this is fatal, not merely suboptimal. In
`examples/val_mfc.m` the plant is a true double integrator; stabilizing position
from a first-order model requires velocity feedback, and the coupled first-order
variant cannot provide it. **It diverges**, and
`tests/golden/1st_coupled_alg.csv` pins that divergence as expected behaviour so
the failure mode stays honest rather than quietly changing.

**Decoupled at first order** has no such restriction: `mfc_feedback_block`
applies $K_d$ explicitly at either order. `val_mfc.m` tunes $K_{d,1} = p$ there
precisely to supply the damping the coupled variant cannot.

> If you want first order **and** you need $K_d$: run decoupled. There is no
> tuning that rescues coupled first order on such a plant.

## First vs second order, side by side

| | first order | second order |
|---|---|---|
| $s$-differentiations | 1 | 2 |
| initial conditions killed | $z(0)$ | $z(0)$, $\dot z(0)$ |
| denominator | $t$ | $t^2$ |
| foldable gains | $K_p$ only | $K_p$ and $K_d$ |
| feedforward | $\dot y^*$ | $\ddot y^*$ |
| state used | `z_km1` | `z_km1`, `z_km2` |

The shared `state` struct carries `z_km2` for both, so the two estimators are
interchangeable at the call site; the first-order one simply never reads it.

## Parameters

| Parameter | Meaning |
|---|---|
| `Ts` | sample time; the backward difference assumes one advance per `Ts` |
| `est_filter_window` | $W$ of the num/den smoother, in samples |
| `est_hold_time` | hold $\hat F$ at 0 while $t \le$ this |
| `b_fold` | $-K_p$ when coupled, $0$ when decoupled |
| `alpha` | ultra-local model gain; must match `mfc_command_block` |

## See also

- [[estimator-algebraic-2nd]] — the full derivation and the smoothing rationale
- [[control-law-coupled-vs-decoupled]] — why the missing $K_d$ matters
- [[ultra-local-model]] — choosing the order in the first place
