# The ultra-local model

The single assumption behind Model-Free Control: over a short time window, any
SISO plant behaves like a pure integrator chain driven by one known gain plus
one unknown lump.

$$
y^{(n)} = F + \alpha u, \qquad n \in \{1, 2\}
$$

- $\alpha$ is a **design parameter**, chosen, not identified.
- $F$ is **everything else** — nonlinearities, unmodelled dynamics,
  disturbances, parameter drift, and the part of the true input gain that
  $\alpha$ got wrong. It is re-estimated at every sample.

This repo implements both orders:

$$
\dot{y} = F + \alpha u \qquad\text{(first order)}
$$
$$
\ddot{y} = F + \alpha u \qquad\text{(second order)}
$$

## Why this is not cheating

$F$ is not assumed small, constant, or structured. It is assumed **estimable
over a short window**, which is much weaker. If you can measure $y$ and you know
what $u$ you applied, you can solve for $F$ over the last few samples — the
estimators in [[estimator-algebraic-2nd]] and [[estimator-sliding-window]] are
two ways of doing exactly that.

Control is then trivial: invert the model.

$$
u = \frac{-\hat{F} + y^{*(n)} - \text{fb}}{\alpha}
$$

If $\hat F \approx F$, substituting back gives $y^{(n)} = y^{*(n)} - \text{fb}$ —
the plant has been *reduced to a pure integrator chain*, and `fb` is free to be
a simple linear law on the error. All the plant's difficulty is absorbed into
$\hat F$ and cancelled. That is the whole method.

In code, one line — `mfc_siso.command`:

```matlab
u_raw = (-F_hat + ff - fb) / alpha;
```

## Choosing $\alpha$

$\alpha$ sets the loop's authority, not its correctness. The estimator absorbs
whatever mismatch is left, so the loop works over a surprisingly wide range.

- **too small** → $1/\alpha$ is large → over-driven, twitchy command
- **too large** → under-driven, sluggish
- a decent start is the true input gain's order of magnitude, if you know it

$\alpha$ is the first knob to try when a loop misbehaves. Both the estimator
blocks and `mfc_command_block` take a live `alpha` input for sweeping or
scheduling it — but **give both the same value**, or $\hat F$ and the inversion
disagree about which model they are talking about.

## Choosing the order

Not "what order is my plant" — it is **how many integrators the command has to
push through**.

| | first order | second order |
|---|---|---|
| model | $\dot y = F + \alpha u$ | $\ddot y = F + \alpha u$ |
| feedforward | $\dot y^*$ | $\ddot y^*$ |
| ideal closed loop | $\dfrac{K_p}{s + K_p}$ | $\dfrac{K_p}{s^2 + K_d s + K_p}$ |
| poles | one, at $-K_p$ | $s^2 + K_d s + K_p$; double pole at $-p$ ⟺ $K_d = 2p,\ K_p = p^2$ |

A second-order plant *can* be run first order — the extra pole just becomes part
of $F$. That works when the estimator is fast relative to the ignored dynamics
and fails when it is not. `examples/val_mfc.m` shows both on a plant with two
unmodelled poles (a velocity integrator and an actuator lag) so the degradation
is visible.

**The one hard restriction:** a first-order *coupled* loop has nowhere to fold a
derivative gain, so it cannot supply damping at any gain. On a true double
integrator it diverges — see [[control-law-coupled-vs-decoupled]]. This is a
structural fact, not a tuning problem.

## What $F$ ends up absorbing

In the `val_mfc` drone plant, $\ddot z = -g - d\dot z + \alpha_p u_{\text{act}}$
with a first-order actuator lag, the controller is told only
$\ddot z = F + \alpha u$. So $F$ absorbs:

| absorbed into $F$ | how well |
|---|---|
| constant gravity $-g$ | perfectly — a constant is the easy case |
| linear drag $-d\dot z$ | well — algebraic in the state |
| actuator lag $\tau$ | **not at all** — it is phase, not a state function |

That last row is the real design constraint. An estimator cannot estimate away a
delay; it can only chase it. So the closed-loop bandwidth must stay below the
actuator's — `val_mfc.m` picks $p = \frac{1}{5\tau}$, one fifth of the actuator
bandwidth, for exactly this reason.

## See also

- [[control-law-coupled-vs-decoupled]] — where the stabilizing dynamics live
- [[estimator-algebraic-2nd]] — the default way of getting $\hat F$
- [[estimator-sliding-window]] — the finite-memory alternative
- [[block-library-signal-flow]] — how the stages wire together
