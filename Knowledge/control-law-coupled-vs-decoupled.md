# Coupled vs decoupled: where the stabilizing dynamics live

Both structures produce the same ideal closed loop. They differ in **who applies
the gains** — the estimator, or the feedback path.

| | **decoupled** | **coupled** |
|---|---|---|
| estimator is driven by | the measurement $y$ | the tracking error $e$ |
| $\hat F$ means | true plant lumped dynamics | plant dynamics *plus* the closed-loop polynomial |
| $K_p$ | explicit in `fb` | **folded** into $\hat F$ |
| $K_d$ | explicit in `fb` | **folded** at 2nd order; *unavailable* at 1st |
| $K_i$ | explicit | explicit (always) |
| works with sliding window | yes | **no** |

The split is decided by coupled/decoupled **alone**. Model order only chooses
whether the feedforward is $\dot y^*$ or $\ddot y^*$.

## Decoupled

The estimator sees only $(y, u, \alpha)$, so $\hat F$ is the honest plant lump.
Stability is then someone else's job — an explicit iPD(I) law:

$$
\text{fb} = K_d \dot e + K_p e + K_i \textstyle\int e
$$

```matlab
% mfc_siso.feedback, decoupled branch
fb = Kd*dot_err + Kp*err + Ki*int_err;
```

Substituting into the inverted model gives $e^{(n)} + K_d \dot e + K_p e = 0$
when $\hat F = F$.

**Use this by default.** $\hat F$ is physically meaningful and can be plotted,
compared against an analytic value, or reused. Every estimator supports it.

## Coupled

Drive the estimator with $e$ instead of $y$, and hand it the closed-loop
polynomial as *known* coefficients. The second-order estimator is told:

$$
\ddot e = F + \alpha u + a_{\text{fold}}\,\dot e + b_{\text{fold}}\, e,
\qquad a_{\text{fold}} = -K_d,\quad b_{\text{fold}} = -K_p
$$

```matlab
% mfc_siso.step
if cfg.coupled
    z_drive = err;                      % error-driven
    a_fold = -cfg.Kd;  b_fold = -cfg.Kp;
end
...
fb = Ki*int_err;                        % Kp and Kd are already inside F_hat
```

Because $s^2 + K_d s + K_p$ is now inside the estimate, the command law needs no
explicit P or D term at all — only the integral, which is never folded.

Why bother? The folded terms are annihilated by the *same* transform that
removes the initial conditions (see [[estimator-algebraic-2nd]]), so they cost
nothing numerically, and $K_d$ is applied without ever forming $\dot e$ — no raw
differentiation of a noisy error signal. That is a real noise advantage.

The price: $\hat F$ is no longer the plant's $F$. It is a hybrid quantity. Do not
read it as a disturbance estimate.

## The two illegal combinations

### Coupled + sliding window — rejected

The sliding-window estimator never sees the error; it is decoupled by
construction, so there is nothing to fold into. `mfc_siso.config` rejects the
combination outright rather than silently mis-running:

```matlab
assert(~(coupled && ~algebraic), ...
    ['mfc_siso.config: the coupled (error-driven, pole-folded) structure is ', ...
     'only defined for the algebraic estimator. ...']);
```

This is the missing corner of the $2\times2\times2$ grid — five supported
variants plus one that diverges, out of eight.

### Coupled + first order — legal, but undampable

The first-order model is $\dot z = F + \alpha u + b_{\text{fold}} z$. There is
**no $\dot z$ term** to fold a derivative gain into. So $K_p$ folds, and $K_d$
has nowhere to go — it is silently unused.

On a plant that needs derivative action this is fatal. In `examples/val_mfc.m`
the plant is a true double integrator (altitude from thrust); stabilizing $z$
from a first-order model fundamentally requires velocity feedback, and the
coupled first-order variant **cannot supply it at any gain**. It diverges, and
`tests/golden/1st_coupled_alg.csv` pins that divergence as expected behaviour.

If you need $K_d$ at first order, **run decoupled** — `mfc_feedback_block`
applies it explicitly there, at either order.

## Getting it wrong

The most common wiring mistake in an assembled loop: leaving `mfc_feedback_block`
on its default `coupled = false` while using a *coupled* estimator. Then $K_p$ is
applied twice — once folded, once explicit — and the loop runs at double the
proportional gain you think you set. Nothing errors; it just behaves oddly.

`mfc_feedback_block`'s `coupled` flag must match the estimator you wired up.

## See also

- [[ultra-local-model]] — the model both structures invert
- [[estimator-algebraic-2nd]] — how folding is actually implemented
- [[estimator-algebraic-1st]] — why there is no `a_fold` at first order
- [[block-library-signal-flow]] — matching the flags across blocks
