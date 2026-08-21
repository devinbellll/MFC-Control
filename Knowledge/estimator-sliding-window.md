# Sliding-window estimator (exact FIR taps)

`functions/mfc_fhat_sliding_window.m` + `mfc_siso.window_kernel` ·
block `mfc_fhat_window_block`

Recovers $F$ from a **fixed-length** window of the last $T_w$ seconds by direct
numerical quadrature, rather than by the operational calculus of
[[estimator-algebraic-2nd]].

## The integrals

**First order** ($\dot y = F + \alpha u$):

$$
F = -\frac{6}{T_w^3}\int_0^{T_w}\Big[(T_w - 2\sigma)\,y(\sigma) + \alpha\,\sigma(T_w - \sigma)\,u(\sigma)\Big]\,d\sigma
$$

**Second order** ($\ddot y = F + \alpha u$):

$$
F = \frac{60}{T_w^5}\int_0^{T_w}\Big[(T_w^2 - 6T_w\sigma + 6\sigma^2)\,y(\sigma) - \frac{\alpha}{2}\sigma^2(T_w - \sigma)^2\,u(\sigma)\Big]\,d\sigma
$$

The polynomial weights are orthogonality kernels: they are constructed so the
unknown initial conditions integrate to zero over the window, which is the same
job the $s$-differentiation does in the algebraic method — done by choice of
weight function instead.

## Implementation

Everything signal-independent is precomputed once, in `mfc_siso.window_kernel`,
as one fixed multiplier — a **tap** — per stored sample:

```matlab
n     = window_samples;          % no rounding; Tw = n*Ts exactly
sigma = (0:n).' * Ts;

if model_order == 1
    cy = [Tw, -2, 0, 0, 0];              % Tw - 2*sigma
    cu = [0, Tw, -1, 0, 0];              % sigma*(Tw - sigma)
    prefactor = -6 / Tw^3;
else
    cy = [Tw^2, -6*Tw, 6, 0, 0];         % Tw^2 - 6*Tw*s + 6*s^2
    cu = [0, 0, -0.5*Tw^2, Tw, -0.5];    % -0.5*s^2*(Tw - s)^2
    prefactor = 60 / Tw^5;
end
```

At run time the estimator is then a pure FIR filter — no division, no
recursion, no state beyond the two windows:

```matlab
integral = kernel.tap_y.' * state.y_buf + alpha * (kernel.tap_u_unit.' * state.u_buf);
F_hat    = integral;      % when valid
```

The prefactor, the sign and the $\tfrac12$ of the second-order input kernel are
folded into the taps; $\alpha$ is deliberately left **out** of `tap_u_unit` so it
can stay a live run-time input multiplying a constant array.

## How the taps are computed: exact, not quadrature

The weighting kernels above are polynomials *we wrote down ourselves*. There is
no reason to approximate them — the only genuine ignorance is what the signals
did **between** samples. So each tap is the exact integral of the kernel against
that sample's interpolation basis:

- **`y` is piecewise linear** (a tent basis — it is sampled, and nothing better
  is known):
  $$\text{tap}_y(i) = \int K_y(\sigma)\,\text{tent}_i(\sigma)\,d\sigma$$
  a cubic (1st order) or quartic (2nd order) integral over the two half-intervals
  either side of node $i$. Closed form, power rule, no quadrature rule involved.
- **`u` is piecewise constant**, which is not a modelling assumption at all but
  the truth: the command reaches the plant through a zero-order hold. So
  `tap_u_unit(j)` is the exact integral of $K_u$ over the one interval that
  sample was held across, and **the input term carries no error whatsoever**.

That also settles the alignment question: `u_buf(end)` is `u_prev`, the command
held over the interval *ending now*, so it owns the last interval, and the oldest
`u` sample — whose interval has fallen out of the window — gets a zero tap.

The whole thing is `mfc_siso.poly_moment` applied a few times:

$$\int_a^b \sigma^j p(\sigma)\,d\sigma = \sum_k c_{k}\,\frac{b^{k+j+1} - a^{k+j+1}}{k+j+1}$$

### What this replaced, and why

Composite Simpson was used here originally. It is a general-purpose rule: it
approximates the *whole product* $K \cdot y$ — including the $K$ we know
exactly — and it weights samples $1, 4, 2, 4, \dots, 1$. Measured against the
DC gain of the estimator (feed it $y = a\sigma^2/2$; it must return exactly $a$):

| $n$ | trapezoid | Simpson | exact taps |
|---|---|---|---|
| 10 | 1.1994 | 1.0024 | 0.99990 |
| 40 | 1.0125 | 1.0000094 | 0.9999996 |

Simpson fixed trapezoid's 20 % gain error, and exact taps improve on Simpson by a
further steady factor of ~24 at any window length. Two more things came free:

- the lumpy $4, 2, 4, 2$ weighting let a few samples dominate, which cost
  **~10 % more output noise** (sum of squared taps: 81231 vs 73723 at $n = 10$)
  at every window size;
- Simpson needed an **even** interval count, so `window_samples` was silently
  rounded up. It no longer is: any $n \ge 2$ is legal and means what it says.

A residual remains at second order and is *correct*: `y` is genuinely only known
at the samples, and a straight line between them is not a parabola.
`tests/octave_sanity.m` §2 pins it at `1e-5` (it was `1e-3` under Simpson), and
§3b pins the three moment conditions on the taps themselves:

```
sum(tap_y)            = 0     blind to a constant y   (no acceleration)
sum(sigma.*tap_y)     = 0     blind to a ramp y       (still none)
tap_y'*(sigma.^2/2)   = 1     unity gain on acceleration
sum(tap_u_unit)       = -1    exact, to machine precision
```

The first two are where the unknown initial position and velocity went: the two
$d/ds$ derivatives of the derivation reappear out here as two vanishing moments
of the tap list.

## Decoupled by construction

This estimator sees only $(y, u, \alpha)$ — never the tracking error. So $\hat F$
is always the true plant lump, there is nothing to fold, and there is no coupled
mode. `mfc_siso.config` **rejects** coupled + sliding window rather than
silently mis-running it.

It also has **no internal smoothing at all** (no numerator/denominator to filter
— there is no ratio). If the estimate is noisy, follow it with an
`mfc_smoother_block` on `F_hat`. That is a genuine post-filter, and here it is
the right tool, unlike in the algebraic case.

## Why you would choose it

| | algebraic (growing) | sliding window |
|---|---|---|
| memory | grows from $t=0$ | fixed $T_w$ |
| time origin | **critical** — weights are powers of $t$ | irrelevant; $t$ only gates the startup |
| long runs | weights grow unbounded, sluggish | unaffected |
| restarts | must reset the clock | safe |
| coupled folding | yes | no |
| internal smoothing | yes (num/den) | none |
| state | 6 scalars | two $(n{+}1)$ buffers |

**Use the sliding window for long or restarting runs**, where the growing window
degrades. Use the algebraic estimators when you want coupled folding or a
smaller state.

## The startup hold

```matlab
valid = t > kernel.Tw;      % hold until the window has filled
```

Simpler than the algebraic hold: there is no near-zero denominator, only a
buffer that has not filled with real data yet. Before $T_w$ the buffers are
still partly zeros and the integral is meaningless, so $\hat F$ is exactly 0.

## An implementation detail worth knowing

```matlab
state.y_buf      = circshift(state.y_buf, -1);
state.y_buf(end) = y;
```

`circshift` + end-assignment, not `[buf(2:end); new]`. Identical result, but it
also stays valid for the $1\times1$ placeholder buffers that algebraic variants
carry — whose dead copy of this code is still compiled under code generation.
See [[codegen-constraints]].

## Parameters

| Parameter | Meaning |
|---|---|
| `model_order` | 1 or 2; picks the kernel pair and prefactor, hence the taps |
| `Ts` | sample time; the buffers shift once per `Ts` |
| `window_samples` | window length in intervals; $T_w = n\,T_s$ exactly, no rounding |
| `alpha` | ultra-local model gain; must match `mfc_command_block` |

## See also

- [[estimator-algebraic-2nd]] — the growing-window alternative
- [[iir-smoother]] — what to put after `F_hat` if it is noisy
- [[codegen-constraints]] — the buffer and struct-shape constraints
