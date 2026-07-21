# Sliding-window estimator (Simpson quadrature)

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

Everything signal-independent is precomputed once, in `mfc_siso.window_kernel`:

```matlab
n     = window_samples + mod(window_samples, 2);   % force even for Simpson
Tw    = n * Ts;
sigma = (0:n).' * Ts;

if model_order == 1
    y_kernel      = Tw - 2*sigma;
    u_kernel_unit = sigma .* (Tw - sigma);
    prefactor     = -6 / Tw^3;
else
    y_kernel      = Tw^2 - 6*Tw*sigma + 6*sigma.^2;
    u_kernel_unit = -0.5 * sigma.^2 .* (Tw - sigma).^2;
    prefactor     = 60 / Tw^5;
end
```

The sign and the $\tfrac12$ of the second-order input kernel are **folded into
`u_kernel_unit`**, so $\alpha$ can stay a live run-time input multiplying a
constant array. Both orders then evaluate through one common expression:

```matlab
integrand = kernel.y_kernel .* state.y_buf + alpha * kernel.u_kernel_unit .* state.u_buf;
integral  = (kernel.Ts/3) * sum(kernel.simpson_weights .* integrand);
F_hat     = kernel.prefactor * integral;      % when valid
```

## Simpson is mandatory, not a refinement

Composite Simpson, with an **even** interval count (`window_samples` is rounded
up, so the realized window is $T_w = n\,T_s$):

```matlab
simpson_weights            = ones(n + 1, 1);
simpson_weights(2:2:end-1) = 4;
simpson_weights(3:2:end-1) = 2;
```

Trapezoidal integration **cannot** be used at second order. Its leakage is
$O(T_s^2/T_w^4)$, which sounds negligible — but the $60/T_w^5$ prefactor
multiplies it, and at practical sample times the result is a ~60× error. Simpson
is exact through cubics, which brings the sliding-window estimate into agreement
with the algebraic ones.

A residual remains at second order and is *correct*: the integrand is
`y_kernel` (quadratic) × `y` (quadratic) = **quartic**, and Simpson is exact only
through cubics. At first order the integrand is quadratic, so Simpson is exact.
`tests/octave_sanity.m` encodes exactly this asymmetry — `1e-6` tolerance for
first order, `1e-3` for second.

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
| `model_order` | 1 or 2; picks the kernel pair and prefactor |
| `Ts` | sample time; the buffers shift once per `Ts` |
| `window_samples` | window length, **rounded up to even** for Simpson |
| `alpha` | ultra-local model gain; must match `mfc_command_block` |

## See also

- [[estimator-algebraic-2nd]] — the growing-window alternative
- [[iir-smoother]] — what to put after `F_hat` if it is noisy
- [[codegen-constraints]] — the buffer and struct-shape constraints
