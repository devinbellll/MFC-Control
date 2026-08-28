---
tags:
  - knowledge
  - control/pd
  - control/model-free
  - control/implementation
  - aero/phi-theory
  - simulation
  - tuning
aliases:
  - pitch break benchmark plant
  - softening Duffing pitch plant
  - SISO relative degree 2 nonlinear test plant
created: 2026-08-05
implementation: "models/plants/siso_pitch_break/pitchbreak_params.m"
model: "models/plants/siso_pitch_break/siso_pitch_break.slx"
---
# Nonlinear Pitch-Break Benchmark Plant for SISO Controller Development

A deliberately minimal SISO test plant for controller development: a softening
Duffing oscillator with quadratic damping, standing in for the pitch axis of a
tailless aircraft that loses static stability past a stall break. It is
relative degree 2 with no zero dynamics, has a **finite basin of attraction
bounded by a saddle**, and has **zero linear damping at the origin**. A 2-DOF
architecture — shaped reference plus exact feedforward, filtered PID feedback —
stabilises it. The achievable bandwidth is set by measurement noise and
actuator rate limit, never by the plant. The saddle itself is an artifact of
proportional feedback and can be passed by a model-free law; the hard envelope
limit is control authority, at $28.5^\circ$. This note fixes the equation, the
parameter meanings, the actuator model, the noise settings and the baseline
gains, so that later MFC results have a fixed comparison point.

> [!warning] These parameters are **chosen design values**, not measured or
> identified from any aircraft or dataset. The plant is a caricature with
> aerodynamically-motivated *structure*; the numbers were picked to place the
> break and the natural frequency at convenient values. Nothing here is an
> experimental fact.

## As built in `siso_pitch_break.slx`

Most of this note is *design intent*. The model implements a subset of it, and
the differences matter when reading a run. What is actually wired, top level:

```
ref (From Workspace: r, rdot, rddot, u_ff)
        │
        ├──────────────► Controller  ──── u ──► Plant ── y ──► Feedback ── yhat ──┐
Step (z, t=5 s) ─────────────────────────────────► Plant  ▲                       │
        └──────────────────────────────────────────────────────────────── yhat ◄──┘
```

- **Plant** (`Plant`): two inports, `u` and `z`. Two integrators realise
  $\ddot y = b u - c\dot y|\dot y| - (ky - \gamma y^3) + z$. Gains are `P.b`,
  `P.c`, `P.k`, `P.gam` — the equation below, exactly, plus the `z` term.
- **`z` is an acceleration disturbance**, driven by a Step of **10** (rad/s²) at
  **t = 5 s**. This is the disturbance-rejection test the note asks for further
  down; it is already in the model and it fires inside every run.
- **Feedback** (`Feedback`): $y$ + noise → 2nd-order Butterworth → `yhat`. The
  noise passes through a **Manual Switch** against Ground, so noise is toggled
  by hand in the diagram, not by a parameter. The pre-filter is a DSP
  **Analog Filter Design** block (Butterworth, lowpass, order 2, $W_{lo}$ =
  `P.wf`); its $W_{hi}$ = 80 field is inert for a lowpass.
- **Controller** is a **variant subsystem** with two label choices, `PID` and
  `HEOL`. **`HEOL` is the active choice as committed.** Both take
  `r, rdot, rddot, uff, yhat` and return `u`. Neither uses `rdot`/`rddot`
  directly — those enter only through the precomputed `u_ff`.
- **Scope** logs `r, u, yhat, y`. Solver is variable-step (auto/ode3),
  **StopTime = 10 s** while the reference timeseries only runs to
  `P.traj.tend = 8` s, so the last 2 s hold the final `ref` sample.

**PID variant.** $e = r - \hat y$ → stock continuous PID (`P.Kp`, `P.Ki`,
`P.Kd`, filter `P.N`) → gain $1/b$ → $+\,u_{ff}$ → $u$.

**HEOL variant.** The model-free law is applied to the *error* dynamics around
the feedforward, not to $y$:

$$
\varepsilon = \hat y - r,\qquad
u_{fb} = \frac{-\hat F - \mathrm{PID}(\varepsilon)}{b},\qquad
u = u_{ff} + u_{fb}
$$

with $\hat F$ from a `mfc_F_est_cont` subsystem reference (continuous algebraic,
decoupled: ports `y, u, alpha, time`). It is fed `y` = $\varepsilon$,
`u` = $u_{fb}$ (**the feedback component only, not the total $u$**),
`alpha` = a Constant `P.b`, and `time` from a global `Goto`/`From` off the
model Clock. Note the sign convention flips between the variants: PID uses
$r - \hat y$, HEOL uses $\hat y - r$, absorbed by the `--` sum ahead of $1/b$.

> [!warning] **Not implemented: the actuator.** There is no position
> saturation, no rate limiter and no first-order lag anywhere in the model —
> `u` reaches the plant directly. `P.tau`, `P.rate` and `P.umax` are defined in
> `pitchbreak_params.m` and **never read**. Everything in this note about the
> $\pm25^\circ$ clamp, the $300^\circ$/s rate limit, the $\tau=0.03$ s phase
> cost, anti-windup and rate-limit limit cycles is therefore *design intent for
> a later revision*, not a property of the current runs.

> [!warning] **Not implemented: discrete time.** The plant, filters, PID and
> estimator are all **continuous**. `P.Ts` feeds only the Band-Limited White
> Noise block's sample time. So "250 Hz controller rate" describes the noise
> colouring alone. Also `pitchbreak_params.m` does `P.Ts = Ts`, reading a bare
> `Ts` it never defines — **`Ts` must already exist in the base workspace** or
> the script errors.

> [!note] **Not implemented: anti-windup.** The PID block has `LimitOutput` off
> and `AntiWindupMode` = `none`. `P.Kb = P.Ki/P.Kd` is computed and unused.
> Harmless while nothing saturates, but it is the first thing to wire when the
> actuator goes in.

> [!note] **Derivative is on the error, not on $\hat y$.** Both variants use the
> stock 1-DOF PID block on a single error input, so its D term differentiates
> $e$ (filtered at `P.N`). The Implementation-notes section below prescribes
> derivative-on-measurement; the model does not do that yet. With a shaped
> septic reference the kick is small, which is why it has not bitten.

## What the plant is meant to represent

The pitch axis of a **swept tailless (flying-wing) aircraft on a single-DOF
pitch rig** — free to rotate about the CG, plunge constrained, at fixed
airspeed. Output $y$ is pitch/incidence angle in radians, input $u$ is elevon
deflection in radians.

Three physical facts are carried over; everything else is discarded.

**1. Second-order rotational dynamics.** Moment in, angular acceleration out.
This is what makes the relative degree exactly 2 and honest — the input enters
the moment equation, and the output is the angle two integrations later.

**2. The pitch break (softening stiffness).** On a swept tailless wing the
outboard panel stalls first. The centre of pressure moves *forward*, the
restoring pitching moment weakens and then reverses, and the aircraft pitches
**up into** the stall rather than nosing down out of it. A cubic $-\gamma y^3$
is the minimal smooth term with that shape: it opposes the linear restoring
term $ky$ and eventually overwhelms it.

**3. Quadratic aerodynamic damping.** Pitch damping comes from dynamic pressure
on the surfaces, which goes as velocity squared. Hence $c\,\dot y|\dot y|$
rather than $c\dot y$ — signed to always oppose motion, so it is passive and
can only remove energy.

What is deliberately **not** modelled: plunge and the phugoid, actuator
dynamics, the $\alpha$-dependence of control effectiveness, unsteady
aerodynamics, any coupling. Those are extensions, not part of the baseline.

## The equation

$$
\boxed{\;\ddot y + c\,\dot y\,|\dot y| + k\,y - \gamma\,y^{3} = b\,u\;}
$$

| Symbol | Value | Units | Meaning |
| --- | --- | --- | --- |
| $y$ | — | rad | pitch angle (the output) |
| $u$ | — | rad | elevon deflection (the input); $\pm 25^\circ$ is the design limit, **not enforced in the model** |
| $c$ | $0.8$ | rad$^{-1}$ | quadratic pitch-damping coefficient |
| $k$ | $60$ | s$^{-2}$ | linear restoring stiffness = static stability |
| $\gamma$ | $490$ | rad$^{-2}$s$^{-2}$ | softening (pitch-break) coefficient |
| $b$ | $-70$ | s$^{-2}$ | control effectiveness; negative = trailing-edge-down elevon gives nose-down moment |

$\gamma$ is not chosen directly. It is **set by where the break is wanted**:

$$
\gamma = \frac{k}{y_{\rm sad}^{2}}
$$

with $y_{\rm sad}$ the break angle. Here $y_{\rm sad} = 0.35$ rad ($20^\circ$)
gives $\gamma = 490$.

```matlab
% state x = [y; ydot]
dx = [x(2); -0.8*x(2)*abs(x(2)) - 60*x(1) + 490*x(1)^3 - 70*u];
```

The `Plant` subsystem adds a second inport `z` summed into $\ddot y$, so what it
integrates is $\ddot y = bu - c\dot y|\dot y| - ky + \gamma y^3 + z$.

## The break is at the saddle, not the inflection

Easy conflation, and one worth stating explicitly because it is a factor of
$\sqrt3$ and it was got wrong once already in developing this plant.

The restoring force is $f(y) = ky - \gamma y^3$, with potential

$$
V(y) = \tfrac{1}{2}ky^{2} - \tfrac{1}{4}\gamma y^{4}
$$

- **Inflection** — where $f'(y)=0$, i.e. stiffness stops growing:
  $y = \sqrt{k/3\gamma} = 0.202$ rad. The restoring force is still *positive*
  here. Nothing dramatic happens.
- **Saddle** — where $f(y)=0$, the hilltop of $V$:
  $$
  y_{\rm sad} = \sqrt{k/\gamma} = \sqrt{60/490} = 0.350\ \text{rad} = 20.0^\circ
  $$
  Inside this the trajectory is trapped in the well; outside it the cubic
  dominates and the state escapes.

> [!important] Divergence happens at the **saddle** $\sqrt{k/\gamma}$, where
> the restoring force changes sign — not at the inflection $\sqrt{k/3\gamma}$.

## Open-loop character

**Small-signal natural frequency.** Linearising about $y=0$ gives
$\ddot y + ky = bu$, so $k$ *is* $\omega_n^2$ by construction:

$$
\omega_n = \sqrt{k} = \sqrt{60} = 7.75\ \text{rad/s} = 1.23\ \text{Hz},\qquad T = 0.81\ \text{s}
$$

**No linear damping at the origin.** $\frac{d}{d\dot y}\big(c\dot y|\dot y|\big)
= 2c|\dot y|$, which is **zero at $\dot y = 0$**. The Jacobian at the
equilibrium is undamped and no linear $\zeta$ exists for this plant. Decay is
amplitude-dependent: large swings bleed energy quickly, small ones barely at
all, and the envelope decays as $1/t$ rather than exponentially.

> [!warning] Expect a response that looks well damped initially and then rings
> on stubbornly at low amplitude. This is the quadratic damping, not a tuning
> fault, and it will corrupt any settling-time metric defined on a tight band.

**Frequency falls with amplitude.** Local stiffness is $k - 3\gamma y^2$, so
oscillations at amplitude $A$ run slower. A describing-function estimate:

$$
\omega(A) \approx \sqrt{k - \tfrac{3}{4}\gamma A^{2}}
$$

| $A$ [rad] | $\omega(A)$ [rad/s] |
| --- | --- |
| $\to 0$ | 7.75 |
| 0.10 | 7.50 |
| 0.20 | 6.73 |
| 0.30 | 5.19 |
| $\to 0.350$ | $\to 0$ |

The period diverging at the separatrix is the standard saddle signature.

*(Both the $1/t$ envelope and the describing-function line above are analytical
readings of the equation, not simulation results — neither has been checked
against a run yet.)*

## Measurement noise

Noise is injected on $y$ with a Simulink **Band-Limited White Noise** block,
whose "Noise power" parameter is a PSD height $P$ [rad²·s] delivering

$$
\sigma^{2} = \frac{P}{T_s}\qquad\Longrightarrow\qquad P = \sigma^{2}T_s
$$

Block sample time is `P.Ts` and the seed is `P.seed = 23341`. Since the rest of
the model is continuous, `P.Ts` sets the noise correlation time only — there is
no controller rate to match it to yet. The noise is summed onto $y$ through a
**Manual Switch** (the other input is Ground), so a noise-free run is a switch
flip in the `Feedback` subsystem.

At $T_s = 0.004$ s (250 Hz):

| Case | $\sigma$ | $\sigma$ [rad] | Noise power $P$ |
| --- | --- | --- | --- |
| Clean estimator, good isolation | $0.1^\circ$ | $1.7\times10^{-3}$ | $1.2\times10^{-8}$ |
| **Baseline — typical small UAV** | $\mathbf{0.3^\circ}$ | $5.2\times10^{-3}$ | $\mathbf{1.1\times10^{-7}}$ |
| Prop vibration, poor isolation | $1.0^\circ$ | $1.7\times10^{-2}$ | $1.2\times10^{-6}$ |

> [!note] White noise is a **modelling convenience here, not a claim about
> EKF2**. Real estimator output error is time-correlated; a white model at
> matched RMS overstates the high-frequency content the derivative filter must
> reject and understates slow drift. Adequate for controller development,
> inadequate as a noise budget. See [[Gauss-Markov Model of EKF2 Output Error]].

## Why the naive filtered derivative is unusable

The derivative path binds. A raw backward difference gives

$$
\sigma_{\dot y} = \frac{\sigma\sqrt2}{T_s} = 1.84\ \text{rad/s}
$$

so a filtered derivative $\dfrac{Ns}{s+N}$ is mandatory. But the filter alone is
not enough, and the reason is a spectrum problem rather than an amplitude one.

$\left|\frac{Ns}{s+N}\right| \to N$ as $\omega\to\infty$: the filtered derivative
**does not roll off**. Above the corner it flattens at gain $N$ and stays there
all the way to Nyquist ($\pi/T_s = 785$ rad/s). So the rate noise is
$\sigma_{\dot y_f}\approx N\sigma \approx 0.16$ rad/s at $N=30$, but its content
is flat and top-heavy across two and a half decades.

The commanded elevon RMS that results ($\approx 2.2^\circ$ at $K_d=17$) is
harmless in itself. What is not harmless is its **rate**. For a flat spectrum
bandlimited at $\omega_N$,

$$
\sigma_{\dot u} \approx \frac{\omega_N}{\sqrt3}\,\sigma_u
= 453 \times 0.038 \approx 17\ \text{rad/s} \approx 970\ ^\circ/\text{s}
$$

> [!warning] That is roughly $3\times$ the rate limit of a realistic servo,
> **from noise alone**. A rate limiter in permanent saturation is not a benign
> clip: it behaves as a large amplitude-dependent phase lag and is a classic
> limit-cycle mechanism. Lowering the gains does not fix this — the fix is to
> roll the measurement off *before* differentiating it.

## The actuator model

> [!warning] **Design intent only — no actuator exists in
> `siso_pitch_break.slx`.** `P.tau`, `P.rate`, `P.umax` are declared and unused;
> the commanded $u$ is the applied $u$. Read this section, the phase budget and
> the rate checks as the specification for the next revision.

A plausible small-UAV elevon servo. Chosen values, not identified from hardware.

| Property | Value | Meaning |
| --- | --- | --- |
| First-order lag $\tau$ | $0.030$ s | $\omega_{\rm act} = 33$ rad/s (5.3 Hz) |
| Rate limit $R$ | $300\ ^\circ$/s $= 5.24$ rad/s | $\approx 0.2$ s per $60^\circ$, loaded |
| Position limit | $\pm25^\circ$ | hard mechanical stop |

$$
\dot u = \mathrm{sat}_{\pm R}\!\left(
\frac{\mathrm{sat}_{\pm 25^\circ}(u_c) - u}{\tau}\right)
$$

**Ordering matters**: position clamp on the *command*, rate limit on the
*slew*, lag last. Applied in that order the rate limiter cannot drive the
internal state past the mechanical stops.

The actuator also sets a hard ceiling on useful loop bandwidth, conventionally
$\omega_n^{\rm CL} \lesssim \omega_{\rm act}/3 \approx 11$ rad/s.

## Measurement pre-filter

A second-order Butterworth on $y$, ahead of the whole controller:

$$
G_f(s) = \frac{\omega_f^{2}}{s^{2} + \sqrt2\,\omega_f s + \omega_f^{2}},
\qquad \omega_f = 60\ \text{rad/s}
$$

$\omega_f$ is placed at $2\times$ the derivative corner $N$ (so it does not
fight the D path) and $\approx0.55\times$ the actuator (so it, not the servo, is
the dominant roll-off). Everything downstream is now bandlimited at 60 rad/s
instead of 785.

Noise-equivalent bandwidth of a 2nd-order Butterworth is $1.11\,\omega_f = 67$
rad/s, so

$$
\sigma_{y_f} = \sigma\sqrt{\frac{67}{785}} = 0.29\,\sigma
= 1.5\times10^{-3}\ \text{rad} = 0.086^\circ
$$

A $3.4\times$ reduction in measured noise — but the reduction in *rate* is the
point, and that comes from the roll-off rather than the amplitude.

> [!important] Two filters, two jobs. $N$ limits how much noise the derivative
> **amplifies**; $\omega_f$ limits how **fast** the resulting command can move.
> Only the second one protects the servo.

## Baseline controller — filtered PID

$$
u_{\rm fb} = -\frac{1}{|b|}\Big(K_p e + K_i\!\!\int\!\! e\,dt + K_d\,\dot e_f\Big),
\qquad e = r - y_f,\qquad \dot e_f = -\frac{Ns}{s+N}\,y_f
$$

The $1/|b|$ is sign/scale bookkeeping for the negative control effectiveness,
not feedback linearisation — $b$ is constant in this baseline. The model writes
this as a single gain `1/P.b` on the PID output (with $b<0$ the two forms agree,
since it pairs with $e = r-\hat y$); there is no separate sign block.

In the model this is a stock **continuous** PID block, parallel form, `P`/`I`/`D`
= `P.Kp`/`P.Ki`/`P.Kd`, filter coefficient `N` = `P.N`, filter on, output
limiting off, anti-windup none — and its derivative acts on the error, not on
$\hat y$.

| Gain | Value | Role |
| --- | --- | --- |
| $K_p$ | 60 | sets $\omega_n = \sqrt{k+K_p} = \sqrt{120} = 11.0$ rad/s |
| $K_i$ | 200 | integral zero at 3.3 rad/s |
| $K_d$ | 15 | $\approx 2\zeta\omega_n$ at $\zeta = 0.7$ |
| $N$ | 30 rad/s | derivative filter corner |
| $\omega_f$ | 60 rad/s | measurement pre-filter |

Retuned down from an earlier noise-only design ($K_p=84$, $K_i=300$, $K_d=17$),
which put $\omega_n$ at 12 rad/s — above the actuator ceiling of 11.

**Small-signal closed loop.** Proportional gain adds directly to stiffness, so
the characteristic polynomial (ignoring filter and actuator poles) is

$$
s^{3} + K_d s^{2} + (k+K_p)s + K_i = s^{3} + 15s^{2} + 120s + 200
$$

which factors into a real pole at $s \approx -2.17$ and a complex pair at
$\omega = 9.6$ rad/s, $\zeta = 0.67$. ✔

**Routh check.** Stability of the cubic requires $K_d(k+K_p) > K_i$:
$15\times120 = 1800 \gg 200$, a margin of $9\times$. ✔

**Noise check.**

- D path: $K_d\,\sigma_{\dot y_f}/|b| = 15\times0.045/70 = 0.55^\circ$ RMS
- P path: $K_p\,\sigma_{y_f}/|b| = 60\times1.5\times10^{-3}/70 = 0.07^\circ$
- Total $\approx 0.56^\circ$ RMS, down from $2.2^\circ$. ✔

**Rate check** — the one that motivated all of this:

$$
\sigma_{\dot u} \approx \frac{\omega_f}{\sqrt3}\,\sigma_u
\approx 34.6 \times 9.7\times10^{-3} \approx 0.34\ \text{rad/s} \approx 19\ ^\circ/\text{s}
$$

About **6% of the rate limit**, with $3\sigma$ peaks near $58\ ^\circ$/s (19%).
The servo is no longer rate-saturated on noise. ✔

> [!note] Useful rule: keep $\sigma_{\dot u}$ under $\sim$50% of $R$ in normal
> operation. A rate limiter that engages even intermittently injects a
> describing-function phase lag right at crossover.

## The phase budget is what you pay for smoothness

> [!note] The actuator row below is not in the model yet (no lag block), so the
> as-built parasitic lag is the pre-filter plus the derivative filter,
> $\approx25^\circ$, and the loop will look better damped than this table
> predicts.

Every one of the above filters costs phase at crossover ($\omega \approx 11$ rad/s):

| Source | Lag at 11 rad/s |
| --- | --- |
| Actuator lag, $\tau = 0.03$ | $\arctan(0.33) = 18.4^\circ$ |
| Pre-filter, $\omega_f = 60$ | $15.0^\circ$ |
| Derivative filter, $N = 30$ | $\approx10^\circ$ effective |
| **Total parasitic** | $\approx\mathbf{43^\circ}$ |

So the $\zeta = 0.7$ from the cubic is optimistic; expect $\zeta \approx
0.45$–$0.5$ once the servo and filters are in the loop. Verify with `margin()`
on the fully linearised loop rather than trusting the cubic. If phase margin
lands below $40^\circ$, take $K_p$ to 45 and $\omega_f$ to 80.

## Reference shaping and feedforward

The reference is passed through its own filter before reaching the controller.
This is a **2-DOF architecture**: the reference defines what the aircraft is
*asked* to do, and the feedback loop handles only the difference.

**A shaped polynomial, not a reference filter.** An earlier draft used a
third-order critically-damped filter. A precomputed rest-to-rest polynomial does
the same job with less machinery and gives exact derivatives rather than filter
states. Normalised on $\tau = (t-t_0)/T$:

$$
\begin{aligned}
\text{quintic:}\quad & s(\tau) = 6\tau^5 - 15\tau^4 + 10\tau^3 \\
\text{septic:}\quad  & s(\tau) = -20\tau^7 + 70\tau^6 - 84\tau^5 + 35\tau^4
\end{aligned}
$$

with $r = y_0 + \Delta y\,s(\tau)$, $\dot r = \Delta y\,s'/T$,
$\ddot r = \Delta y\,s''/T^2$.

**Septic by default.** Since $u_{ff}\propto\ddot r$, the elevon *rate* follows
jerk. A quintic has a jerk step at both endpoints, so $u_{ff}$ starts with a
rate discontinuity — the same defect the second-order filter had, relocated. The
septic zeroes jerk at both ends and $u_{ff}$ eases in from zero rate.

| Peak factors | $\dot r$ | $\ddot r$ | $\dddot r$ |
| --- | --- | --- | --- |
| quintic | $1.875\,\Delta y/T$ | $5.77\,\Delta y/T^2$ | $60\,\Delta y/T^3$ |
| septic | $2.188\,\Delta y/T$ | $7.51\,\Delta y/T^2$ | $52.5\,\Delta y/T^3$ |

17% more peak velocity for a much cleaner command. Sizing $T$ is then explicit
rather than implicit in a filter bandwidth.

**As built.** `rest_to_rest` in `pitchbreak_params.m` builds the septic
(`P.traj.order = 7`) on the `P.Ts` grid: $y_0 = 0$, $y_1 = 10^\circ$,
$t_0 = 1$ s, $T = 2$ s, out to `tend` = 8 s. The four columns
$[r,\ \dot r,\ \ddot r,\ u_{ff}]$ are packed into the `ref` timeseries and read
by a From Workspace block, then demuxed. So $u_{ff}$ is **precomputed offline**,
not evaluated in the diagram — changing a plant gain in `P` without re-running
the script silently leaves a stale feedforward.

**Feedforward is exact.** Relative degree 2 with no zero dynamics means the
plant inverts along the reference:

$$
\boxed{\;u_{ff} = \frac{1}{b}\Big(\ddot r + c\,\dot r|\dot r| + k r - \gamma r^{3}\Big)\;}
\qquad u_c = u_{ff} + u_{\rm fb}
$$

This leaves feedback handling only disturbance and model error. It is the
flatness-based 2-DOF structure on a toy problem, and the MFC comparison will
want the same architecture.

The command is dominated by the **static trim term** $kr/b$, not the transient —
the shaping is what made the transient small. For a $10^\circ$ hold the trim
alone is $6.4^\circ$ of elevon.

> [!warning] A slow, smooth reference makes a badly damped loop look fine,
> because it never excites the modes. Keep a **disturbance-rejection test** — an
> impulse on $\ddot y$, or a step injected at $u$ — as the real damping check.
> That path does not pass through the reference at all. **This is already in the
> model**: the `z` inport on `Plant`, driven by a Step of 10 rad/s² at t = 5 s,
> i.e. 3 s after the manoeuvre completes.

## Control authority, not the saddle, is the true envelope limit

The saddle $y_{\rm sad}^{\rm CL}$ is an **artifact of proportional feedback**,
not a property of the plant: it exists because $K_p$ contributes a fixed linear
stiffness that the cubic eventually outruns. A control law with no fixed
stiffness — model-free control, where the cubic lands in $F$ and is cancelled by
$\hat F$ — has error dynamics

$$
\ddot e + K_d\dot e + K_p e = -(F - \hat F)
$$

in which $\gamma$ does not appear. Its basin is set by estimation error, not by
stiffness balance, so it is not bounded by the saddle.

What no control law escapes is **statics**. Holding any equilibrium $y^*$
requires

$$
u^*(y^*) = \frac{k y^* - \gamma y^{*3}}{b}
$$

| $y^*$ | $\vert u^*\vert$ | |
| --- | --- | --- |
| 0.10 rad ($5.7^\circ$) | $4.5^\circ$ | |
| 0.202 rad ($11.6^\circ$) | $6.6^\circ$ | **peak demand**, then reverses |
| 0.350 rad ($20.0^\circ$) | $0^\circ$ | the saddle — trims hands-off |
| 0.45 rad ($25.8^\circ$) | $14.4^\circ$ | **opposite sign** |
| 0.498 rad ($28.5^\circ$) | $25^\circ$ | saturated |

> [!note] The **inflection** $\sqrt{k/3\gamma}$ earns its keep after all: it is
> where $du^*/dy$ vanishes, i.e. where required trim peaks and **reverses**.
> Below it, more elevon holds more incidence; above it, progressively less; past
> the saddle, opposite elevon is needed to hold the nose down. That reversal
> *is* the pitch break as a pilot would meet it.

Solving $|ky - \gamma y^3| = |b|u_{\max}$ gives

$$
\boxed{\;y_{\max} = 0.498\ \text{rad} = 28.5^\circ\;}
$$

> [!warning] $y_{\max} = 28.5^\circ$ and the PD saddle $y_{\rm sad}^{\rm CL} =
> 28.4^\circ$ are numerically almost identical at the baseline gains, because
> $K_p y \approx |b|u_{\max}$ near that point is roughly the condition for both.
> **In this parameterisation MFC buys essentially no steady-state envelope over
> PD.** To make the comparison show anything, either raise $u_{\max}$ or lower
> $K_p$ so the PD basin is visibly the smaller of the two. This is a coincidence
> of the chosen numbers, not a general result.

Where MFC should still win, unchanged by the above:

- **Transient excursions past the saddle.** Passing through $25^\circ$ and
  returning is not the same as holding it. PD diverges once past its saddle
  regardless of available authority; MFC tracks through.
- **State-dependent effectiveness $b(y)$** — the extension a fixed $1/b$
  inversion cannot handle.

**The binding constraint is estimator bandwidth against escape rate.** Past the
saddle the linearisation has an unstable eigenvalue

$$
\lambda = \sqrt{3\gamma y^{2} - k}
$$

$13.2$ rad/s at $y=0.40$, $17.5$ rad/s at $y=0.50$. A window-based $\hat F$ lags
by roughly $T_w/2$, so tracking requires $2/T_w \gg \lambda$ — of order
$T_w \lesssim 20$ ms, i.e. **5 samples at 250 Hz**. Short windows are exactly
where measurement noise bites, since $\hat F$ needs a second derivative of $y$.

> [!important] The real experiment is not "can MFC beat the saddle" — noise-free,
> it trivially can. It is **how far past the saddle the estimator holds before
> noise-driven $\hat F$ error exceeds the divergence rate**: a sweep of window
> length against $\sigma_y$ against target $y$.

For the baseline, set $\alpha = b = -70$ — implemented as a Constant `P.b` on the
estimator's `alpha` port, so the mismatch study is a one-block edit. Note also
that the as-built estimator sees the tracking error $\varepsilon = \hat y - r$
and the feedback command $u_{fb}$, not $y$ and the total $u$: $\hat F$ is the
residual of the *error* dynamics after the feedforward has taken the nominal
plant out. The true gain is known here, so the
baseline should isolate estimator error from gain error; the $\alpha$-mismatch
study belongs separately — see
[[Sliding-Window Estimator Gain and the Alpha Mismatch]].

## Proportional gain buys basin, at a square-root rate

Because $K_p$ adds to stiffness while $\gamma$ is untouched, the closed-loop
saddle moves out:

$$
y_{\rm sad}^{\rm CL} = \sqrt{\frac{k+K_p}{\gamma}}
$$

| $K_p$ | $y_{\rm sad}^{\rm CL}$ | |
| --- | --- | --- |
| 0 | 0.350 rad = $20.0^\circ$ | open loop |
| 60 | 0.495 rad = $28.4^\circ$ | **baseline, actuator-limited** |
| 84 | 0.542 rad = $31.1^\circ$ | noise-limited only, no actuator |
| 300 | 0.870 rad = $49.8^\circ$ | ideal, neither constraint |

Three consequences worth keeping:

- **Damping does not help.** $K_d$ cannot move $y_{\rm sad}^{\rm CL}$ — the
  saddle is a property of the stiffness balance alone. Only $K_p$ moves it.
- **Realism is paid for in stall margin.** Noise costs $19^\circ$ of basin
  (300 → 84); actuator dynamics cost a further $2.7^\circ$ (84 → 60). A clean
  two-entry ledger for a report — though at the baseline gains the ledger runs
  out exactly where control authority does, $28.5^\circ$.
- Stabilisation here is **local**, never global. The characteristic failure is a
  loop that tracks beautifully and then falls off a cliff at a finite angle.

## Implementation notes

The four bullets below are prescriptions; as noted at the top, the first two are
**not yet honoured by the model** (stock PID on the error, no saturation and no
anti-windup).

- **Derivative on the filtered measurement**, never on the error. With $K_d=15$
  a reference step differentiated directly would slam the elevon; and the
  reference filter already supplies $\dot r$ if a rate term on $r$ is wanted.
- **Anti-windup must see the true surface**, i.e. post-position-limit *and*
  post-rate-limit. With a rate limiter present, command and actual surface can
  diverge for tens of milliseconds even when nothing is position-saturated, and
  the integrator will wind on that gap. Back-calculation with
  $T_t \approx K_d/K_i \approx 0.08$ s, or clamping. Same principle as feeding
  **applied** rather than commanded $u$ to an estimator — see
  [[Estimator Windup Under Actuator Saturation]].
- **The integrator carries the low-amplitude tail.** Quadratic damping vanishes
  as $\dot y\to0$, so near-zero error is worked off almost entirely by $K_i$.
  Expect a slow sub-degree crawl that looks worse than the design $\zeta$ implies.
- **Hold the environment fixed across controllers.** $\omega_f$, $\tau$, $R$ and
  the reference filter are part of the test rig, not the design. For an
  MFC-vs-PID comparison to mean anything, only the control law should vary — the
  estimator window length plays the same bandwidth-vs-noise role that $N$ and
  $\omega_f$ play here.

## Extensions, in order of usefulness

1. **State-dependent control effectiveness.** Elevons sit on the outboard
   trailing edge, exactly where separation starts, so effectiveness fades into
   the break: $b(y) = b_0(1 - \eta\,y^2/y_{\rm sad}^2)$, clipped away from zero
   so relative degree 2 survives. This is the case a fixed $1/b$ PID cannot
   handle and an ultra-local estimator should. Note it also breaks the exact
   feedforward inversion above, which then needs $b(r)$ rather than $b$.
2. **Escape probability under noise.** Sweep $P$ across two decades against
   initial condition relative to $y_{\rm sad}^{\rm CL}$. Near the hilltop noise
   does not merely degrade tracking, it randomly kicks the state over, making
   the basin boundary probabilistic. Likely a sharper discriminator between
   controllers than any step-response metric.
3. **Model-free control past the saddle**, with the estimator window swept
   against $\sigma_y$ and target incidence. Raise $u_{\max}$ or lower $K_p$
   first, or the PD baseline and the authority limit coincide and the comparison
   shows nothing.
4. **Rate-limit-induced limit cycle.** Deliberately re-raise $\omega_f$ or $K_d$
   until the servo rate-saturates on noise, and characterise the resulting
   oscillation. Establishes the failure boundary the current design sits inside.

> [!todo] `siso_pitch_break.slx` exists and runs, but no results from it are
> recorded here, and its missing actuator means it cannot yet test the claims
> that depend on one. Every gain, pole location, noise
> figure, trim value and rate estimate in this note is an analytical prediction
> and wants a run to confirm it — particularly the phase-budget claim and the
> $T_w \lesssim 20$ ms estimator bound, which are the two most likely to be wrong.

> [!todo] `report_margins` in `pitchbreak_params.m` currently warns at 80% of the
> closed-loop saddle. It should also hard-error on a reference beyond
> $y_{\max} = 0.498$ rad: that is a physics limit, not a tuning preference, and
> no controller can be blamed for failing to track it.

## Related notes
- [[Estimator Windup Under Actuator Saturation]]
- [[Gauss-Markov Model of EKF2 Output Error]]
- [[Ultra-Local Model with Actuator Lag]]
- [[Sliding-Window Estimator Gain and the Alpha Mismatch]]
