# Model-Free Control (MFC) — SISO

Theory and implementation notes for the MFC SISO controller in this repo.

Written to render in **both** GitHub and Obsidian: `$…$` / `$$…$$` for math,
`[[wikilinks]]` between notes, no external images.

## Reading order

1. [[ultra-local-model]] — what MFC actually assumes, and what $F$ absorbs
2. [[control-law-coupled-vs-decoupled]] — where the stabilizing dynamics live
3. One estimator, depending on what you are doing:
   - [[estimator-algebraic-2nd]] — the default
   - [[estimator-algebraic-1st]] — and why it cannot be damped when coupled
   - [[estimator-sliding-window]] — finite memory, no time origin
   - [[riachy-trick]] — a second-order iPD with no derivative anywhere
4. [[iir-smoother]] — the one filter used in four places (two of them internal)
5. [[block-library-signal-flow]] — the blocks, their ports, and the two wiring rules
6. [[codegen-constraints]] — why the code is shaped the way it is

## The pipeline

```
setpoint ──► [1 input smoother] ──┬── sp_filt ──► (err) ─────────────┐
                                  └── dot_sp / ddot_sp ── ff ──┐     │
                                                               │     ▼
measurement ──┬──────────────────────────────► [2 F-hat estimator] ──┤
              │                                        │             │
              │                                   (3 F-hat filter)   │
              │                                        │             ▼
              │                                        └──► [4 command] ◄── [3 feedback]
              │                                                  │
              │                                                  ▼
              └──◄── plant ◄──────── [5 filter + saturation] ◄────┘
                                          │
                                     (unit delay) ──► u_prev, back to the estimator

Stage 3's F-hat filter is internal to the algebraic estimators. In a composed
loop stage 3's feedback is a stock Discrete PID (or Ground, when the estimator is
coupled) and stage 5 is a plain clamp on the command block; only `mfc_siso_core`
implements the full five stages in one block.
```

| Stage | Function | Block |
|---|---|---|
| 1 input smoother | `mfc_siso.ref_traj` → `mfc_iir_smoother` | `mfc_smoother_block` |
| 2 F-hat estimator | `mfc_fhat_algebraic_first_order` / `_second_order` / `mfc_fhat_sliding_window` | `mfc_fhat_alg{1,2}_{coupled,decoupled}_block`, `mfc_fhat_window_block` |
| 3 F-hat filter | `mfc_iir_smoother` (**inside** the algebraic estimators) | — (or `mfc_smoother_block` as a post-filter) |
| 3 feedback | `mfc_siso.feedback` | — (stock Discrete PID, or Ground when coupled) |
| 4 command | `mfc_siso.command` | `mfc_command_block` (clamp optional) |
| 5 command filter | `mfc_siso.limit` | — (core only; stock Discrete Filter otherwise) |
| *all of it* | `mfc_siso.step` | `mfc_siso_core` |

Stages with **—** deliberately have no block: either the math must stay inside
its estimator (3, the num/den filter), or Simulink already ships the block
(3 feedback, 5). See [[block-library-signal-flow]].

## Where the code lives

```
functions/    the math, plain functions + the mfc_siso namespace
  mfc_siso.m                          config, init, step, the 5 stages, window_kernel
  mfc_iir_smoother.m                  the shared 2nd-order IIR
  mfc_fhat_algebraic_first_order.m    growing-window estimator, 1st order
  mfc_fhat_algebraic_second_order.m   growing-window estimator, 2nd order
  mfc_fhat_sliding_window.m           finite-window Simpson estimator

blocks/       matlab.System wrappers -- no math, only ports/state/masks
  mfc_siso_core.m         the assembled all-in-one controller
  mfc_*_block.m           one block per pipeline stage

library/      build_mfc_lib.m generates the Simulink library from blocks/
examples/     val_mfc.m (all variants), val_mfc_composed.m (monolith vs stages)
tests/        see tests/README.md
```

**Two ways in.** Use `mfc_siso_core` when you want a controller that works. Use
the stage blocks when you want to take the method apart — they call exactly the
same functions, which [[block-library-signal-flow#Is the decomposition honest|is tested]].

## Notation

| Symbol | Code | Meaning |
|---|---|---|
| $y$ | `measure`, `y` | plant measurement |
| $y^*$ | `sp_filt` | filtered reference trajectory |
| $e = y - y^*$ | `err` | tracking error |
| $F$ | `F_hat` | lumped unknown dynamics (everything not $\alpha u$) |
| $\alpha$ | `alpha` | ultra-local model input gain (a design parameter) |
| $u$ | `u` | command |
| $T_s$ | `Ts` | sample time |
| $T_w$ | `kernel.Tw` | sliding-window length |
| $W$ | `window` | IIR smoother memory, in samples |
