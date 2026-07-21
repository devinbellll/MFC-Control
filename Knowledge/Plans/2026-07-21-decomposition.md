# Decompose MFC SISO into a hackable block library + Knowledge docs

## Context

`functions/` currently reads like a finished MathWorks add-on: one manager System
object (`mfc_siso_core.m`) sitting on a 3-file lifecycle (`mfc_siso_config.m` →
`mfc_siso_init.m` → `mfc_siso_step.m`) plus two helper-math files
(`mfc_window_kernel.m`, `mfc_iir_smoother.m`) and three estimators. Everything a
variant can be is decided by mask dropdowns funnelled through one config struct
and one `if/else` dispatch inside `mfc_siso_step.m:83-111`.

That is great for shipping and bad for research. To try a new F-hat estimator, a
different smoother placement, or a different feedback law right now you have to
edit the dispatch, add a config field, add a state field, add a DiscreteState
property, and update `getDiscreteStateSpecificationImpl` — five files for one
experiment. The pipeline stages (input smoother → F-hat estimator → F-hat
filtering → feedback law → command computation → command filter) exist only as
numbered comment blocks inside one function; they cannot be re-ordered, swapped,
or probed individually.

The goal is to expose those stages as first-class, individually usable Simulink
blocks and MATLAB objects, collapse the lifecycle scaffolding into a single hub
file, and write the theory + implementation down in `Knowledge/` — **without
regressing** the code-generation-compatible `mfc_siso_core` block, which stays as
the assembled all-in-one controller.

Note also: `tests/test_estimators.m` and `tests/README.md` are **stale** — they
call `mfc_siso_non_algebraic`, `mfc_siso_decoupled`, `mfc_siso_run` and
`mfc_siso_core('int_window',…,'kp',…)`, all deleted or renamed in commit
`228cf5d`. Nothing in `tests/` currently runs. Fixing that is part of this work,
and the regression harness is the safety net for the whole refactor.

## Constraints

- `mfc_siso_core.m` must keep working, keep its mask, and keep compiling under
  MATLAB Coder. Its codegen-driven quirks are load-bearing and documented in
  place: single concrete `cfg` type (`mfc_siso_config.m:115-126`), identical
  `dbg` struct shape across all three estimators, `resetImpl` typing the discrete
  states with a compile-time-constant buffer length
  (`mfc_siso_core.m:273-292`), `circshift` on 1×1 placeholder buffers
  (`mfc_fhat_sliding_window.m:48-55`).
- No MATLAB in this environment — only Octave (`/usr/bin/octave`). Anything
  requiring Simulink or `matlab.System` must be shipped as a script you run once.
- Repo convention: models are text-format OPC `.mdl`, not `.slx`.
- The numeric behaviour of every existing variant must be bit-identical after the
  refactor.

## Design

### Layer 0 — kernels (plain functions, math unchanged)

Stay exactly as they are; they are the single source of truth that both the
monolith and the new stage blocks call.

- `functions/mfc_iir_smoother.m`
- `functions/mfc_fhat_algebraic_first_order.m`
- `functions/mfc_fhat_algebraic_second_order.m`
- `functions/mfc_fhat_sliding_window.m`

### Layer 1 — hub file (the simplification)

Fold `mfc_siso_config.m` + `mfc_siso_init.m` + `mfc_siso_step.m` +
`mfc_window_kernel.m` into **one** file `functions/mfc_siso.m`, a `classdef
mfc_siso` with only static methods and no properties (MATLAB permits one public
function per `.m` file, so a classdef namespace is the only way to get several
entry points into one file):

```
mfc_siso.config(varargin)              <- was mfc_siso_config.m
mfc_siso.init(cfg)                     <- was mfc_siso_init.m
mfc_siso.step(sp,y,t,u_prev,alpha,cfg,state)   <- was mfc_siso_step.m
mfc_siso.window_kernel(order,win,Ts)   <- was mfc_window_kernel.m
mfc_siso.ref_traj(sp, state, cfg)      <- extracted stage 1
mfc_siso.feedback(err, state, cfg)     <- extracted stage 3 (fb, int_err, dot_err)
mfc_siso.command(F_hat, ff, fb, alpha) <- extracted stage 3b
mfc_siso.limit(u_raw, u_prev, cfg)     <- extracted stage 4 (EMA + sat + windup)
```

`mfc_siso.step` becomes a short readable pipeline that calls the four stage
statics plus one estimator — the same numbers, but each stage now has a name and
a call site that a stage block can reuse.

`functions/` goes 8 files → 5.

**Risk to retire first:** MATLAB Coder supports static-method calls on a
codegen-compatible classdef, and Octave's classdef support is partial. Step 1
below verifies both before any other work; the fallback if either balks is the
`+mfc/` package folder (same split, plain function files, zero risk).

### Layer 2 — stage blocks (`blocks/`)

Each is a thin `matlab.System` wrapper (~60–120 lines) over a Layer-0/1 function.
No new math anywhere.

| Block | Wraps | In → Out |
|---|---|---|
| `mfc_smoother_block` | `mfc_iir_smoother` | `x` → `x_filt` (+ optional `dot_x`, `ddot_x`) |
| `mfc_fhat_alg1_block` | `mfc_fhat_algebraic_first_order` | `z, u_prev, [alpha], t` → `F_hat, valid` (+ optional `num_raw, den_raw`) |
| `mfc_fhat_alg2_block` | `mfc_fhat_algebraic_second_order` | same, plus `a_fold` param |
| `mfc_fhat_window_block` | `mfc_fhat_sliding_window` + `mfc_siso.window_kernel` in `setupImpl` | `y, u_prev, [alpha], t` → `F_hat, valid` (+ optional `integral`) |
| `mfc_fhat_divide_block` | the divide + hold logic | `num, den, t` → `F_hat, valid` |
| `mfc_feedback_block` | `mfc_siso.feedback` | `err, [freeze]` → `fb, int_err, dot_err` |
| `mfc_command_block` | `mfc_siso.command` | `F_hat, ff, fb, [alpha]` → `u_raw` |
| `mfc_command_filter_block` | `mfc_siso.limit` | `u_raw, [u_prev]` → `u, sat` |
| `mfc_siso_core` | `mfc_siso.step` | unchanged ports, unchanged mask |

`mfc_smoother_block` is deliberately one class in four roles — reference-trajectory
filter, numerator filter, denominator filter, F-hat post-filter — so tweaking the
smoother math is one edit that propagates everywhere.

**F-hat filtering, dissected mode.** The algebraic estimators smooth numerator and
denominator *separately, before* dividing (`mfc_fhat_algebraic_second_order.m:78-87`)
— that is mathematically essential and stays the default, self-contained
behaviour. Each algebraic block also gets a `expose_raw` checkbox and an
`internal_filter` checkbox, so the chain can be rebuilt externally and mutated:

```
default:     [alg2] ──────────────────────────────────► F_hat

dissected:   [alg2 raw] ─ num_raw ─► [smoother] ─┐
                        ─ den_raw ─► [smoother] ─┴─► [divide+hold] ─► F_hat

post-filter: F_hat ─► [smoother] ─► F_hat_filt      (mainly for the sliding
                                                     window, which has none)
```

**Two consequences of splitting, to be documented not hidden:**

1. *Loop closure.* The monolith carries `u_prev` internally (`state.u_km1`). A
   loop assembled from stage blocks must route `u` back through an explicit Unit
   Delay to the estimator and command-filter inputs. Same signal, now visible.
2. *Anti-windup.* The monolith freezes the integrator in the same sample it
   saturates (`mfc_siso_step.m:131-137`). Split apart, `mfc_command_filter_block`
   emits a `sat` flag into `mfc_feedback_block`'s optional `freeze` port, which
   needs a Unit Delay to avoid an algebraic loop — so the composed loop freezes
   **one sample later** than the monolith. This is an intentional, documented
   difference; the monolith remains the reference implementation.

### Layer 3 — library and examples

- `library/build_mfc_lib.m` — run once in MATLAB; programmatically creates,
  masks, positions and saves `library/mfc_lib.mdl` (all nine blocks + a
  pre-wired "assembled loop" demo subsystem showing the dissected chain).
- `library/slblocks.m` — puts `mfc_lib` in the Simulink Library Browser.
- `examples/val_mfc.m` — moved out of `functions/` (it is a script, not a
  function); unchanged otherwise.
- `examples/val_mfc_composed.m` — new: builds the same controller twice, once as
  `mfc_siso_core` and once from stage objects in plain MATLAB code, and overlays
  them. This is the proof that the decomposition is faithful and the template for
  how you tweak guts without Simulink.
- `setup.m` — add `blocks`, `library`, `examples` to the path.

### Layer 4 — `Knowledge/`

GitHub- and Obsidian-compatible Markdown: `$…$` / `$$…$$` math (renders in both),
`[[wikilinks]]` between docs, fenced ```matlab blocks, relative image-free.

```
Knowledge/
  README.md                          index, block map, reading order
  ultra-local-model.md               F + alpha*u, 1st vs 2nd order, what F absorbs
  estimator-algebraic-1st.md         Laplace annihilator derivation → discretization → mfc_fhat_algebraic_first_order.m walkthrough
  estimator-algebraic-2nd.md         same, twice-differentiated, incl. a_fold/b_fold algebra
  estimator-sliding-window.md        Eq.11/Eq.16 integral kernels, Simpson vs trapezoid (why 60/Tw^5 forbids trapezoid)
  iir-smoother.md                    unity-DC critically-damped 2nd-order IIR, pole at W/(W+1)
  control-law-coupled-vs-decoupled.md  the fold-vs-explicit table, why 1st-order coupled cannot fold Kd
  block-library-signal-flow.md       the nine blocks, ports, params, delay/anti-windup rules above
  codegen-constraints.md             every codegen quirk and why the code looks like that
  Plans/
    2026-07-21-decomposition.md      this plan, checked in
```

Each algorithm doc follows: **Theory** (LaTeX) → **Discretization** → **MATLAB
implementation** (annotated excerpt, file:line pointers, parameter meanings) →
**Gotchas** → **See also** wikilinks. Most of this content already exists as
excellent header comments in the source and gets promoted, not rewritten.

## Execution order

1. **Retire the hub-file risk.** Write a 20-line Octave probe for
   `classdef` + static methods, and confirm the Coder pattern in MATLAB. If
   either fails, switch Layer 1 to a `+mfc/` package folder and continue
   unchanged.
2. **Golden traces (before touching anything).** New `tests/golden_capture.m`
   drives today's `mfc_siso_step` over all six supported variants on the
   `val_mfc.m` drone plant and writes `tests/golden/*.csv`
   (`u, F_hat, sp_filt, err, u_raw, valid`). Runs in Octave — `mfc_siso_step` is a
   plain function today. This is the no-regression contract.
3. **Layer 1 fold.** Create `functions/mfc_siso.m`; delete `mfc_siso_config.m`,
   `mfc_siso_init.m`, `mfc_siso_step.m`, `mfc_window_kernel.m`. Repoint
   `mfc_siso_core.m` (`buildConfig`, `stepImpl`) and `mfc_fhat_sliding_window`'s
   callers. Re-run golden traces → must match to `eps`.
4. **Layer 2 blocks.** Add `blocks/` one class at a time, each a wrapper over an
   already-tested function. `mfc_siso_core.m` moves to `blocks/` too.
5. **Layer 3.** `library/build_mfc_lib.m`, `library/slblocks.m`, move + add
   examples, update `setup.m`.
6. **Fix `tests/`.** Rewrite `tests/test_estimators.m` against the real current
   API (it currently calls three deleted classes), rewrite `tests/octave_sanity.m`
   to call `functions/mfc_fhat_*` directly instead of carrying a duplicate copy of
   the math, and rewrite `tests/README.md`. Add `tests/test_golden.m` (Octave) and
   `tests/test_composed.m` (MATLAB, stage objects vs monolith).
7. **Layer 4 docs**, written last so they describe what actually shipped.

## Files touched

| Action | Path |
|---|---|
| new | `functions/mfc_siso.m` |
| delete | `functions/mfc_siso_config.m`, `mfc_siso_init.m`, `mfc_siso_step.m`, `mfc_window_kernel.m` |
| unchanged | `functions/mfc_iir_smoother.m`, `functions/mfc_fhat_*.m` (3) |
| move + repoint | `functions/mfc_siso_core.m` → `blocks/mfc_siso_core.m` |
| new (8) | `blocks/mfc_smoother_block.m`, `mfc_fhat_alg1_block.m`, `mfc_fhat_alg2_block.m`, `mfc_fhat_window_block.m`, `mfc_fhat_divide_block.m`, `mfc_feedback_block.m`, `mfc_command_block.m`, `mfc_command_filter_block.m` |
| new | `library/build_mfc_lib.m`, `library/slblocks.m` |
| move | `functions/val_mfc.m` → `examples/val_mfc.m` |
| new | `examples/val_mfc_composed.m` |
| rewrite | `tests/test_estimators.m`, `tests/octave_sanity.m`, `tests/README.md` |
| new | `tests/golden_capture.m`, `tests/test_golden.m`, `tests/test_composed.m`, `tests/golden/*.csv` |
| edit | `setup.m` (path additions) |
| new | `Knowledge/` (9 docs + `Plans/`) |

## Verification

**In this sandbox (Octave, no Simulink):**

```bash
octave --no-gui -q tests/golden_capture.m   # step 2, before refactor
octave --no-gui -q tests/test_golden.m      # after every subsequent step
octave --no-gui -q tests/octave_sanity.m    # estimator math vs analytic F
```

`test_golden.m` must report zero mismatch across all six variants — that is the
gate on Layer 1 and the definition of "did not regress".

**In your MATLAB (I cannot run these):**

```matlab
>> setup
>> test_estimators        % System objects via step(), closed loop
>> test_composed          % stage objects vs mfc_siso_core, same trace
>> val_mfc                % the six-variant overlay plot, unchanged from today
>> val_mfc_composed       % monolith vs hand-assembled loop, overlaid
>> build_mfc_lib          % generates library/mfc_lib.mdl
```

**Codegen (the explicit no-regression check on the shipped block):**

```matlab
>> cfg = coder.config('lib');
>> codegen -config cfg -report mfc_siso_core -args {0,0,0}
```

must still succeed for at least the 2nd-order coupled algebraic and 2nd-order
decoupled sliding-window configurations — the two that between them exercise both
estimator branches, the window buffers, and the `cfg`-struct type constraints.

Expected differences from today, both intentional and documented in
`Knowledge/block-library-signal-flow.md`: none in `mfc_siso_core`; one sample of
anti-windup latency in a loop assembled from stage blocks.
