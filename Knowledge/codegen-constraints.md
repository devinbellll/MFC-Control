# Code generation constraints

Several things in this codebase look over-careful until you try to run
MATLAB Coder on them. This note records **why** each one is written that way, so
nobody "simplifies" it back into a build failure.

All of these were established in commit `efabee4`
(*"make MFC SISO MATLAB System block compatible with code generation"*).

## The rule behind most of them

Under code generation, **both branches of every `if` are compiled**, even ones
that can never execute for a given configuration. `cfg.algebraic` is a run-time
struct field inside `mfc_siso.step`, not a compile-time constant, so the
sliding-window branch is compiled even for a purely algebraic controller — and
it must be *type-valid*, not merely unreachable.

Nearly every oddity below follows from that.

## 1. `cfg.kernel` always exists

```matlab
% mfc_siso.config
if algebraic
    win = max(2, round(p.est_filter_window));
else
    win = p.est_filter_window;
end
kernel = mfc_siso.window_kernel(p.model_order, win, p.Ts);
```

The quadrature kernel is built for **every** variant, including algebraic ones
that never evaluate it, so `cfg` has a single concrete type. The algebraic
window value is sanitized (`max(2, round(...))`) only so that any positive
smoother memory — including a non-integer one — remains an acceptable input to
`window_kernel`'s integer validation.

## 2. `cfg` is assembled in one pass

```matlab
cfg           = p;
cfg.coupled   = coupled;
cfg.algebraic = algebraic;
cfg.kernel    = kernel;
```

MATLAB Coder forbids adding struct fields after a struct has been read. Every
field must exist before first use — hence the single assembly block at the end
of `mfc_siso.config` rather than fields accreting as they are computed.

## 3. All three estimators return the same `dbg` struct

```matlab
% mfc_fhat_algebraic_*: 'integral' is zero-filled
dbg = struct('num_raw', num_raw, 'den_raw', den_raw, ...
             'num_filt', num_filt, 'den_filt', den_filt, ...
             'integral', 0, 'valid', valid);

% mfc_fhat_sliding_window: num/den are zero-filled
dbg = struct('num_raw', 0, 'den_raw', 0, 'num_filt', 0, 'den_filt', 0, ...
             'integral', integral, 'valid', valid);
```

The dispatch in `mfc_siso.step` assigns `est_dbg` from either branch, so both
must produce **identically shaped** structs. The zero-filled fields are not
padding for tidiness; they are what makes the dispatch compile.

## 4. Placeholder window buffers, and `circshift`

Algebraic variants still carry `y_buf`/`u_buf`, sized $1\times1$:

```matlab
% mfc_siso.init
if cfg.algebraic
    n_buf = 1;                             % unused placeholder
else
    n_buf = cfg.kernel.n_intervals + 1;
end
```

And the sliding-window estimator shifts them with `circshift`, not slicing:

```matlab
state.y_buf      = circshift(state.y_buf, -1);
state.y_buf(end) = y;
```

`[buf(2:end); new]` gives an identical result for a real buffer but is invalid
for a $1\times1$ one — and that dead $1\times1$ code path is still compiled for
algebraic variants. `circshift` is valid for both.

## 5. `resetImpl` types the discrete states, so buffer length must be constant

```matlab
% mfc_siso_core.resetImpl -- NOT mfc_siso.init
n_buf = bufferLength(obj);
obj.y_buf = zeros(n_buf, 1);
```

`mfc_siso_core` deliberately zeroes its state inline instead of calling
`mfc_siso.init`. Code generation infers the discrete states' types and sizes
from these assignments, so the buffer length must be a **compile-time constant**
— and `mfc_siso.init` sizes them from a run-time `cfg` value.

`bufferLength` is therefore built only from `Nontunable` properties, and is
shared with `getDiscreteStateSpecificationImpl` so the two cannot diverge:

```matlab
function n_buf = bufferLength(obj)
    if strncmp(obj.estimator_type, 'Algebraic', 9)
        n_buf = 1;
    else
        n = obj.est_filter_window;
        n = n + mod(n, 2);
        n_buf = n + 1;
    end
end
```

`mfc_fhat_window_block` follows the same pattern.

## 6. `est_filter_window` is Nontunable

It sizes the window buffers and the quadrature kernel. Making it tunable would
mean re-sizing a discrete state at run time, which is not possible under
codegen. `Ts` is Nontunable for the same class of reason — it fixes the block's
sample time.

Gains (`Kp`, `Kd`, `Ki`, `alpha`, `command_filter`, `u_min`, `u_max`) *are*
tunable, and `mfc_siso_core.stepImpl` copies them into its cached config every
sample so live changes take effect:

```matlab
c = obj.cfg;
c.alpha = obj.alpha;
c.Kp    = obj.Kp;
...
```

## 7. Fixed discrete rate, never inherited

```matlab
function sts = getSampleTimeImpl(obj)
    sts = createSampleTime(obj, 'Type', 'Discrete', 'SampleTime', obj.Ts);
end
```

Every stateful block does this. The backward differences, trapezoidal integral,
EMA and window shifts all assume they advance **exactly once per `Ts`**.
Inheriting a rate would silently invalidate every recursion.

## The `mfc_siso` classdef namespace

`config`, `init`, `step`, the pipeline stages and `window_kernel` are static
methods of a `classdef` with **no properties and no instances**. MATLAB allows
one public function per `.m` file, so this is the only way to get several entry
points into one file.

Coder supports static-method calls on a codegen-compatible classdef, and Octave
11 handles them too (verified: `varargin`, dynamic field assignment,
`validateattributes` and static-to-static calls all work). If a future toolchain
objects, the mechanical fallback is a `+mfc/` package folder — same split, plain
function files, no classdef.

## Verifying

```matlab
>> cfg = coder.config('lib');
>> codegen -config cfg -report mfc_siso_core -args {0,0,0}
```

Check at minimum **2nd-order coupled algebraic** and **2nd-order decoupled
sliding window**. Between them they exercise both estimator branches, the window
buffers, and the `cfg` struct-type constraints above.

## See also

- [[block-library-signal-flow]] — the blocks these constraints shape
- [[estimator-sliding-window]] — where the buffers and kernel come from
