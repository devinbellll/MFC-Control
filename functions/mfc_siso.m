classdef mfc_siso
%MFC_SISO Model-Free Control SISO controller: configuration, state and stages.
%
%   A namespace of static methods -- there is no object to construct. This
%   one file holds the whole controller lifecycle and every pipeline stage,
%   so the guts of the method can be read and edited in one place:
%
%     LIFECYCLE
%       cfg   = mfc_siso.config('Name', Value, ...)   build + validate a variant
%       state = mfc_siso.init(cfg)                    zero the controller state
%       [out, state] = mfc_siso.step(sp, y, t, u_prev, alpha, cfg, state)
%
%     PIPELINE STAGES (each is also wrapped by a block in blocks/)
%       [sp_filt, dot_sp, ddot_sp] = mfc_siso.ref_traj(...)   1. input smoother
%       (F-hat estimation dispatches to the mfc_fhat_* functions)  2. estimator
%       [fb, int_err, dot_err]     = mfc_siso.feedback(...)   3. feedback law
%       u_raw                      = mfc_siso.command(...)    4. command
%       [u, frozen]                = mfc_siso.limit(...)      5. filter + sat
%
%     HELPERS
%       kernel = mfc_siso.window_kernel(order, window_samples, Ts)
%       m      = mfc_siso.poly_moment(coeffs, j, a, b)
%
%   The stage methods take loose scalars rather than the state struct, so a
%   stage block that owns its own DiscreteState can call them directly
%   without assembling a full controller state. mfc_siso.step is then a
%   short pipeline over the same five calls -- the single source of truth
%   shared by the all-in-one block and the individual stage blocks.
%
%   The numerics live in plain functions alongside this file:
%   mfc_iir_smoother, mfc_fhat_algebraic_first_order,
%   mfc_fhat_algebraic_second_order, mfc_fhat_sliding_window.
%
%   Theory and derivations: see Knowledge/.
%
%   See also MFC_SISO_CORE, MFC_IIR_SMOOTHER, MFC_FHAT_ALGEBRAIC_FIRST_ORDER,
%   MFC_FHAT_ALGEBRAIC_SECOND_ORDER, MFC_FHAT_SLIDING_WINDOW.

methods (Static)

% =====================================================================
%  LIFECYCLE
% =====================================================================

function cfg = config(varargin)
%MFC_SISO.CONFIG Build and validate a controller configuration.
%
%   cfg = mfc_siso.config('Name', Value, ...)
%
%   Applies defaults, validates the variant selection, and precomputes the
%   sliding-window quadrature kernel. This is the single place where a
%   variant is fully specified.
%
%   Variant selection (2 x 2 x 2 grid, minus the undefined corner)
%     'model_order'  : 1 | 2 (default 2)
%         1: ultra-local model  dot_y = F + alpha*u
%         2: ultra-local model ddot_y = F + alpha*u
%     'structure'    : 'coupled' | 'decoupled' (default 'coupled')
%         coupled  : the F estimator is driven by the tracking ERROR. At
%                    2nd order the closed-loop polynomial s^2 + Kd*s + Kp is
%                    folded into the estimate; at 1st order the closed-loop
%                    pole s + Kp is folded (Kd has no derivative room to
%                    fold at 1st order, so it is unused there).
%         decoupled: the F estimator is driven by the pure MEASUREMENT
%                    (true-plant F); feedback is an explicit iP/iPD(I) law,
%                    at either model order.
%     'estimator'    : 'algebraic' | 'sliding_window' (default 'algebraic')
%         algebraic     : growing-window operational-calculus recursion.
%         sliding_window: fixed-length FIR window integral, exact taps
%                         (decoupled only -- coupled+sliding_window errors).
%
%   Channel count
%     'n'            : number of channels (default 1). n = 1 is the scalar
%         SISO controller and is the DEFAULT PATH in every sense -- the
%         command is a division, the anti-windup handshake is a scalar
%         test. n > 1 makes every signal n-by-1, alpha a square n-by-n
%         matrix inverted by MFC_SISO.COMMAND_MIMO, the gains scalars or
%         n-by-n matrices, and the saturation freeze per-channel. Nothing
%         else changes: the estimator kernels are vector-safe already,
%         because their numerator is per-element while their denominator
%         (t, t^2, or the FIR taps) is scalar and shared.
%
%   Tuning
%     'Ts'                : sample time [s] (default 0.01)
%     'alpha'             : ultra-local model input gain (default 1;
%                           square n-by-n when n > 1)
%     'Kp', 'Kd', 'Ki'    : feedback gains (defaults 25, 10, 0).
%         2nd order: characteristic polynomial s^2 + Kd*s + Kp
%                    (double pole at -p  <=>  Kd = 2p, Kp = p^2).
%         1st order: single pole at -Kp.
%         Kp is folded into F_hat whenever coupled (either order); explicit
%         whenever decoupled (either order). Kd is folded into F_hat when
%         coupled at 2nd order, unused when coupled at 1st order (no
%         derivative room to fold), and explicit whenever decoupled.
%         Ki is always applied explicitly (0 disables integral action).
%     'ref_filter_window' : reference trajectory filter memory [samples]
%                           (default 10; see MFC_IIR_SMOOTHER)
%     'est_filter_window' : estimator memory [samples] (default 10)
%                           algebraic: num/den IIR smoother memory;
%                           sliding_window: window length (rounded to even)
%     'est_hold_time'     : algebraic estimator held at 0 for t <= this [s]
%                           (default 0.1; ignored by sliding_window, which
%                           holds until its window fills)
%     'command_filter'    : output EMA constant, u = (raw + (c-1)*u_prev)/c;
%                           1 disables (default 1)
%     'use_ref_filter'    : true/false (default true). When false the raw
%                           setpoint passes through; feedforward derivatives
%                           are then finite differences of the raw setpoint.
%     'use_control_sat'   : true/false (default false); clamp to [u_min, u_max]
%                           with integrator freeze anti-windup
%     'u_min', 'u_max'    : saturation limits (defaults -600, 600)
%
%   Output: cfg struct with the fields above (structure/estimator resolved
%   to logicals cfg.coupled / cfg.algebraic) plus cfg.kernel.

    p = struct( ...
        'n',                 1, ...
        'model_order',       2, ...
        'structure',         'coupled', ...
        'estimator',         'algebraic', ...
        'Ts',                0.01, ...
        'alpha',             1, ...
        'Kp',                25, ...
        'Kd',                10, ...
        'Ki',                0, ...
        'ref_filter_window', 10, ...
        'est_filter_window', 10, ...
        'est_hold_time',     0.1, ...
        'command_filter',    1, ...
        'use_ref_filter',    true, ...
        'use_control_sat',   false, ...
        'u_min',             -600, ...
        'u_max',             600);

    assert(mod(numel(varargin), 2) == 0, 'mfc_siso.config: name-value pairs expected.');
    for i = 1:2:numel(varargin)
        name = varargin{i};
        assert(isfield(p, name), 'mfc_siso.config: unknown option ''%s''.', name);
        p.(name) = varargin{i+1};
    end

    % --- validate the variant selection ---------------------------------
    assert(p.model_order == 1 || p.model_order == 2, ...
        'mfc_siso.config: model_order must be 1 or 2.');
    assert(any(strcmp(p.structure, {'coupled', 'decoupled'})), ...
        'mfc_siso.config: structure must be ''coupled'' or ''decoupled''.');
    assert(any(strcmp(p.estimator, {'algebraic', 'sliding_window'})), ...
        'mfc_siso.config: estimator must be ''algebraic'' or ''sliding_window''.');

    coupled   = strcmp(p.structure, 'coupled');
    algebraic = strcmp(p.estimator, 'algebraic');

    assert(~(coupled && ~algebraic), ...
        ['mfc_siso.config: the coupled (error-driven, pole-folded) structure is ', ...
         'only defined for the algebraic estimator. Use structure=''decoupled'' ', ...
         'with estimator=''sliding_window''.']);

    % --- validate tuning -------------------------------------------------
    validateattributes(p.n, {'numeric'}, {'scalar', 'integer', 'positive'}, '', 'n');
    if p.n > 1
        assert(isequal(size(p.alpha), [p.n p.n]), ...
            'mfc_siso.config: with n = %d, alpha must be %d-by-%d.', p.n, p.n, p.n);
    end
    validateattributes(p.Ts, {'numeric'}, {'scalar', 'positive'}, '', 'Ts');
    validateattributes(p.command_filter, {'numeric'}, {'scalar', '>=', 1}, '', 'command_filter');
    validateattributes(p.ref_filter_window, {'numeric'}, {'scalar', 'nonnegative'}, '', 'ref_filter_window');
    validateattributes(p.est_filter_window, {'numeric'}, {'scalar', 'positive'}, '', 'est_filter_window');
    % all(): the limits may be n-by-1 when n > 1 (per-channel actuator range)
    assert(all(p.u_max > p.u_min), 'mfc_siso.config: u_max must exceed u_min.');

    % --- precompute the sliding-window kernel ----------------------------
    % Computed for EVERY variant: under code generation both estimator
    % branches of step() are compiled (cfg.algebraic is a run-time struct
    % field there), so cfg.kernel must always be a struct with one concrete
    % field set. Algebraic variants never evaluate it; their window value is
    % sanitized so any positive smoother memory remains accepted.
    if algebraic
        win = max(2, round(p.est_filter_window));
    else
        win = p.est_filter_window;
    end
    kernel = mfc_siso.window_kernel(p.model_order, win, p.Ts);

    % Assemble cfg in one pass: MATLAB Coder forbids adding struct fields
    % after the struct has been read, so every field must exist before use.
    cfg           = p;
    cfg.coupled   = coupled;
    cfg.algebraic = algebraic;
    cfg.kernel    = kernel;
end


function state = init(cfg)
%MFC_SISO.INIT Zero-initialize the controller state.
%
%   state = mfc_siso.init(cfg)
%
%   Window buffers are sized from the precomputed quadrature kernel for
%   sliding-window variants and collapse to scalar placeholders for
%   algebraic variants. Every field is sized from cfg.n, so the same state
%   struct serves the scalar controller (n = 1) and the vector one; the
%   denominator history stays SCALAR at any n, because the estimators share
%   one integration window across the channels.
%
%   State fields
%     .sp_filt_km1/2   reference trajectory filter history
%     .z_km1, .z_km2   algebraic estimator drive-signal history
%     .num_filt_km1/2  algebraic estimator smoothed numerator history
%     .den_filt_km1/2  algebraic estimator smoothed denominator history
%     .y_buf, .u_buf   sliding-window buffers [(win+1) x cfg.n], newest
%                      last, one COLUMN per channel
%     .err_km1         previous tracking error (derivative / integral)
%     .int_err         trapezoidal error integral
%     .u_km1           previously issued command (internal feedback path)

    if cfg.algebraic
        n_buf = 1;                             % unused placeholder
    else
        n_buf = cfg.kernel.n_intervals + 1;
    end

    nc = cfg.n;
    state = struct( ...
        'sp_filt_km1',  zeros(nc, 1), ...
        'sp_filt_km2',  zeros(nc, 1), ...
        'z_km1',        zeros(nc, 1), ...
        'z_km2',        zeros(nc, 1), ...
        'num_filt_km1', zeros(nc, 1), ...
        'num_filt_km2', zeros(nc, 1), ...
        'den_filt_km1', 0, ...            % t / t^2: scalar, shared by all channels
        'den_filt_km2', 0, ...
        'y_buf',        zeros(n_buf, nc), ...
        'u_buf',        zeros(n_buf, nc), ...
        'err_km1',      zeros(nc, 1), ...
        'int_err',      zeros(nc, 1), ...
        'u_km1',        zeros(nc, 1));
end


function [out, state] = step(setpoint, measure, t, u_prev, alpha, cfg, state)
%MFC_SISO.STEP One sample of the MFC SISO law (all variants).
%
%   [out, state] = mfc_siso.step(setpoint, measure, t, u_prev, alpha, cfg, state)
%
%   The complete, Simulink-independent controller: reference trajectory
%   filtering, F-hat estimation (dispatched to the variant selected in cfg),
%   feedback, command generation, output EMA filtering and saturation with
%   anti-windup. Each stage is one call to the corresponding stage method
%   below, so the all-in-one block and the individual stage blocks execute
%   the same code.
%
%   WIDTH. With cfg.n > 1 every signal is n-by-1, alpha is a square n-by-n
%   matrix and the command becomes a linear solve (MFC_SISO.COMMAND_MIMO);
%   the saturation freeze is then per-channel. The estimator dispatch and
%   the feedback law are unchanged -- the kernels are vector-safe, and the
%   gains may be scalars or n-by-n matrices. n = 1 executes exactly the
%   scalar code it always did.
%
%   Variant dispatch (cfg from mfc_siso.config):
%
%                        | algebraic (growing window) | sliding window
%     -------------------+-----------------------------+------------------
%     2nd order coupled  | error-driven, poles folded  |   (undefined)
%     2nd order decoupled| measurement-driven + iPD(I) | same, FIR window
%     1st order coupled  | error-driven, pole folded   |   (undefined)
%     1st order decoupled| measurement-driven + iP(I)  | same, FIR window
%
%   Command law (the fold vs explicit split is decided by coupled/decoupled
%   alone -- model order only selects the feedforward derivative order):
%     2nd order: u_raw = ( -F_hat + ddot_sp - fb ) / alpha
%     1st order: u_raw = ( -F_hat + dot_sp  - fb ) / alpha
%     fb = Kd*dot_err + Kp*err + Ki*int_err   (decoupled, either order)
%     fb =                       Ki*int_err   (coupled, either order: Kp
%                                              folded always, Kd folded only
%                                              at 2nd order)
%
%   Inputs
%     setpoint : reference input
%     measure  : plant measurement y
%     t        : current time [s] (growing-window estimators and hold logic)
%     u_prev   : command applied to the plant over the LAST sample. Pass
%                state.u_km1 for the standard internal feedback path, or the
%                externally measured/clamped command when another source
%                modifies what this controller issued.
%     alpha    : ultra-local model input gain for THIS sample (live value;
%                pass cfg.alpha if it does not vary)
%     cfg      : configuration from mfc_siso.config
%     state    : controller state from mfc_siso.init (updated in place)
%
%   Output struct 'out'
%     .u  .u_raw  .F_hat  .sp_filt  .err  .dot_sp  .ddot_sp
%     .est_valid    true once the estimator is past its startup hold
%     .est_dbg      estimator debug struct (see the mfc_fhat_* functions)

    Ts = cfg.Ts;

    % --- 1) Reference trajectory filter (input smoother) -----------------
    [sp_filt, dot_sp, ddot_sp] = mfc_siso.ref_traj( ...
        setpoint, state.sp_filt_km1, state.sp_filt_km2, Ts, ...
        cfg.ref_filter_window, cfg.use_ref_filter);

    err = measure - sp_filt;

    % --- 2) F-hat estimation, dispatched to the selected variant ---------
    if cfg.algebraic
        if cfg.coupled
            z_drive = err;              % error-driven (reference dynamics absorbed)
        else
            z_drive = measure;          % measurement-driven (pure-plant F)
        end
        if cfg.model_order == 2
            if cfg.coupled              % fold s^2 + Kd*s + Kp into the estimate
                a_fold = -cfg.Kd;  b_fold = -cfg.Kp;
            else
                a_fold = 0;        b_fold = 0;
            end
            [F_hat, state, est_dbg] = mfc_fhat_algebraic_second_order( ...
                z_drive, u_prev, alpha, t, Ts, cfg.est_filter_window, ...
                cfg.est_hold_time, a_fold, b_fold, state);
        else
            if cfg.coupled              % fold s + Kp into the estimate
                b_fold = -cfg.Kp;
            else
                b_fold = 0;
            end
            [F_hat, state, est_dbg] = mfc_fhat_algebraic_first_order( ...
                z_drive, u_prev, alpha, t, Ts, cfg.est_filter_window, ...
                cfg.est_hold_time, b_fold, state);
        end
    else                                % sliding window: decoupled only
        [F_hat, state, est_dbg] = mfc_fhat_sliding_window( ...
            measure, u_prev, alpha, t, cfg.kernel, state);
    end

    % --- 3) Explicit feedback --------------------------------------------
    [feedback, int_err] = mfc_siso.feedback( ...
        err, state.err_km1, state.int_err, Ts, ...
        cfg.Kp, cfg.Kd, cfg.Ki, cfg.coupled);

    % --- 4) Command generation -------------------------------------------
    if cfg.model_order == 2
        ff = ddot_sp;
    else
        ff = dot_sp;
    end
    if cfg.n > 1
        u_raw = mfc_siso.command_mimo(F_hat, ff, feedback, alpha);   % linear solve
    else
        u_raw = mfc_siso.command(F_hat, ff, feedback, alpha);        % division
    end

    % --- 5) Output EMA filter, then saturation with anti-windup ----------
    [u, frozen] = mfc_siso.limit(u_raw, u_prev, cfg.command_filter, ...
                                 cfg.use_control_sat, cfg.u_min, cfg.u_max);
    if cfg.n > 1
        % Per-channel freeze: one saturated actuator must not stop the
        % other channels' integrators. The scalar branch below is kept
        % literally as it was, so the SISO path is untouched.
        int_err(frozen) = state.int_err(frozen);
    elseif frozen
        int_err = state.int_err;        % freeze integrator (anti-windup)
    end

    % --- 6) Advance the shared controller state --------------------------
    state.sp_filt_km2 = state.sp_filt_km1;
    state.sp_filt_km1 = sp_filt;
    state.err_km1     = err;
    state.int_err     = int_err;
    state.u_km1       = u;

    out = struct( ...
        'u',         u, ...
        'u_raw',     u_raw, ...
        'F_hat',     F_hat, ...
        'sp_filt',   sp_filt, ...
        'err',       err, ...
        'dot_sp',    dot_sp, ...
        'ddot_sp',   ddot_sp, ...
        'est_valid', logical(est_dbg.valid), ...
        'est_dbg',   est_dbg);
end


% =====================================================================
%  PIPELINE STAGES
%  Loose scalars in, loose scalars out -- no state struct -- so a stage
%  block owning its own DiscreteState can call these directly.
% =====================================================================

function [sp_filt, dot_sp, ddot_sp] = ref_traj(setpoint, sp_km1, sp_km2, Ts, window, use_filter)
%MFC_SISO.REF_TRAJ Stage 1: input smoother and feedforward derivatives.
%
%   [sp_filt, dot_sp, ddot_sp] = mfc_siso.ref_traj(setpoint, sp_km1, sp_km2, ...
%                                                  Ts, window, use_filter)
%
%   Shapes the raw setpoint into a twice-differentiable reference trajectory
%   with MFC_IIR_SMOOTHER and returns its first two derivatives as backward
%   finite differences of the stored trajectory history.
%
%   With use_filter = false the raw setpoint passes through unchanged; the
%   derivatives are then finite differences of the RAW setpoint, which for a
%   step input are impulsive. That is exactly why the smoother exists: the
%   command law feeds ddot_sp (2nd order) straight through 1/alpha.
%
%   Inputs
%     setpoint   : raw reference at step k
%     sp_km1/2   : filtered trajectory history (the caller's state)
%     Ts         : sample time [s]
%     window     : smoother memory [samples]; 0 is exact pass-through
%     use_filter : false bypasses the smoother
%
%   Wrapped by blocks/mfc_smoother_block.

    if use_filter
        sp_filt = mfc_iir_smoother(setpoint, sp_km1, sp_km2, window);
    else
        sp_filt = setpoint;
    end
    dot_sp  = (sp_filt - sp_km1) / Ts;
    ddot_sp = (sp_filt - 2*sp_km1 + sp_km2) / Ts^2;
end


function [fb, int_err, dot_err] = feedback(err, err_km1, int_err_km1, Ts, Kp, Kd, Ki, coupled)
%MFC_SISO.FEEDBACK Stage 3: the explicit feedback law.
%
%   [fb, int_err, dot_err] = mfc_siso.feedback(err, err_km1, int_err_km1, ...
%                                              Ts, Kp, Kd, Ki, coupled)
%
%   Everything the estimator did NOT already absorb:
%
%     coupled   : fb = Ki*int_err
%                 Kp is folded into F_hat at either model order, Kd is
%                 folded at 2nd order and unused at 1st, so only the
%                 integral term remains explicit here.
%     decoupled : fb = Kd*dot_err + Kp*err + Ki*int_err
%                 the full iPD(I) law, at either model order.
%
%   The integral is trapezoidal and is ALWAYS returned advanced, whether or
%   not Ki is zero and whether or not the coupled branch uses it. The caller
%   is responsible for anti-windup: if MFC_SISO.LIMIT reports a frozen
%   sample, discard this int_err and keep int_err_km1.
%
%   dot_err is a raw backward difference and is the noise-sensitive term of
%   the whole controller -- it is returned even in the coupled case (where
%   fb ignores it) so a block can log or filter it.
%
%   No block wraps this: in an assembled Simulink loop the explicit law is
%   a stock Discrete PID Controller into mfc_command_block's fb input (and
%   with a coupled estimator, fb is simply Ground). It stays here because
%   mfc_siso.step -- and therefore mfc_siso_core -- still uses it.

    int_err = int_err_km1 + (err + err_km1)/2 * Ts;   % trapezoidal integral
    dot_err = (err - err_km1) / Ts;                   % noise-sensitive term

    if coupled
        fb = Ki*int_err;                              % Kp/Kd already in F_hat
    else
        fb = Kd*dot_err + Kp*err + Ki*int_err;        % explicit iPD(I)
    end
end


function u_raw = command(F_hat, ff, fb, alpha)
%MFC_SISO.COMMAND Stage 4: the ultra-local model inversion.
%
%   u_raw = mfc_siso.command(F_hat, ff, fb, alpha)
%
%   The one line that makes this model-free control: invert the ultra-local
%   model, cancelling the estimated lumped dynamics and injecting the
%   reference feedforward and the feedback.
%
%       u_raw = ( -F_hat + ff - fb ) / alpha
%
%   ff is the feedforward derivative of the reference trajectory: ddot_sp
%   for a 2nd-order ultra-local model, dot_sp for 1st order. That is the
%   ONLY place model order enters the command law -- the fold-vs-explicit
%   split is decided by coupled/decoupled, not by order.
%
%   alpha is the input gain of the ultra-local model. It is a design
%   parameter, not a plant identification: too small over-drives the
%   command, too large under-drives it, and the estimator absorbs the
%   mismatch into F_hat either way.
%
%   Wrapped by blocks/mfc_command_block.

    u_raw = (-F_hat + ff - fb) / alpha;
end


function u_raw = command_mimo(F_hat, ff, fb, alpha)
%MFC_SISO.COMMAND_MIMO Stage 4 with a MATRIX input gain.
%
%   u_raw = mfc_siso.command_mimo(F_hat, ff, fb, alpha)
%
%   Same inversion as MFC_SISO.COMMAND, but F_hat, ff, fb and u_raw are
%   n-by-1 vectors and alpha is a square n-by-n matrix, so the scalar
%   division becomes a linear solve:
%
%       alpha * u_raw = -F_hat + ff - fb
%
%   alpha must be invertible; it is a design parameter, not a plant
%   identification, so choose it well-conditioned -- the estimator absorbs
%   the mismatch into F_hat, but not the loss of rank.
%
%   Wrapped by blocks/mfc_command_mimo_block.

    u_raw = alpha \ (-F_hat + ff - fb);
end


function [u, frozen] = limit(u_raw, u_prev, command_filter, use_sat, u_min, u_max)
%MFC_SISO.LIMIT Stage 5: output EMA filter and saturation with anti-windup.
%
%   [u, frozen] = mfc_siso.limit(u_raw, u_prev, command_filter, use_sat, u_min, u_max)
%
%   Exponential moving average on the raw command,
%
%       u = ( u_raw + (c - 1)*u_prev ) / c,      c = command_filter
%
%   (c = 1 is an exact pass-through), then an optional clamp to
%   [u_min, u_max].
%
%   'frozen' is the anti-windup handshake: true when the clamp actually bit
%   this sample, meaning the caller must discard the advanced integral and
%   keep the previous one. This method does not own the integrator, so it
%   cannot freeze it itself -- it only reports.
%
%   No block wraps this: mfc_siso_core is the only user, and it is the only
%   place the anti-windup handshake exists at all. mfc_command_block offers
%   a plain clamp with no freeze (nothing to freeze -- it owns no
%   integrator), so if you need integral action against a real actuator
%   limit, use mfc_siso_core or a PID block with its own anti-windup.
%   See Knowledge/block-library-signal-flow.md.

    u      = (u_raw + (command_filter - 1)*u_prev) / command_filter;
    frozen = false;
    if use_sat
        u_sat = min(u_max, max(u_min, u));
        frozen = (u_sat ~= u);
        u = u_sat;
    end
end


% =====================================================================
%  HELPERS
% =====================================================================

function kernel = window_kernel(model_order, window_samples, Ts)
%MFC_SISO.WINDOW_KERNEL Precompute the sliding-window FIR taps.
%
%   kernel = mfc_siso.window_kernel(model_order, window_samples, Ts)
%
%   Builds the constant part of the sliding-window (non-algebraic) F
%   estimators used by MFC_FHAT_SLIDING_WINDOW: one fixed multiplier -- a
%   TAP -- for every stored sample of y and of u. Everything that does not
%   depend on the signals is computed once here, so the estimator itself is
%   a multiply-accumulate over the two windows.
%
%   model_order = 1  (ultra-local model dot_y = F + alpha*u, Eq. 11):
%       F = -(6/Tw^3) * int_0^Tw [ (Tw - 2*sigma) * y(sigma)
%                                  + alpha*sigma*(Tw - sigma) * u(sigma) ] dsigma
%
%   model_order = 2  (ultra-local model ddot_y = F + alpha*u, Eq. 16):
%       F = (60/Tw^5) * int_0^Tw [ (Tw^2 - 6*Tw*sigma + 6*sigma^2) * y(sigma)
%                                  - (alpha/2)*sigma^2*(Tw - sigma)^2 * u(sigma) ] dsigma
%
%   Both are returned in the common form used by MFC_FHAT_SLIDING_WINDOW:
%
%       F = tap_y' * y_buf  +  alpha * ( tap_u_unit' * u_buf )
%
%   (the prefactor, the sign and the 1/2 of the second-order input kernel
%   are folded into the taps; alpha is left out of tap_u_unit so it can stay
%   a live run-time input).
%
%   HOW THE TAPS ARE COMPUTED. The weighting kernels above are POLYNOMIALS
%   we wrote down ourselves, so there is no reason to approximate them --
%   the only genuine ignorance is what the signals did BETWEEN samples.
%   Each tap is therefore the exact integral of the kernel against that
%   sample's interpolation basis:
%
%     y  is modelled as piecewise LINEAR (a "tent" basis, since y is
%        sampled and nothing better is known), so
%            tap_y(i) = int K_y(sigma) * tent_i(sigma) dsigma
%        which is a cubic (1st order) or quartic (2nd order) integral over
%        the two half-intervals either side of node i -- closed form, power
%        rule, no quadrature rule involved.
%
%     u  is modelled as piecewise CONSTANT, which is not a model at all but
%        the truth: the command reaches the plant through a zero-order
%        hold. So tap_u_unit(j) is the exact integral of K_u over the one
%        interval that sample was held across, and the input term carries
%        NO quadrature error whatsoever.
%
%   This also fixes the alignment: u_buf(end) is u_prev, the command held
%   over the interval ENDING now, so it owns the last interval and the
%   oldest u sample (whose interval fell out of the window) gets a zero tap.
%
%   Composite Simpson was used here previously. It is a general-purpose
%   rule that approximates the whole product K*y, including the K we know
%   exactly, and it weights samples 1,4,2,4,...,1 -- a lumpiness that costs
%   both accuracy (0.24% vs 0.010% gain error on an 11-sample window; a
%   steady factor of ~24 at any size) and about 10% more noise for nothing.
%   It also forced window_samples to be rounded UP to even. None of that
%   applies any more: any window length >= 2 intervals is valid and means
%   exactly what it says.
%
%   Sanity properties the taps must have, all checked in tests:
%     sum(tap_y) == 0          blind to a constant y (no acceleration)
%     sum(sigma.*tap_y) == 0   blind to a ramp y (no acceleration either)
%     tap_y' * (sigma.^2/2) == 1 (2nd order)   unity gain on acceleration
%
%   Output: kernel struct with fields
%     .model_order      1 or 2
%     .n_intervals      interval count (buffers hold n_intervals+1 samples)
%     .Tw               window length = n_intervals*Ts [s]
%     .Ts               sample time [s]
%     .sigma            [(n+1)x1] window-local time nodes, 0..Tw
%     .tap_y            [(n+1)x1] measurement taps, prefactor folded in
%     .tap_u_unit       [(n+1)x1] input taps for alpha = 1, prefactor folded in

    validateattributes(model_order, {'numeric'}, {'scalar'});
    assert(model_order == 1 || model_order == 2, ...
           'mfc_siso.window_kernel: model_order must be 1 or 2.');
    validateattributes(window_samples, {'numeric'}, {'scalar', 'integer', '>=', 2});
    validateattributes(Ts, {'numeric'}, {'scalar', 'positive'});

    n     = window_samples;
    Tw    = n * Ts;
    sigma = (0:n).' * Ts;                              % window-local time, 0..Tw

    % Kernel polynomials in ascending powers of sigma, degree 4 either way
    % so the two orders share one code path.
    if model_order == 1
        cy = [Tw, -2, 0, 0, 0];                        % Tw - 2*sigma
        cu = [0, Tw, -1, 0, 0];                        % sigma*(Tw - sigma)
        prefactor = -6 / Tw^3;
    else
        cy = [Tw^2, -6*Tw, 6, 0, 0];                   % Tw^2 - 6*Tw*s + 6*s^2
        cu = [0, 0, -0.5*Tw^2, Tw, -0.5];              % -0.5*s^2*(Tw - s)^2
        prefactor = 60 / Tw^5;
    end

    tap_y      = zeros(n + 1, 1);
    tap_u_unit = zeros(n + 1, 1);

    for i = 1:n+1
        s_i = sigma(i);
        acc = 0;

        % Rising half of the tent, (sigma - a)/Ts on [s_i - Ts, s_i]
        if i > 1
            a   = s_i - Ts;
            acc = acc + (mfc_siso.poly_moment(cy, 1, a, s_i) ...
                         - a*mfc_siso.poly_moment(cy, 0, a, s_i)) / Ts;
        end

        % Falling half of the tent, (b - sigma)/Ts on [s_i, s_i + Ts]
        if i < n+1
            b   = s_i + Ts;
            acc = acc + (b*mfc_siso.poly_moment(cy, 0, s_i, b) ...
                         - mfc_siso.poly_moment(cy, 1, s_i, b)) / Ts;
        end

        tap_y(i) = prefactor * acc;

        % Zero-order hold: sample i was held over the interval ENDING at
        % s_i. The oldest sample's interval is outside the window, so its
        % tap stays 0.
        if i > 1
            tap_u_unit(i) = prefactor * mfc_siso.poly_moment(cu, 0, sigma(i-1), s_i);
        end
    end

    kernel = struct( ...
        'model_order',  model_order, ...
        'n_intervals',  n, ...
        'Tw',           Tw, ...
        'Ts',           Ts, ...
        'sigma',        sigma, ...
        'tap_y',        tap_y, ...
        'tap_u_unit',   tap_u_unit);
end


function m = poly_moment(c, j, a, b)
%MFC_SISO.POLY_MOMENT Exact integral of sigma^j * p(sigma) over [a, b].
%
%   m = mfc_siso.poly_moment(c, j, a, b)
%
%   p is given by its coefficients in ASCENDING powers,
%   p(sigma) = c(1) + c(2)*sigma + ... , and the integral is evaluated by
%   the power rule, term by term:
%
%       int_a^b sigma^j * p(sigma) dsigma
%           = sum_k c(k+1) * ( b^(k+j+1) - a^(k+j+1) ) / (k+j+1)
%
%   Exact, not approximate -- this is the whole reason MFC_SISO.WINDOW_KERNEL
%   needs no quadrature rule. Only j = 0 and j = 1 are used there.

    m = 0;
    for k = 0:numel(c)-1
        p = k + j + 1;
        m = m + c(k+1) * (b^p - a^p) / p;
    end
end

end
end
