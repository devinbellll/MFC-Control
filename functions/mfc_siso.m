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
%         sliding_window: fixed-length Simpson-quadrature integral
%                         (decoupled only -- coupled+sliding_window errors).
%
%   Tuning
%     'Ts'                : sample time [s] (default 0.01)
%     'alpha'             : ultra-local model input gain (default 1)
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
    validateattributes(p.Ts, {'numeric'}, {'scalar', 'positive'}, '', 'Ts');
    validateattributes(p.command_filter, {'numeric'}, {'scalar', '>=', 1}, '', 'command_filter');
    validateattributes(p.ref_filter_window, {'numeric'}, {'scalar', 'nonnegative'}, '', 'ref_filter_window');
    validateattributes(p.est_filter_window, {'numeric'}, {'scalar', 'positive'}, '', 'est_filter_window');
    assert(p.u_max > p.u_min, 'mfc_siso.config: u_max must exceed u_min.');

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
%   algebraic variants.
%
%   State fields
%     .sp_filt_km1/2   reference trajectory filter history
%     .z_km1, .z_km2   algebraic estimator drive-signal history
%     .num_filt_km1/2  algebraic estimator smoothed numerator history
%     .den_filt_km1/2  algebraic estimator smoothed denominator history
%     .y_buf, .u_buf   sliding-window buffers [(n+1)x1], newest last
%     .err_km1         previous tracking error (derivative / integral)
%     .int_err         trapezoidal error integral
%     .u_km1           previously issued command (internal feedback path)

    if cfg.algebraic
        n_buf = 1;                             % unused placeholder
    else
        n_buf = cfg.kernel.n_intervals + 1;
    end

    state = struct( ...
        'sp_filt_km1',  0, ...
        'sp_filt_km2',  0, ...
        'z_km1',        0, ...
        'z_km2',        0, ...
        'num_filt_km1', 0, ...
        'num_filt_km2', 0, ...
        'den_filt_km1', 0, ...
        'den_filt_km2', 0, ...
        'y_buf',        zeros(n_buf, 1), ...
        'u_buf',        zeros(n_buf, 1), ...
        'err_km1',      0, ...
        'int_err',      0, ...
        'u_km1',        0);
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
%   Variant dispatch (cfg from mfc_siso.config):
%
%                        | algebraic (growing window) | sliding window
%     -------------------+-----------------------------+------------------
%     2nd order coupled  | error-driven, poles folded  |   (undefined)
%     2nd order decoupled| measurement-driven + iPD(I) | same, Simpson
%     1st order coupled  | error-driven, pole folded   |   (undefined)
%     1st order decoupled| measurement-driven + iP(I)  | same, Simpson
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
    u_raw = mfc_siso.command(F_hat, ff, feedback, alpha);

    % --- 5) Output EMA filter, then saturation with anti-windup ----------
    [u, frozen] = mfc_siso.limit(u_raw, u_prev, cfg.command_filter, ...
                                 cfg.use_control_sat, cfg.u_min, cfg.u_max);
    if frozen
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
%   Wrapped by blocks/mfc_feedback_block.

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
%   Wrapped by blocks/mfc_command_filter_block, whose 'sat' output is this
%   flag. NOTE: in a loop assembled from separate blocks that flag must pass
%   through a unit delay to reach the feedback block without forming an
%   algebraic loop, so a composed loop freezes ONE SAMPLE LATER than
%   mfc_siso.step does. See Knowledge/block-library-signal-flow.md.

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
%MFC_SISO.WINDOW_KERNEL Precompute the sliding-window quadrature kernel.
%
%   kernel = mfc_siso.window_kernel(model_order, window_samples, Ts)
%
%   Builds the constant part of the sliding-window (non-algebraic) F
%   estimators used by MFC_FHAT_SLIDING_WINDOW. The estimators evaluate a
%   weighted integral over a window of length Tw; everything that does not
%   depend on the signals (time nodes, measurement kernel, unit input
%   kernel, Simpson weights, prefactor) is computed once here.
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
%       F = prefactor * (Ts/3) * sum( simpson_weights .* ...
%             ( y_kernel .* y_buf + alpha * u_kernel_unit .* u_buf ) )
%
%   (the sign and the 1/2 of the second-order input kernel are folded into
%   u_kernel_unit, so alpha can stay a live run-time input).
%
%   QUADRATURE: composite Simpson with an EVEN interval count
%   (window_samples is rounded up to even; realized window Tw =
%   n_intervals*Ts). Trapezoidal integration is NOT usable for the
%   second-order kernel: the 60/Tw^5 prefactor amplifies its O(Ts^2/Tw^4)
%   leakage to ~60x error at practical sample times. Simpson is exact
%   through cubics, so the sliding-window estimate matches the algebraic
%   variants.
%
%   Output: kernel struct with fields
%     .model_order      1 or 2
%     .n_intervals      even interval count (buffers hold n_intervals+1 samples)
%     .Tw               realized window length = n_intervals*Ts [s]
%     .Ts               sample time [s]
%     .sigma            [(n+1)x1] window-local time nodes, 0..Tw
%     .y_kernel         [(n+1)x1] measurement weighting
%     .u_kernel_unit    [(n+1)x1] input weighting for alpha = 1
%     .simpson_weights  [(n+1)x1] composite-Simpson weights (1 4 2 ... 4 1)
%     .prefactor        -6/Tw^3 (1st order) or 60/Tw^5 (2nd order)

    validateattributes(model_order, {'numeric'}, {'scalar'});
    assert(model_order == 1 || model_order == 2, ...
           'mfc_siso.window_kernel: model_order must be 1 or 2.');
    validateattributes(window_samples, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(Ts, {'numeric'}, {'scalar', 'positive'});

    n     = window_samples + mod(window_samples, 2);   % force even for Simpson
    Tw    = n * Ts;
    sigma = (0:n).' * Ts;                              % window-local time, 0..Tw

    if model_order == 1
        y_kernel      = Tw - 2*sigma;
        u_kernel_unit = sigma .* (Tw - sigma);
        prefactor     = -6 / Tw^3;
    else
        y_kernel      = Tw^2 - 6*Tw*sigma + 6*sigma.^2;
        u_kernel_unit = -0.5 * sigma.^2 .* (Tw - sigma).^2;
        prefactor     = 60 / Tw^5;
    end

    simpson_weights            = ones(n + 1, 1);
    simpson_weights(2:2:end-1) = 4;
    simpson_weights(3:2:end-1) = 2;

    kernel = struct( ...
        'model_order',     model_order, ...
        'n_intervals',     n, ...
        'Tw',              Tw, ...
        'Ts',              Ts, ...
        'sigma',           sigma, ...
        'y_kernel',        y_kernel, ...
        'u_kernel_unit',   u_kernel_unit, ...
        'simpson_weights', simpson_weights, ...
        'prefactor',       prefactor);
end

end
end
