function [out, state] = mfc_siso_step(setpoint, measure, t, u_prev, alpha, cfg, state)
%MFC_SISO_STEP One sample of the Model-Free Control SISO law (all variants).
%
%   [out, state] = MFC_SISO_STEP(setpoint, measure, t, u_prev, alpha, cfg, state)
%
%   Executes one controller sample: reference trajectory filtering, F-hat
%   estimation (dispatched to the variant selected in cfg), command
%   generation, output EMA filtering and saturation with anti-windup. This
%   is the complete, Simulink-independent controller core; the block
%   MFC_SISO_CORE is a thin wrapper around it.
%
%   Variant dispatch (cfg from MFC_SISO_CONFIG):
%
%                        | algebraic (growing window) | sliding window
%     -------------------+-----------------------------+------------------
%     2nd order coupled  | error-driven, poles folded  |   (undefined)
%     2nd order decoupled| measurement-driven + iPD(I) | same, Simpson
%     1st order coupled  | error-driven + explicit iP  |   (undefined)
%     1st order decoupled| measurement-driven + iP(I)  | same, Simpson
%
%   Command law (feedback is explicit except what is folded into F_hat):
%     2nd order: u_raw = ( -F_hat + ddot_sp - fb ) / alpha
%                fb = Kd*dot_err + Kp*err + Ki*int_err   (decoupled)
%                fb =                       Ki*int_err   (coupled: Kp, Kd folded)
%     1st order: u_raw = ( -F_hat + dot_sp - fb ) / alpha
%                fb = Kp*err + Ki*int_err                (both structures)
%
%   Then u = EMA(u_raw) via cfg.command_filter, clamped to
%   [cfg.u_min, cfg.u_max] when cfg.use_control_sat (integrator frozen on
%   saturation).
%
%   Inputs
%     setpoint : reference input
%     measure  : plant measurement y
%     t        : current time [s] (growing-window estimators and hold logic)
%     u_prev   : command applied to the plant over the LAST sample. Pass
%                state.u_km1 for the standard internal feedback path, or the
%                externally measured/clamped command when another source
%                (e.g. an actuator model or a supervisor) modifies what this
%                controller issued. Used by the estimator, the EMA filter,
%                and the window buffers.
%     alpha    : ultra-local model input gain for THIS sample (live value;
%                pass cfg.alpha if it does not vary)
%     cfg      : configuration from MFC_SISO_CONFIG
%     state    : controller state from MFC_SISO_INIT (updated in place)
%
%   Output struct 'out'
%     .u            issued command (after EMA filter and saturation)
%     .u_raw        command before EMA filter and saturation
%     .F_hat        current F estimate
%     .sp_filt      filtered reference trajectory
%     .err          tracking error, measure - sp_filt
%     .dot_sp       first derivative of the filtered reference
%     .ddot_sp      second derivative of the filtered reference
%     .est_valid    true once the estimator is past its startup hold
%     .est_dbg      estimator debug struct (see the mfc_fhat_* functions)
%
%   See also MFC_SISO_CONFIG, MFC_SISO_INIT, MFC_SISO_CORE,
%   MFC_FHAT_ALGEBRAIC_FIRST_ORDER, MFC_FHAT_ALGEBRAIC_SECOND_ORDER,
%   MFC_FHAT_SLIDING_WINDOW.

Ts = cfg.Ts;

% 1) Reference trajectory filter (or raw pass-through). The feedforward
%    derivatives always come from the stored trajectory history, so with the
%    filter disabled they are finite differences of the raw setpoint.
if cfg.use_ref_filter
    sp_filt = mfc_iir_smoother(setpoint, state.sp_filt_km1, state.sp_filt_km2, ...
                               cfg.ref_filter_window);
else
    sp_filt = setpoint;
end
dot_sp  = (sp_filt - state.sp_filt_km1) / Ts;
ddot_sp = (sp_filt - 2*state.sp_filt_km1 + state.sp_filt_km2) / Ts^2;

err = measure - sp_filt;

% 2) F-hat estimation, dispatched to the selected variant.
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
        [F_hat, state, est_dbg] = mfc_fhat_algebraic_first_order( ...
            z_drive, u_prev, alpha, t, Ts, cfg.est_filter_window, ...
            cfg.est_hold_time, state);
    end
else                                % sliding window: decoupled only
    [F_hat, state, est_dbg] = mfc_fhat_sliding_window( ...
        measure, u_prev, alpha, t, cfg.kernel, state);
end

% 3) Explicit feedback and command generation.
int_err = state.int_err + (err + state.err_km1)/2 * Ts;   % trapezoidal integral

if cfg.model_order == 2
    if cfg.coupled
        feedback = cfg.Ki*int_err;              % Kp, Kd already folded in F_hat
    else
        dot_err  = (err - state.err_km1) / Ts;  % noise-sensitive term
        feedback = cfg.Kd*dot_err + cfg.Kp*err + cfg.Ki*int_err;
    end
    u_raw = (-F_hat + ddot_sp - feedback) / alpha;
else
    % 1st order: single closed-loop pole; Kp always explicit (no fold room).
    feedback = cfg.Kp*err + cfg.Ki*int_err;
    u_raw = (-F_hat + dot_sp - feedback) / alpha;
end

% 4) Output EMA filter, then saturation with integrator-freeze anti-windup.
u = (u_raw + (cfg.command_filter - 1)*u_prev) / cfg.command_filter;
if cfg.use_control_sat
    u_sat = min(cfg.u_max, max(cfg.u_min, u));
    if u_sat ~= u
        int_err = state.int_err;                % freeze integrator (anti-windup)
    end
    u = u_sat;
end

% 5) Advance the shared controller state.
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
