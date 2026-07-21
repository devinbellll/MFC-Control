function [F_hat, state, dbg] = mfc_fhat_algebraic_first_order( ...
    z, u_prev, alpha, t, Ts, filter_window, hold_time, b_fold, state)
%MFC_FHAT_ALGEBRAIC_FIRST_ORDER Algebraic (growing-window) F estimator, 1st order.
%
%   [F_hat, state, dbg] = MFC_FHAT_ALGEBRAIC_FIRST_ORDER( ...
%       z, u_prev, alpha, t, Ts, filter_window, hold_time, b_fold, state)
%
%   Estimates F in the first-order ultra-local model
%
%       dot_z = F + alpha * u + b_fold * z
%
%   using the algebraic (operational-calculus) method: the model is written
%   in the Laplace domain, differentiated once w.r.t. s to annihilate the
%   unknown initial condition, and mapped back to the time domain, where
%   multiplication by t replaces -d/ds. The resulting growing-window
%   expressions (time-weighted from the start of the run) are discretized
%   with backward differences:
%
%       num[k] = -z[k] + ( t*z[k] - (t-Ts)*z[k-1] ) / Ts - t*alpha*u_prev - t*b_fold*z[k]
%       den[k] = t
%
%   b_fold*z is a known-coefficient term added to the model on the same
%   footing as alpha*u (both order-0, known-coefficient terms), so it is
%   annihilated by the identical t^1 transform.
%
%   Both are smoothed with the shared second-order IIR (MFC_IIR_SMOOTHER,
%   memory = filter_window samples) before the division, and the estimate is
%   held at zero until t > hold_time so early transients cannot blow up the
%   near-zero denominator.
%
%   COUPLED vs DECOUPLED use (selected by the caller through z):
%     * decoupled: z = y_measured, b_fold = 0  -> F_hat is the pure-plant
%       lumped dynamics dot_y - alpha*u; stabilizing feedback must be added
%       explicitly in the command law.
%     * coupled:   z = tracking error, b_fold = -Kp  -> the closed-loop
%       characteristic polynomial s + Kp is folded into the estimate, so the
%       command law needs no explicit proportional feedback.
%
%   Inputs
%     z             : estimator drive signal at step k (measurement or error)
%     u_prev        : command actually applied over the last sample
%     alpha         : ultra-local model input gain
%     t             : current time [s] (from run/sim start; growing window)
%     Ts            : sample time [s]
%     filter_window : IIR smoother memory [samples]
%     hold_time     : F_hat is forced to 0 while t <= hold_time [s]
%     b_fold        : folded z coefficient (coupled: -Kp, decoupled: 0)
%     state         : struct, fields used/updated here:
%                       .z_km1, .z_km2            signal history
%                       .num_filt_km1/2           smoothed numerator history
%                       .den_filt_km1/2           smoothed denominator history
%
%   Outputs
%     F_hat : estimate of F (0 while invalid)
%     state : updated state struct
%     dbg   : debug struct: num_raw, den_raw, num_filt, den_filt, valid
%
%   See also MFC_FHAT_ALGEBRAIC_SECOND_ORDER, MFC_FHAT_SLIDING_WINDOW,
%   MFC_SISO.STEP.

% Growing-window annihilator, discretized (backward differences)
num_raw = -z + (t*z - (t - Ts)*state.z_km1)/Ts - t*alpha*u_prev - t*b_fold*z;
den_raw = t;

% Smooth numerator and denominator identically before dividing
num_filt = mfc_iir_smoother(num_raw, state.num_filt_km1, state.num_filt_km2, filter_window);
den_filt = mfc_iir_smoother(den_raw, state.den_filt_km1, state.den_filt_km2, filter_window);

valid = (den_filt ~= 0) && (t > hold_time);
if valid
    F_hat = num_filt / den_filt;
else
    F_hat = 0;
end

% Advance estimator state
state.z_km2        = state.z_km1;
state.z_km1        = z;
state.num_filt_km2 = state.num_filt_km1;
state.num_filt_km1 = num_filt;
state.den_filt_km2 = state.den_filt_km1;
state.den_filt_km1 = den_filt;

% 'integral' is zero-filled so dbg has the same struct type as the other
% mfc_fhat_* variants (required for code generation of the dispatch).
dbg = struct('num_raw', num_raw, 'den_raw', den_raw, ...
             'num_filt', num_filt, 'den_filt', den_filt, ...
             'integral', 0, 'valid', valid);
end
