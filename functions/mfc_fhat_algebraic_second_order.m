function [F_hat, state, dbg] = mfc_fhat_algebraic_second_order( ...
    z, u_prev, alpha, t, Ts, filter_window, hold_time, a_fold, b_fold, state)
%MFC_FHAT_ALGEBRAIC_SECOND_ORDER Algebraic (growing-window) F estimator, 2nd order.
%
%   [F_hat, state, dbg] = MFC_FHAT_ALGEBRAIC_SECOND_ORDER( ...
%       z, u_prev, alpha, t, Ts, filter_window, hold_time, a_fold, b_fold, state)
%
%   Estimates F in the second-order ultra-local model
%
%       ddot_z = F + alpha * u + a_fold * dot_z + b_fold * z
%
%   using the algebraic (operational-calculus) method: the model is written
%   in the Laplace domain, differentiated twice w.r.t. s to annihilate both
%   unknown initial conditions, and mapped back to the time domain, where
%   multiplication by t^n replaces (-d/ds)^n. Discretized with backward
%   differences over the growing window (t counted from run/sim start):
%
%       s_dz    = -( t*z - (t-Ts)*z[k-1] ) / Ts                    %  d/ds (sZ)
%       s2_d2z  =  ( t^2*z - 2(t-Ts)^2*z[k-1] + (t-2Ts)^2*z[k-2] ) / Ts^2
%       s_d2z   =  ( t^2*z - (t-Ts)^2*z[k-1] ) / Ts
%       dz      = -t*z
%       d2z     =  t^2*z
%       d2u     =  t^2*u_prev
%
%       num[k]  = 2z + 4*s_dz + s2_d2z - a_fold*(2*dz + s_d2z) ...
%                 - b_fold*d2z - alpha*d2u
%       den[k]  = t^2
%
%   Both are smoothed with the shared second-order IIR (MFC_IIR_SMOOTHER,
%   memory = filter_window samples) before the division, and the estimate is
%   held at zero until t > hold_time so early transients cannot blow up the
%   near-zero denominator.
%
%   COUPLED vs DECOUPLED use (selected by the caller):
%     * decoupled: z = y_measured, a_fold = b_fold = 0
%       -> F_hat is the pure-plant lumped dynamics ddot_y - alpha*u; the
%          stabilizing PD(I) feedback must be added explicitly in the
%          command law.
%     * coupled:   z = tracking error, a_fold = -Kd, b_fold = -Kp
%       -> the closed-loop characteristic polynomial s^2 + Kd*s + Kp is
%          folded into the estimate, so the command law needs no explicit
%          proportional/derivative feedback (u = (-F_hat + ddot_ref)/alpha).
%
%   Inputs
%     z             : estimator drive signal at step k (measurement or error)
%     u_prev        : command actually applied over the last sample
%     alpha         : ultra-local model input gain
%     t             : current time [s] (from run/sim start; growing window)
%     Ts            : sample time [s]
%     filter_window : IIR smoother memory [samples]
%     hold_time     : F_hat is forced to 0 while t <= hold_time [s]
%     a_fold        : folded dot_z coefficient  (coupled: -Kd, decoupled: 0)
%     b_fold        : folded z coefficient      (coupled: -Kp, decoupled: 0)
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
%   See also MFC_FHAT_ALGEBRAIC_FIRST_ORDER, MFC_FHAT_SLIDING_WINDOW,
%   MFC_SISO_STEP.

% Operational-calculus terms, discretized (backward differences)
s_dz   = -(t*z - (t - Ts)*state.z_km1) / Ts;
s2_d2z =  (t^2*z - 2*(t - Ts)^2*state.z_km1 + (t - 2*Ts)^2*state.z_km2) / Ts^2;
s_d2z  =  (t^2*z - (t - Ts)^2*state.z_km1) / Ts;
dz     = -t*z;
d2z    =  t^2*z;
d2u    =  t^2*u_prev;

num_raw = 2*z + 4*s_dz + s2_d2z - a_fold*(2*dz + s_d2z) - b_fold*d2z - alpha*d2u;
den_raw = t^2;

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
