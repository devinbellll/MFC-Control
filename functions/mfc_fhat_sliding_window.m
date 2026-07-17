function [F_hat, state, dbg] = mfc_fhat_sliding_window(y, u_prev, alpha, t, kernel, state)
%MFC_FHAT_SLIDING_WINDOW Sliding-window (non-algebraic) F estimator, 1st/2nd order.
%
%   [F_hat, state, dbg] = MFC_FHAT_SLIDING_WINDOW(y, u_prev, alpha, t, kernel, state)
%
%   Estimates F in the ultra-local model
%
%       dot_y  = F + alpha*u        (kernel.model_order = 1, Eq. 11)
%       ddot_y = F + alpha*u        (kernel.model_order = 2, Eq. 16)
%
%   by evaluating a fixed-length weighted integral of the measurement and
%   applied-input histories over the last Tw seconds (composite Simpson
%   quadrature; see MFC_WINDOW_KERNEL for the exact formulas and why
%   trapezoidal integration is not usable here):
%
%       F_hat = prefactor * (Ts/3) * sum( simpson_weights .* ...
%                 ( y_kernel .* y_buf + alpha * u_kernel_unit .* u_buf ) )
%
%   This estimator is DECOUPLED by construction: it sees only (y, u, alpha),
%   never the tracking error, so F_hat is the TRUE plant lumped dynamics
%   (dot_y or ddot_y) - alpha*u, and the stabilizing feedback must be added
%   explicitly in the command law. Unlike the algebraic (growing-window)
%   variants it has finite memory and no growing time weights, so it is
%   insensitive to the choice of time origin.
%
%   The estimate is held at zero until the window has filled (t > Tw).
%
%   Inputs
%     y      : current measurement
%     u_prev : command actually applied over the last sample (newest window
%              entry -- the current command is not known yet)
%     alpha  : ultra-local model input gain (may vary at run time; the
%              alpha-independent kernel is precomputed)
%     t      : current time [s], only used for the window-fill hold
%     kernel : precomputed quadrature kernel from MFC_WINDOW_KERNEL
%     state  : struct, fields used/updated here:
%                .y_buf [(n+1)x1] measurement window, newest last
%                .u_buf [(n+1)x1] applied-input window, newest last
%
%   Outputs
%     F_hat : estimate of F (0 while the window is filling)
%     state : updated state struct
%     dbg   : debug struct: integral (before prefactor), valid
%
%   See also MFC_WINDOW_KERNEL, MFC_FHAT_ALGEBRAIC_FIRST_ORDER,
%   MFC_FHAT_ALGEBRAIC_SECOND_ORDER, MFC_SISO_STEP.

% Slide the window buffers (newest sample last)
state.y_buf = [state.y_buf(2:end); y];
state.u_buf = [state.u_buf(2:end); u_prev];

% Composite Simpson quadrature of the weighted integrand
integrand = kernel.y_kernel .* state.y_buf + alpha * kernel.u_kernel_unit .* state.u_buf;
integral  = (kernel.Ts/3) * sum(kernel.simpson_weights .* integrand);

valid = t > kernel.Tw;                 % hold until the window has filled
if valid
    F_hat = kernel.prefactor * integral;
else
    F_hat = 0;
end

dbg = struct('integral', integral, 'valid', valid);
end
