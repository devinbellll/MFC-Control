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
%   applied-input histories over the last Tw seconds. That integral is
%   precomputed into one fixed multiplier -- a TAP -- per stored sample
%   (see MFC_SISO.WINDOW_KERNEL for the kernels and how the taps are
%   derived), so the estimator is a pure FIR filter:
%
%       F_hat = tap_y' * y_buf  +  alpha * ( tap_u_unit' * u_buf )
%
%   No poles, no recursion, no stability question: a sample enters, is
%   weighted for exactly Tw seconds, and then leaves completely.
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
%     kernel : precomputed taps from MFC_SISO.WINDOW_KERNEL
%     state  : struct, fields used/updated here:
%                .y_buf [(n+1)x1] measurement window, newest last
%                .u_buf [(n+1)x1] applied-input window, newest last
%
%   Outputs
%     F_hat : estimate of F (0 while the window is filling)
%     state : updated state struct
%     dbg   : debug struct: integral (the raw tap sum, i.e. F_hat before
%             the startup hold is applied), valid
%
%   See also MFC_SISO.WINDOW_KERNEL, MFC_FHAT_ALGEBRAIC_FIRST_ORDER,
%   MFC_FHAT_ALGEBRAIC_SECOND_ORDER, MFC_SISO.STEP.

% Slide the window buffers (newest sample last). circshift + end-assignment
% instead of [buf(2:end); new]: identical result, but also valid for the 1x1
% placeholder buffers of algebraic variants, whose (dead) copy of this code
% is still compiled under code generation.
state.y_buf        = circshift(state.y_buf, -1);
state.y_buf(end)   = y;
state.u_buf        = circshift(state.u_buf, -1);
state.u_buf(end)   = u_prev;

% One multiply-accumulate per window. The taps already carry the
% prefactor; alpha is kept out of tap_u_unit so it can vary at run time.
integral = kernel.tap_y.' * state.y_buf ...
           + alpha * (kernel.tap_u_unit.' * state.u_buf);

valid = t > kernel.Tw;                 % hold until the window has filled
if valid
    F_hat = integral;
else
    F_hat = 0;
end

% num/den fields are zero-filled so dbg has the same struct type as the
% algebraic mfc_fhat_* variants (required for code generation of the dispatch).
dbg = struct('num_raw', 0, 'den_raw', 0, 'num_filt', 0, 'den_filt', 0, ...
             'integral', integral, 'valid', valid);
end
