function kernel = mfc_window_kernel(model_order, window_samples, Ts)
%MFC_WINDOW_KERNEL Precompute the sliding-window F-estimator quadrature kernel.
%
%   kernel = MFC_WINDOW_KERNEL(model_order, window_samples, Ts)
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
%   QUADRATURE: composite Simpson with an EVEN interval count (window_samples
%   is rounded up to even; realized window Tw = n_intervals*Ts). Trapezoidal
%   integration is NOT usable for the second-order kernel: the 60/Tw^5
%   prefactor amplifies its O(Ts^2/Tw^4) leakage to ~60x error at practical
%   sample times. Simpson is exact through cubics, so the sliding-window
%   estimate matches the algebraic variants.
%
%   Inputs
%     model_order    : 1 or 2 (order of the ultra-local model)
%     window_samples : requested window length [samples]; rounded up to even
%     Ts             : sample time [s]
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
%
%   See also MFC_FHAT_SLIDING_WINDOW, MFC_SISO_STEP.

validateattributes(model_order, {'numeric'}, {'scalar'});
assert(model_order == 1 || model_order == 2, ...
       'mfc_window_kernel: model_order must be 1 or 2.');
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
