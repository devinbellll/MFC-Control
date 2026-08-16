function octave_sanity()
%OCTAVE_SANITY Estimator and smoother math vs analytically known answers.
%
%   octave --no-gui -q tests/octave_sanity.m     (also runs in MATLAB)
%
%   Checks the numerics against plants whose F is known in closed form:
%   feed each estimator the EXACT sampled signals of a plant built to have a
%   chosen F, and demand it recovers that F with the right magnitude AND the
%   right sign.
%
%   This file used to carry its own hand-written copy of the estimator math
%   because the algorithms lived inside matlab.System classes that Octave
%   cannot instantiate. They now live in plain functions, so it calls the
%   REAL code -- there is no second copy to keep in sync any more.
%
%   Sign guards. Two historical bugs are pinned here explicitly:
%     - 2nd-order sliding-window prefactor: -60/Tw^5 returns -F
%     - 1st-order sliding-window input term: a negated u kernel returns
%       F + (extra)*u
%   Both show up as a sign flip, not a magnitude error, so every check below
%   asserts the sign separately from the tolerance.
%
%   Sections 7-8 cover the MIMO (matrix-alpha) pair directly: den_raw/den_filt
%   staying scalar while num_raw/num_filt are per-channel, a diagonal alpha
%   not cross-coupling channels, and the command block's matrix solve
%   reducing to the scalar division. The closed-loop equivalent lives in
%   tests/test_golden_mimo.m.
%
%   See also TEST_GOLDEN, TEST_GOLDEN_MIMO, TEST_COMPOSED, TEST_ESTIMATORS.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'functions'));

global N_PASS N_FAIL;  N_PASS = 0;  N_FAIL = 0;
Ts = 1e-3;

% =====================================================================
% 1) mfc_iir_smoother: pass-through and unity DC gain
% =====================================================================
check('smoother: W=0 is an exact pass-through', ...
      mfc_iir_smoother(3.7, 1.1, 2.2, 0) == 3.7);

x1 = 0; x2 = 0;                       % settle a constant through the filter
for k = 1:2000
    xf = mfc_iir_smoother(5, x1, x2, 10);  x2 = x1;  x1 = xf;
end
check(sprintf('smoother: unity DC gain (settled %.10f vs 5)', x1), abs(x1 - 5) < 1e-9);

x1 = 0; x2 = 0;                       % monotone, no overshoot (critically damped)
prev = -inf;  mono = true;
for k = 1:400
    xf = mfc_iir_smoother(1, x1, x2, 10);  x2 = x1;  x1 = xf;
    mono = mono && (xf >= prev - 1e-15) && (xf <= 1 + 1e-12);
    prev = xf;
end
check('smoother: step response monotone, no overshoot', mono);

% =====================================================================
% 2) Sliding window, 2nd order: ddot_y = F + alpha*u
%    y = 0.5*a*t^2  =>  ddot_y = a  =>  F_true = a - alpha*u0
% =====================================================================
% u0 chosen so F_true is comfortably non-zero: a sign check against F_true = 0
% would be vacuous. Tolerance is 1e-3, not eps: the 2nd-order integrand is
% y_kernel (quadratic) x y (quadratic) = quartic, and Simpson is exact only
% through cubics, so a small quadrature residual is expected and correct.
a = 3; alpha = 2; u0 = 0.5;  Ftrue = a - alpha*u0;
F = run_window(2, 40, Ts, alpha, @(tt) 0.5*a*tt.^2, u0);
check(sprintf('window 2nd: F=%.6f vs true %.6f', F, Ftrue), abs(F - Ftrue) < 1e-3);
check('window 2nd: F sign correct (guards -60/Tw^5)', sign(F) == sign(Ftrue));

a = 3; alpha = 2; u0 = 3;  Ftrue = a - alpha*u0;      % force F_true < 0
F = run_window(2, 40, Ts, alpha, @(tt) 0.5*a*tt.^2, u0);
check(sprintf('window 2nd: negative F=%.6f vs true %.6f', F, Ftrue), abs(F - Ftrue) < 1e-3);
check('window 2nd: negative F sign correct', F < 0);

% =====================================================================
% 3) Sliding window, 1st order: dot_y = F + alpha*u
%    y = v*t  =>  dot_y = v  =>  F_true = v - alpha*u0
% =====================================================================
% 1st order: y_kernel (linear) x y (linear) = quadratic, which Simpson
% integrates exactly -- hence the much tighter tolerance than 2nd order.
v = 4; alpha = 2; u0 = 0.5;  Ftrue = v - alpha*u0;
F = run_window(1, 40, Ts, alpha, @(tt) v*tt, u0);
check(sprintf('window 1st: F=%.6f vs true %.6f', F, Ftrue), abs(F - Ftrue) < 1e-6);
check('window 1st: F sign correct (guards the u-kernel sign)', sign(F) == sign(Ftrue));

v = 1; alpha = 2; u0 = 3;  Ftrue = v - alpha*u0;      % force F_true < 0
F = run_window(1, 40, Ts, alpha, @(tt) v*tt, u0);
check(sprintf('window 1st: negative F=%.6f vs true %.6f', F, Ftrue), abs(F - Ftrue) < 1e-6);

% =====================================================================
% 4) Algebraic estimators, decoupled (no folding), constant F
%    Backward-difference discretization, so O(Ts) error is expected.
% =====================================================================
a = 3; alpha = 2; u0 = 0.5;  Ftrue = a - alpha*u0;
F = run_alg(2, Ts, alpha, @(tt) 0.5*a*tt.^2, u0);
check(sprintf('algebraic 2nd: F=%.4f vs true %.4f', F, Ftrue), abs(F - Ftrue) < 2e-2);
check('algebraic 2nd: F sign correct', sign(F) == sign(Ftrue));

v = 4; alpha = 2; u0 = 0.5;  Ftrue = v - alpha*u0;
F = run_alg(1, Ts, alpha, @(tt) v*tt, u0);
check(sprintf('algebraic 1st: F=%.4f vs true %.4f', F, Ftrue), abs(F - Ftrue) < 1e-2);
check('algebraic 1st: F sign correct', sign(F) == sign(Ftrue));

% =====================================================================
% 5) Cross-check: the two estimator families must agree on the same plant
% =====================================================================
a = 3; alpha = 2; u0 = 0.5;  Ftrue = a - alpha*u0;
Fw = run_window(2, 40, Ts, alpha, @(tt) 0.5*a*tt.^2, u0);
Fa = run_alg(2, Ts, alpha, @(tt) 0.5*a*tt.^2, u0);
check(sprintf('2nd order: window %.4f vs algebraic %.4f agree (both -> %.1f)', Fw, Fa, Ftrue), ...
      abs(Fw - Fa) < 2e-2 && abs(Fw - Ftrue) < 2e-2);

% =====================================================================
% 7) MIMO algebraic estimator (2nd order, decoupled): den_raw is SCALAR and
%    SHARED across the vector while num_raw is PER-ELEMENT. That asymmetry
%    is the whole reason one integration window can serve a vector channel,
%    and it is what a C port is most likely to get wrong -- so it is pinned
%    directly here, not only via the closed-loop golden traces in
%    tests/test_golden_mimo.m.
% =====================================================================
a1 = 3; a2 = -2; alpha_v = 2; u0 = 0.5;      % two DIFFERENT channel dynamics
Ftrue1 = a1 - alpha_v*u0;  Ftrue2 = a2 - alpha_v*u0;
[F2, dbg2] = run_alg_vec(2, Ts, alpha_v*eye(2), ...
    @(tt) [0.5*a1*tt.^2; 0.5*a2*tt.^2], [u0; u0]);
check(sprintf('MIMO algebraic 2nd: F=[%.4f %.4f] vs true [%.4f %.4f]', ...
      F2(1), F2(2), Ftrue1, Ftrue2), ...
      abs(F2(1) - Ftrue1) < 2e-2 && abs(F2(2) - Ftrue2) < 2e-2);
check('MIMO algebraic 2nd: den_raw/den_filt are scalar (shared)', ...
      isscalar(dbg2.den_raw) && isscalar(dbg2.den_filt));
check('MIMO algebraic 2nd: num_raw/num_filt are per-channel (2x1)', ...
      isequal(size(dbg2.num_raw), [2 1]) && isequal(size(dbg2.num_filt), [2 1]));

% Diagonal alpha must not cross-couple the channels: running two DIFFERENT
% signals through the vector estimator with alpha = eye(2) must match
% running each one through the scalar estimator independently, sample for
% sample (bit-identical -- multiplying by the identity's 1s and 0s
% introduces no rounding).
[Fv, ~] = run_alg_vec_trace(2, Ts, eye(2), ...
    @(tt) [0.5*a1*tt.^2; 0.5*a2*tt.^2], [u0; u0], 200);
Fs1 = run_alg_trace(Ts, 1, @(tt) 0.5*a1*tt.^2, u0, 200);
Fs2 = run_alg_trace(Ts, 1, @(tt) 0.5*a2*tt.^2, u0, 200);
check('MIMO algebraic 2nd: diagonal alpha does not cross-couple channels', ...
      isequal(Fv(:,1), Fs1) && isequal(Fv(:,2), Fs2));

% =====================================================================
% 8) MIMO command (matrix alpha): diagonal reduces to the scalar division;
%    off-diagonal genuinely exercises the linear solve (guards a silent
%    regression to elementwise division, which would be wrong whenever
%    alpha is not diagonal).
% =====================================================================
Fh = [1.5; -2.0];  ffv = [0.3; 0.1];  fbv = [0.2; -0.1];
alpha_diag2 = diag([2, 3]);
u_mimo  = mfc_siso.command_mimo(Fh, ffv, fbv, alpha_diag2);
u_scalar = [mfc_siso.command(Fh(1), ffv(1), fbv(1), 2); ...
            mfc_siso.command(Fh(2), ffv(2), fbv(2), 3)];
check(sprintf('MIMO command: diagonal alpha reduces to scalar command (max|d|=%.3g)', ...
      max(abs(u_mimo - u_scalar))), max(abs(u_mimo - u_scalar)) < 1e-9);

alpha_off = [1.0 0.4; -0.3 1.1];
u_solve   = mfc_siso.command_mimo(Fh, ffv, fbv, alpha_off);
u_naive   = (-Fh + ffv - fbv) ./ diag(alpha_off);   % the wrong, elementwise "port"
check('MIMO command: off-diagonal alpha genuinely couples channels (solve != elementwise divide)', ...
      max(abs(u_solve - u_naive)) > 1e-3);

% =====================================================================
% 9) Startup hold: both families must output exactly 0 before they are valid
% =====================================================================
st = zero_state(1);
[F0, ~, dbg0] = mfc_fhat_algebraic_second_order(1, 1, 1, 0, Ts, 10, 0.1, 0, 0, st);
check('algebraic: F is exactly 0 during the startup hold', F0 == 0 && ~dbg0.valid);

kern = mfc_siso.window_kernel(2, 40, Ts);
stw  = struct('y_buf', zeros(kern.n_intervals+1, 1), 'u_buf', zeros(kern.n_intervals+1, 1));
[Fw0, ~, dbgw] = mfc_fhat_sliding_window(1, 1, 1, 0, kern, stw);
check('window: F is exactly 0 while the window fills', Fw0 == 0 && ~dbgw.valid);

fprintf('\n%d passed, %d failed\n', N_PASS, N_FAIL);
if N_FAIL > 0
    error('octave_sanity: %d check(s) failed', N_FAIL);
end
end

% ===== helpers =====

function check(name, cond)
    global N_PASS N_FAIL;
    if cond
        N_PASS = N_PASS + 1;  fprintf('  PASS  %s\n', name);
    else
        N_FAIL = N_FAIL + 1;  fprintf('  FAIL  %s\n', name);
    end
end

function st = zero_state(n_buf)
    st = struct('z_km1', 0, 'z_km2', 0, ...
                'num_filt_km1', 0, 'num_filt_km2', 0, ...
                'den_filt_km1', 0, 'den_filt_km2', 0, ...
                'y_buf', zeros(n_buf, 1), 'u_buf', zeros(n_buf, 1));
end

function F = run_window(order, win, Ts, alpha, yfun, u0)
% Drive the sliding-window estimator with exact samples of yfun and a
% constant input, long enough for the window to fill, and return the final F.
    kern = mfc_siso.window_kernel(order, win, Ts);
    st   = struct('y_buf', zeros(kern.n_intervals+1, 1), ...
                  'u_buf', zeros(kern.n_intervals+1, 1));
    N = round(4*kern.Tw/Ts);
    for k = 1:N
        t = (k-1)*Ts;
        [F, st] = mfc_fhat_sliding_window(yfun(t), u0, alpha, t, kern, st);
    end
end

function F = run_alg(order, Ts, alpha, yfun, u0)
% Drive an algebraic estimator (decoupled: no folding) with exact samples of
% yfun and a constant input, and return the final F.
    st = zero_state(1);
    N  = 4000;
    for k = 1:N
        t = (k-1)*Ts;
        if order == 2
            [F, st] = mfc_fhat_algebraic_second_order(yfun(t), u0, alpha, t, Ts, 10, 0.1, 0, 0, st);
        else
            [F, st] = mfc_fhat_algebraic_first_order(yfun(t), u0, alpha, t, Ts, 10, 0.1, 0, st);
        end
    end
end

function [F, dbg] = run_alg_vec(n, Ts, alpha, yfun, u0)
% Drive the 2nd-order algebraic estimator (decoupled) with a vector z and a
% square alpha, for 4000 samples, and return the final F and its dbg struct
% (used to inspect num_raw/den_raw shapes).
    st = struct('z_km1', zeros(n,1), 'z_km2', zeros(n,1), ...
                'num_filt_km1', zeros(n,1), 'num_filt_km2', zeros(n,1), ...
                'den_filt_km1', 0, 'den_filt_km2', 0);
    N = 4000;
    for k = 1:N
        t = (k-1)*Ts;
        [F, st, dbg] = mfc_fhat_algebraic_second_order(yfun(t), u0, alpha, t, Ts, 10, 0.1, 0, 0, st);
    end
end

function [Ftrace, st] = run_alg_vec_trace(n, Ts, alpha, yfun, u0, N)
% Like run_alg_vec, but returns F at every sample as an [N x n] trace
% instead of only the final value.
    st = struct('z_km1', zeros(n,1), 'z_km2', zeros(n,1), ...
                'num_filt_km1', zeros(n,1), 'num_filt_km2', zeros(n,1), ...
                'den_filt_km1', 0, 'den_filt_km2', 0);
    Ftrace = zeros(N, n);
    for k = 1:N
        t = (k-1)*Ts;
        [F, st] = mfc_fhat_algebraic_second_order(yfun(t), u0, alpha, t, Ts, 10, 0.1, 0, 0, st);
        Ftrace(k, :) = F.';
    end
end

function Ftrace = run_alg_trace(Ts, alpha, yfun, u0, N)
% Like run_alg (2nd order, decoupled), but returns F at every sample as an
% [N x 1] trace instead of only the final value.
    st = zero_state(1);
    Ftrace = zeros(N, 1);
    for k = 1:N
        t = (k-1)*Ts;
        [F, st] = mfc_fhat_algebraic_second_order(yfun(t), u0, alpha, t, Ts, 10, 0.1, 0, 0, st);
        Ftrace(k) = F;
    end
end
