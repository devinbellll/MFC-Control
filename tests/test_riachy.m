function test_riachy()
%TEST_RIACHY Riachy's trick: the transform, the estimate, and the loop.
%
%   octave --no-gui -q --path tests --path functions --eval test_riachy
%   (also runs in MATLAB)
%
%   Riachy's trick replaces the derivative feedback of a second-order iPD
%   by an INTEGRAL of the measurement:
%
%       Y = y + Kd*int y     =>     ddot_Y = Fk + alpha*u,   Fk = F + Kd*dot_y
%
%   Three things have to hold, and this file checks each of them directly:
%
%     1. the discrete transform is the trapezoidal integral it claims to be
%     2. a standard 2nd-order estimator driven by Y returns Fk, not F --
%        i.e. the derivative term really does come back out of the estimate
%     3. the composed loop (Riachy estimator + a P/PI feedback with NO
%        derivative term + ff = ddot_sp + Kd*dot_sp) tracks a step, and
%        tracks it like the ordinary decoupled iPD it is meant to replace
%
%   Section 3 is the load-bearing one: it is the only place the wiring
%   invariant "D = 0 on the PID, Kd on the feedforward" is asserted rather
%   than merely documented.
%
%   Plain functions only, so it runs without MATLAB. The block layer
%   (mfc_fhat_riachy2_block: ports, mask, reset, estimator selector) is
%   covered by tests/test_estimators.m section 6.
%
%   See also MFC_RIACHY_TRANSFORM, MFC_FHAT_RIACHY2_BLOCK, OCTAVE_SANITY,
%   TEST_COMPOSED.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'functions'));

global N_PASS N_FAIL;  N_PASS = 0;  N_FAIL = 0;
Ts = 1e-3;

% =====================================================================
% 1) The transform is a trapezoidal integral of y
% =====================================================================
% Constant y = c, stepped N times from t = 0. The first sample only counts
% half (y_km1 starts at 0), so int[N] = c*(N - 1/2)*Ts: the continuous
% integral shifted by half a sample. That offset is HARMLESS and expected --
% the lower limit c of int_c^t is free in the derivation, and a constant
% offset in the integral does not change ddot_Y at all.
c = 2.5;  Kd = 4;  N = 100;
rs = struct('int_km1', 0, 'y_km1', 0);
for k = 1:N
    [Y, rs] = mfc_riachy_transform(c, Kd, Ts, rs);
end
int_expect = c*(N - 0.5)*Ts;
check(sprintf('transform: int y = %.6f, expected %.6f (trapezoid, half-sample offset)', ...
      rs.int_km1, int_expect), abs(rs.int_km1 - int_expect) < 1e-12);
check('transform: Y = y + Kd*int y', abs(Y - (c + Kd*int_expect)) < 1e-12);

% Linear y = v*t is integrated EXACTLY by the trapezoid, up to the same
% half-sample offset -- this is why the transform is trapezoidal and not
% forward Euler.
v = 3;  rs = struct('int_km1', 0, 'y_km1', 0);
for k = 1:N
    t = (k-1)*Ts;
    [~, rs] = mfc_riachy_transform(v*t, Kd, Ts, rs);
end
t_end = (N-1)*Ts;
check(sprintf('transform: linear y integrated exactly (%.3e)', ...
      abs(rs.int_km1 - 0.5*v*t_end^2)), abs(rs.int_km1 - 0.5*v*t_end^2) < 1e-12);

% =====================================================================
% 2) An estimator driven by Y returns Fk = F + Kd*dot_y
% =====================================================================
% Plant with an analytically known F: y = a*t^2/2 under a constant input,
% so ddot_y = a, F = a - alpha*u0, and dot_y = a*t. Fk therefore RAMPS,
% which is the point -- if the estimator returned F the check would fail by
% the whole Kd*a*t term.
a = 3;  alpha = 2;  u0 = 0.5;  Kd = 4;  t_eval = 0.5;
F_true  = a - alpha*u0;
Fk_true = F_true + Kd*a*t_eval;

% Fk RAMPS, so the estimate lags by half the window: a fixed-length window
% returns the average over [t-Tw, t], which for a linear Fk is exactly
% Fk(t - Tw/2). That is a property of the quadrature, not of the trick, so
% it is checked as an equality rather than hidden in a tolerance.
Tw       = 40*Ts;
Fk_centre = F_true + Kd*a*(t_eval - Tw/2);

Fw = run_riachy(false, Ts, alpha, Kd, a, u0, t_eval);
check(sprintf('window Riachy: Fk=%.4f vs Fk(t-Tw/2)=%.4f (F alone would be %.4f)', ...
      Fw, Fk_centre, F_true), abs(Fw - Fk_centre) < 1e-3);
check(sprintf('window Riachy: within half a window of Fk(t)=%.4f', Fk_true), ...
      abs(Fw - Fk_true) < Kd*a*Tw);
check('window Riachy: the Kd*dot_y term is genuinely present', ...
      abs(Fw - F_true) > 1);

Fa = run_riachy(true, Ts, alpha, Kd, a, u0, t_eval);
% Algebraic: a GROWING window on a ramping Fk, so it lags more than the
% fixed-memory one and by an amount that depends on t. Hence the loose
% tolerance here -- and hence the block help recommending the sliding
% window for this trick.
check(sprintf('algebraic Riachy: Fk=%.4f vs true %.4f', Fa, Fk_true), ...
      abs(Fa - Fk_true) < 0.5);
check('algebraic Riachy: the Kd*dot_y term is genuinely present', ...
      abs(Fa - F_true) > 1);

% =====================================================================
% 3) The composed loop: no derivative of y anywhere, and it still tracks
% =====================================================================
% Plant ddot_y = -a1*dot_y - a0*y + b*u, closed with
%   Riachy estimator -> F_hat = Fk
%   fb = Kp*err                 (P only: D = 0, Kd lives in Y and in ff)
%   ff = ddot_sp + Kd*dot_sp
%   u  = (-F_hat + ff - fb)/alpha
Ts3 = 1e-3;  a0 = 4;  a1 = 0.5;  b = 1.2;  alpha3 = b;
p = 8;  Kp = p^2;  Kd3 = 2*p;                 % double pole at -p
t3 = (0:Ts3:3)';  ref = double(t3 >= 0.1);

[y_r, u_r] = sim_riachy(t3, Ts3, a0, a1, b, alpha3, Kp, Kd3, ref);
check(sprintf('riachy loop: tracks the step (final y = %.4f)', y_r(end)), ...
      abs(y_r(end) - 1) < 2e-2);
check('riachy loop: stays bounded', all(isfinite(y_r)) && max(abs(y_r)) < 3);
check('riachy loop: command stays bounded', all(isfinite(u_r)) && max(abs(u_r)) < 500);

% Same plant, same gains, the ORDINARY decoupled iPD (explicit Kd*dot_err
% from the feedback stage). Riachy's trick is supposed to reproduce this
% response without ever differentiating y.
[y_c, ~] = sim_classic(t3, Ts3, a0, a1, b, alpha3, Kp, Kd3, ref);
check(sprintf('riachy loop: matches the classic iPD response (max|dy| = %.4f)', ...
      max(abs(y_r - y_c))), max(abs(y_r - y_c)) < 8e-2);

% The wiring invariant, stated as a test: dropping Kd*dot_sp from the
% feedforward must visibly cost tracking on a MOVING setpoint (it is
% invisible on a settled step, which is exactly why it gets left out).
ramp = min(1, max(0, (t3 - 0.1)/2));
[y_ff,  ~] = sim_riachy(t3, Ts3, a0, a1, b, alpha3, Kp, Kd3, ramp);
[y_nff, ~] = sim_riachy(t3, Ts3, a0, a1, b, alpha3, Kp, Kd3, ramp, false);
e_ff  = max(abs(y_ff(end-500:end)  - ramp(end-500:end)));
e_nff = max(abs(y_nff(end-500:end) - ramp(end-500:end)));
check(sprintf('riachy loop: ff must carry Kd*dot_sp (err %.2e with, %.2e without)', ...
      e_ff, e_nff), e_nff > 5*e_ff);

fprintf('\ntest_riachy: %d passed, %d FAILED\n', N_PASS, N_FAIL);
if N_FAIL > 0, error('test_riachy: %d check(s) failed.', N_FAIL); end
end


% ===== helpers =======================================================

function check(name, cond)
global N_PASS N_FAIL;
if cond
    N_PASS = N_PASS + 1;  fprintf('  ok   %s\n', name);
else
    N_FAIL = N_FAIL + 1;  fprintf('  FAIL %s\n', name);
end
end

function F = run_riachy(algebraic, Ts, alpha, Kd, a, u0, t_eval)
% Drive a Riachy estimator with the exact samples of y = a*t^2/2, u = u0.
kern = mfc_siso.window_kernel(2, 40, Ts);
rs   = struct('int_km1', 0, 'y_km1', 0);
est  = struct('z_km1', 0, 'z_km2', 0, ...
              'num_filt_km1', 0, 'num_filt_km2', 0, ...
              'den_filt_km1', 0, 'den_filt_km2', 0, ...
              'y_buf', zeros(kern.n_intervals+1, 1), ...
              'u_buf', zeros(kern.n_intervals+1, 1));
F = 0;
for k = 1:round(t_eval/Ts) + 1
    t = (k-1)*Ts;
    [Y, rs] = mfc_riachy_transform(0.5*a*t^2, Kd, Ts, rs);
    if algebraic
        [F, est] = mfc_fhat_algebraic_second_order(Y, u0, alpha, t, Ts, 10, 0.1, 0, 0, est);
    else
        [F, est] = mfc_fhat_sliding_window(Y, u0, alpha, t, kern, est);
    end
end
end

function [y, U] = sim_riachy(t, Ts, a0, a1, b, alpha, Kp, Kd, ref, use_kd_ff)
% Composed Riachy loop (sliding-window estimator) on ddot_y = -a1*dot_y - a0*y + b*u.
if nargin < 10, use_kd_ff = true; end
kern = mfc_siso.window_kernel(2, 40, Ts);
rs   = struct('int_km1', 0, 'y_km1', 0);
est  = struct('z_km1', 0, 'z_km2', 0, ...
              'num_filt_km1', 0, 'num_filt_km2', 0, ...
              'den_filt_km1', 0, 'den_filt_km2', 0, ...
              'y_buf', zeros(kern.n_intervals+1, 1), ...
              'u_buf', zeros(kern.n_intervals+1, 1));
W_REF = 100;                              % 0.1 s reference smoothing
sp1 = 0;  sp2 = 0;  u_prev = 0;
yk = 0;  dyk = 0;
N = numel(t);  y = zeros(N, 1);  U = zeros(N, 1);
for k = 1:N
    y(k) = yk;
    [sp_filt, dot_sp, ddot_sp] = mfc_siso.ref_traj(ref(k), sp1, sp2, Ts, W_REF, true);
    sp2 = sp1;  sp1 = sp_filt;

    [Y, rs]  = mfc_riachy_transform(yk, Kd, Ts, rs);
    [F, est] = mfc_fhat_sliding_window(Y, u_prev, alpha, t(k), kern, est);

    err = yk - sp_filt;
    fb  = Kp*err;                          % P only -- no derivative anywhere
    if use_kd_ff
        ff = ddot_sp + Kd*dot_sp;
    else
        ff = ddot_sp;
    end
    u    = mfc_siso.command(F, ff, fb, alpha);
    U(k) = u;

    ddy = -a1*dyk - a0*yk + b*u;           % ZOH on u, forward Euler
    yk  = yk + Ts*dyk;
    dyk = dyk + Ts*ddy;
    u_prev = u;
end
end

function [y, U] = sim_classic(t, Ts, a0, a1, b, alpha, Kp, Kd, ref)
% The ordinary decoupled iPD on the same plant: estimator on y, explicit
% Kp*err + Kd*dot_err feedback. The response Riachy's trick must reproduce.
kern = mfc_siso.window_kernel(2, 40, Ts);
est  = struct('y_buf', zeros(kern.n_intervals+1, 1), ...
              'u_buf', zeros(kern.n_intervals+1, 1));
W_REF = 100;
sp1 = 0;  sp2 = 0;  u_prev = 0;  err_km1 = 0;  int_km1 = 0;
yk = 0;  dyk = 0;
N = numel(t);  y = zeros(N, 1);  U = zeros(N, 1);
for k = 1:N
    y(k) = yk;
    [sp_filt, ~, ddot_sp] = mfc_siso.ref_traj(ref(k), sp1, sp2, Ts, W_REF, true);
    sp2 = sp1;  sp1 = sp_filt;

    [F, est] = mfc_fhat_sliding_window(yk, u_prev, alpha, t(k), kern, est);

    err = yk - sp_filt;
    [fb, int_km1] = mfc_siso.feedback(err, err_km1, int_km1, Ts, Kp, Kd, 0, false);
    err_km1 = err;

    u    = mfc_siso.command(F, ddot_sp, fb, alpha);
    U(k) = u;

    ddy = -a1*dyk - a0*yk + b*u;
    yk  = yk + Ts*dyk;
    dyk = dyk + Ts*ddy;
    u_prev = u;
end
end
