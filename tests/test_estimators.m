function test_estimators()
%TEST_ESTIMATORS Closed-loop tests of the MFC System objects (MATLAB only).
%
%   >> cd <repo root>;  setup;  test_estimators
%
%   Drives the REAL matlab.System objects through step() in closed loop.
%   Requires MATLAB -- Octave cannot instantiate matlab.System. For a
%   MATLAB-free check of the same math use tests/octave_sanity.m,
%   tests/test_golden.m and tests/test_composed.m, which cover the numerics;
%   this file covers the OBJECT layer: port wiring, state reset, mask
%   parameter plumbing, and the stage blocks agreeing with mfc_siso_core.
%
%   Covers
%     1. mfc_siso_core across the supported variant grid: step tracking and
%        steady-state F accuracy/sign against an analytic value
%     2. rejection of the undefined coupled + sliding-window combination
%     3. sign-flip regression on an open-loop-UNSTABLE plant that only stays
%        bounded when F is estimated with the correct sign
%     4. the stage blocks, wired by hand, reproducing mfc_siso_core exactly
%     5. reset() genuinely restoring the initial state
%
%   See also MFC_SISO_CORE, OCTAVE_SANITY, TEST_GOLDEN, TEST_COMPOSED.

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..', 'functions'));
    addpath(fullfile(here, '..', 'blocks'));

    global N_PASS N_FAIL;  N_PASS = 0;  N_FAIL = 0;
    Ts = 1e-3;

    ORD = {'First order (dot_y = F + alpha*u)', 'Second order (ddot_y = F + alpha*u)'};
    STR = {'Coupled (error-driven estimator, poles folded)', ...
           'Decoupled (measurement-driven estimator, explicit iP/iPD)'};
    EST = {'Algebraic (growing window)', 'Sliding window (Simpson quadrature)'};

    % ---------------------------------------------------------------------
    % 1) Every supported 2nd-order variant tracks a step on a stable plant.
    %    Plant: yddot = -a1*yd - a0*y + b*u,  alpha = b (matched).
    %    Steady state: ydd = yd = 0 -> u_ss = a0*ref/b, F_true = -alpha*u_ss.
    % ---------------------------------------------------------------------
    a0 = 1; a1 = 1.4; b = 2; alpha = 2; ref = 1;
    Ftrue = -alpha * (a0*ref/b);

    variants = { ...
        ORD{2}, STR{1}, EST{1}, '2nd coupled algebraic'; ...
        ORD{2}, STR{2}, EST{1}, '2nd decoupled algebraic'; ...
        ORD{2}, STR{2}, EST{2}, '2nd decoupled window'};

    for v = 1:size(variants, 1)
        c = mfc_siso_core('model_order', variants{v,1}, ...
                          'controller_structure', variants{v,2}, ...
                          'estimator_type', variants{v,3}, ...
                          'Ts', Ts, 'alpha', alpha, 'Kp', 25, 'Kd', 10, 'Ki', 0, ...
                          'ref_filter_window', 10, 'est_filter_window', 10);
        [t, y, ~, F] = sim2(c, a0, a1, b, Ts, 6, @(tt) ref*(tt >= 0.2));
        fin = t > 5;
        nm  = variants{v,4};
        check(sprintf('%s: step tracking', nm), abs(mean(y(fin)) - ref) < 0.05);
        check(sprintf('%s: bounded',       nm), max(abs(y)) < 5);
        check(sprintf('%s: steady F=%.3f vs true %.3f', nm, mean(F(fin)), Ftrue), ...
              abs(mean(F(fin)) - Ftrue) < 0.15);
        check(sprintf('%s: F sign correct', nm), mean(F(fin)) < 0);
    end

    % ---------------------------------------------------------------------
    % 2) The undefined corner of the grid must be REJECTED, not silently
    %    mis-run: folding has no meaning for an estimator that never sees
    %    the tracking error.
    % ---------------------------------------------------------------------
    threw = false;
    try
        cbad = mfc_siso_core('model_order', ORD{2}, 'controller_structure', STR{1}, ...
                             'estimator_type', EST{2}, 'Ts', Ts);
        step(cbad, 0, 0, 0);
    catch
        threw = true;
    end
    check('coupled + sliding window is rejected', threw);

    % ---------------------------------------------------------------------
    % 3) SIGN-FLIP REGRESSION: open-loop-unstable plant with weak PD, so
    %    correct F cancellation is the ONLY thing keeping it bounded.
    %    yddot = +4*y + 2*u.  Correct (+60/Tw^5): bounded. Bug: diverges.
    % ---------------------------------------------------------------------
    cu = mfc_siso_core('model_order', ORD{2}, 'controller_structure', STR{2}, ...
                       'estimator_type', EST{2}, 'Ts', Ts, 'alpha', 2, ...
                       'Kp', 1, 'Kd', 2, 'Ki', 0, ...
                       'ref_filter_window', 10, 'est_filter_window', 10);
    [~, yu] = sim2(cu, -4, 0, 2, Ts, 6, @(tt) 1.0*(tt >= 0.2));
    check(sprintf('unstable plant stays bounded (max|y|=%.3g)', max(abs(yu))), ...
          max(abs(yu)) < 10);

    % ---------------------------------------------------------------------
    % 4) First-order decoupled variants on a first-order plant.
    %    Plant: ydot = -a0*y + b*u,  alpha = b.  F_true = -alpha*u_ss.
    % ---------------------------------------------------------------------
    a0 = 1; b = 2; alpha = 2; ref = 1;
    F1true = -alpha * (a0*ref/b);
    for v = 1:2
        c1 = mfc_siso_core('model_order', ORD{1}, 'controller_structure', STR{2}, ...
                           'estimator_type', EST{v}, 'Ts', Ts, 'alpha', alpha, ...
                           'Kp', 5, 'Kd', 0, 'Ki', 5, ...
                           'ref_filter_window', 10, 'est_filter_window', 10);
        [t1, y1, ~, F1] = sim1(c1, a0, b, Ts, 6, @(tt) ref*(tt >= 0.2));
        fin1 = t1 > 5;
        check(sprintf('1st decoupled (%s): step tracking', EST{v}), ...
              abs(mean(y1(fin1)) - ref) < 0.05);
        check(sprintf('1st decoupled (%s): steady F=%.3f vs true %.3f', ...
              EST{v}, mean(F1(fin1)), F1true), abs(mean(F1(fin1)) - F1true) < 0.15);
    end

    % ---------------------------------------------------------------------
    % 5) DECOMPOSITION at the object level: the stage blocks wired by hand
    %    must reproduce mfc_siso_core sample for sample. This is the twin of
    %    tests/test_composed.m, which proves the same for the stage maths.
    %    No saturation here, so there is no anti-windup delay difference and
    %    the match is exact.
    % ---------------------------------------------------------------------
    a0 = 1; a1 = 1.4; b = 2; alpha = 2; Kp = 25; Kd = 10; Ki = 0;

    mono = mfc_siso_core('model_order', ORD{2}, 'controller_structure', STR{2}, ...
                         'estimator_type', EST{1}, 'Ts', Ts, 'alpha', alpha, ...
                         'Kp', Kp, 'Kd', Kd, 'Ki', Ki, ...
                         'ref_filter_window', 10, 'est_filter_window', 10);

    smoo = mfc_smoother_block('Ts', Ts, 'window', 10, 'output_derivatives', true);
    est2 = mfc_fhat_alg2_block('Ts', Ts, 'alpha', alpha, 'a_fold', 0, 'b_fold', 0, ...
                               'est_filter_window', 10, 'est_hold_time', 0.1);
    fbk  = mfc_feedback_block('Ts', Ts, 'Kp', Kp, 'Kd', Kd, 'Ki', Ki, 'coupled', false);
    cmd  = mfc_command_block('alpha', alpha);
    filt = mfc_command_filter_block('Ts', Ts, 'command_filter', 1, 'use_control_sat', false);

    tv = 0:Ts:3;
    ym = 0; dym = 0;  yc = 0; dyc = 0;  u_km1 = 0;  worst = 0;
    for k = 1:numel(tv)
        r = 1.0*(tv(k) >= 0.2);

        um = step(mono, r, ym, tv(k));               % all-in-one

        [sp, ~, ddsp] = step(smoo, r);               % assembled from stages
        e   = yc - sp;
        Fh  = step(est2, yc, u_km1, tv(k));          % u_km1 = the explicit unit delay
        fb  = step(fbk, e);
        ur  = step(cmd, Fh, ddsp, fb);
        uc  = step(filt, ur);
        u_km1 = uc;

        worst = max(worst, abs(um - uc));

        ydd = -a1*dym - a0*ym + b*um;  dym = dym + ydd*Ts;  ym = ym + dym*Ts;
        ydd = -a1*dyc - a0*yc + b*uc;  dyc = dyc + ydd*Ts;  yc = yc + dyc*Ts;
    end
    check(sprintf('stage blocks reproduce mfc_siso_core (max|du|=%.3g)', worst), ...
          worst < 1e-12);

    % ---------------------------------------------------------------------
    % 6) reset() must genuinely restore the initial state: the same input
    %    sequence after a reset must give the same output sequence.
    % ---------------------------------------------------------------------
    cr = mfc_siso_core('model_order', ORD{2}, 'controller_structure', STR{2}, ...
                       'estimator_type', EST{1}, 'Ts', Ts, 'alpha', 2, ...
                       'ref_filter_window', 10, 'est_filter_window', 10);
    first = zeros(1, 50);
    for k = 1:50, first(k) = step(cr, 1, 0.1*k, (k-1)*Ts); end
    reset(cr);
    again = zeros(1, 50);
    for k = 1:50, again(k) = step(cr, 1, 0.1*k, (k-1)*Ts); end
    check('reset() restores the initial state exactly', isequal(first, again));

    fprintf('\n%d passed, %d failed\n', N_PASS, N_FAIL);
    if N_FAIL > 0
        error('test_estimators: %d test(s) failed', N_FAIL);
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

function [t, y, U, F] = sim2(ctrl, a0, a1, b, Ts, Tend, ref)
% Closed loop on yddot = -a1*yd - a0*y + b*u (ZOH on u, forward Euler).
    t = 0:Ts:Tend;  y = 0;  yd = 0;
    Y = zeros(size(t));  U = Y;  F = Y;
    for k = 1:numel(t)
        [u, f] = step(ctrl, ref(t(k)), y, t(k));   % measure y(t_k), compute u_k
        Y(k) = y;  U(k) = u;  F(k) = f;
        ydd = -a1*yd - a0*y + b*u;                 % apply u_k over [t_k, t_k+Ts]
        yd  = yd + ydd*Ts;
        y   = y  + yd*Ts;
    end
    y = Y;
end

function [t, y, U, F] = sim1(ctrl, a0, b, Ts, Tend, ref)
% Closed loop on ydot = -a0*y + b*u (ZOH on u, forward Euler).
    t = 0:Ts:Tend;  y = 0;
    Y = zeros(size(t));  U = Y;  F = Y;
    for k = 1:numel(t)
        [u, f] = step(ctrl, ref(t(k)), y, t(k));
        Y(k) = y;  U(k) = u;  F(k) = f;
        yd = -a0*y + b*u;
        y  = y + yd*Ts;
    end
    y = Y;
end
