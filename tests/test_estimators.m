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
%     4. the first-order decoupled variants on a first-order plant
%     5. the stage blocks, wired by hand, reproducing mfc_siso_core exactly
%     6. reset() genuinely restoring the initial state
%     7. mfc_fhat_riachy2_block: ports, the algebraic/window selector, the
%        optional Y output, and it being a thin wrapper over
%        mfc_riachy_transform + the selected estimator
%     8. mfc_fhat_riachy2_mimo_block: the same wrapper claim with matrix Kd
%        and matrix alpha, plus port widths, the mask size checks, and the
%        reduction to the SISO block when both matrices are diagonal
%     9. the rest of the n-channel estimator surface -- the 1st-order
%        decoupled, the two coupled and the sliding-window MIMO blocks --
%        each against the kernel it claims to wrap, plus the n-channel
%        smoother
%    10. mfc_mimo_core: n = 1 reproducing mfc_siso_core bit for bit across
%        the decoupled grid, and n = 2 reproducing a hand-composed loop
%    11. mfc_fhat_decoupled_dev_block: every one of its eight settings
%        reproducing the specific block it stands in for, exactly
%
%   (the numbering matches the section comments in the body)
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
    EST = {'Algebraic (growing window)', 'Sliding window (FIR taps)'};

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
    %    must reproduce mfc_siso_core sample for sample, for BOTH controller
    %    structures. This is the twin of tests/test_composed.m, which proves
    %    the same for the stage maths.
    %
    %    Since the feedback block was removed, a composed loop puts the
    %    explicit law in a stock Discrete PID; mfc_siso.feedback stands in
    %    for it here so the comparison stays exact. With Ki = 0 the coupled
    %    case needs no feedback at all -- fb is literally Ground.
    % ---------------------------------------------------------------------
    a0 = 1; a1 = 1.4; b = 2; alpha = 2; Kp = 25; Kd = 10; Ki = 0;

    for coupled = [false true]
        if coupled, structure = STR{1}; else, structure = STR{2}; end

        mono = mfc_siso_core('model_order', ORD{2}, 'controller_structure', structure, ...
                             'estimator_type', EST{1}, 'Ts', Ts, 'alpha', alpha, ...
                             'Kp', Kp, 'Kd', Kd, 'Ki', Ki, ...
                             'ref_filter_window', 10, 'est_filter_window', 10);

        smoo = mfc_smoother_block('Ts', Ts, 'window', 10, 'output_derivatives', true);
        cmd  = mfc_command_block('alpha', alpha);
        if coupled
            % Kp and Kd live INSIDE the estimator; nothing explicit is left.
            est2 = mfc_fhat_alg2_coupled_block('Ts', Ts, 'alpha', alpha, ...
                       'Kp', Kp, 'Kd', Kd, ...
                       'est_filter_window', 10, 'est_hold_time', 0.1);
        else
            est2 = mfc_fhat_alg2_decoupled_block('Ts', Ts, 'alpha', alpha, ...
                       'est_filter_window', 10, 'est_hold_time', 0.1);
        end

        tv = 0:Ts:3;
        ym = 0; dym = 0;  yc = 0; dyc = 0;  u_km1 = 0;  worst = 0;
        e_km1 = 0; int_e = 0;
        for k = 1:numel(tv)
            r = 1.0*(tv(k) >= 0.2);

            um = step(mono, r, ym, tv(k));               % all-in-one

            [sp, ~, ddsp] = step(smoo, r);               % assembled from stages
            e = yc - sp;
            if coupled
                Fh = step(est2, e, u_km1, tv(k));        % error-driven
            else
                Fh = step(est2, yc, u_km1, tv(k));       % measurement-driven
            end
            % the stock Discrete PID an assembled loop would use
            [fb, int_e] = mfc_siso.feedback(e, e_km1, int_e, Ts, Kp, Kd, Ki, coupled);
            e_km1 = e;

            uc = step(cmd, Fh, ddsp, fb);
            u_km1 = uc;                                  % the explicit unit delay

            worst = max(worst, abs(um - uc));

            ydd = -a1*dym - a0*ym + b*um;  dym = dym + ydd*Ts;  ym = ym + dym*Ts;
            ydd = -a1*dyc - a0*yc + b*uc;  dyc = dyc + ydd*Ts;  yc = yc + dyc*Ts;
        end

        if coupled, nm = 'coupled'; else, nm = 'decoupled'; end
        check(sprintf('stage blocks reproduce mfc_siso_core, %s (max|du|=%.3g)', ...
                      nm, worst), worst < 1e-12);
    end

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

    % ---------------------------------------------------------------------
    % 7) mfc_fhat_riachy2_block (Riachy's trick): the OBJECT layer only --
    %    ports, the estimator selector, the optional Y output, reset. The
    %    numerics and the loop wiring are tests/test_riachy.m.
    %
    %    The equality asserted here is that the block is a thin wrapper:
    %    mfc_riachy_transform followed by the selected estimator, with
    %    nothing folded. If it ever starts computing something of its own,
    %    this fails.
    % ---------------------------------------------------------------------
    Kd_r = 4;  alpha_r = 2;
    for kind = {'Algebraic (growing window)', 'Sliding window (FIR)'}
        algebraic = strncmp(kind{1}, 'Algebraic', 9);
        rb = mfc_fhat_riachy2_block('estimator', kind{1}, 'Ts', Ts, ...
                 'alpha', alpha_r, 'Kd', Kd_r, ...
                 'est_filter_window', 10, 'est_hold_time', 0.1, ...
                 'window_samples', 40, 'output_Y', true);

        check(sprintf('riachy block (%s): ports are y,u_prev,t -> F_hat,Y', kind{1}), ...
              isequal(cellstr(rb.getInputNames()),  {'y'; 'u_prev'; 't'}) && ...
              isequal(cellstr(rb.getOutputNames()), {'F_hat'; 'Y'}));

        % Hand-rolled reference: the two functions the block claims to call.
        kern = mfc_siso.window_kernel(2, 40, Ts);
        rs   = struct('int_km1', 0, 'y_km1', 0);
        est  = struct('z_km1', 0, 'z_km2', 0, ...
                      'num_filt_km1', 0, 'num_filt_km2', 0, ...
                      'den_filt_km1', 0, 'den_filt_km2', 0, ...
                      'y_buf', zeros(kern.n_intervals+1, 1), ...
                      'u_buf', zeros(kern.n_intervals+1, 1));
        worst = 0;  worstY = 0;  Fb = zeros(1, 200);
        for k = 1:200
            tk = (k-1)*Ts;
            yk = 0.3*sin(5*tk) + 0.2;
            uk = cos(2*tk);
            [Fb(k), Yb] = step(rb, yk, uk, tk);

            [Yr, rs] = mfc_riachy_transform(yk, Kd_r, Ts, rs);
            if algebraic
                [Fr, est] = mfc_fhat_algebraic_second_order( ...
                    Yr, uk, alpha_r, tk, Ts, 10, 0.1, 0, 0, est);
            else
                [Fr, est] = mfc_fhat_sliding_window(Yr, uk, alpha_r, tk, kern, est);
            end
            worst  = max(worst,  abs(Fb(k) - Fr));
            worstY = max(worstY, abs(Yb - Yr));
        end
        check(sprintf('riachy block (%s): equals transform + estimator (max|dF|=%.3g)', ...
              kind{1}, worst), worst < 1e-12 && worstY < 1e-12);
        check(sprintf('riachy block (%s): F_hat is not identically zero', kind{1}), ...
              any(Fb ~= 0));

        reset(rb);
        again = zeros(1, 200);
        for k = 1:200
            tk = (k-1)*Ts;
            again(k) = step(rb, 0.3*sin(5*tk) + 0.2, cos(2*tk), tk);
        end
        check(sprintf('riachy block (%s): reset() restores the initial state', kind{1}), ...
              isequal(Fb, again));
    end

    % The selector must actually select: the two estimators disagree on the
    % same input (a guard against the dispatch silently collapsing to one).
    ra = mfc_fhat_riachy2_block('estimator', 'Algebraic (growing window)', ...
             'Ts', Ts, 'alpha', alpha_r, 'Kd', Kd_r, 'window_samples', 40);
    rw = mfc_fhat_riachy2_block('estimator', 'Sliding window (FIR)', ...
             'Ts', Ts, 'alpha', alpha_r, 'Kd', Kd_r, 'window_samples', 40);
    gap = 0;
    for k = 1:200
        tk = (k-1)*Ts;
        gap = max(gap, abs(step(ra, 0.3*sin(5*tk) + 0.2, cos(2*tk), tk) - ...
                           step(rw, 0.3*sin(5*tk) + 0.2, cos(2*tk), tk)));
    end
    check(sprintf('riachy block: the estimator selector selects (gap=%.3g)', gap), gap > 1e-6);

    % ---------------------------------------------------------------------
    % 8) mfc_fhat_riachy2_mimo_block: the vector twin, OBJECT layer only.
    %    Same thin-wrapper claim (mfc_riachy_transform with a matrix Kd,
    %    then the selected estimator), plus the two things only the vector
    %    block can get wrong: port widths and the n-vs-matrix mask checks.
    %    The numerics and the loop wiring are tests/test_riachy.m section 4.
    % ---------------------------------------------------------------------
    n_m     = 2;
    Kd_m    = [4 1; 0 6];
    alpha_m = [2 0.3; -0.1 1.5];
    for kind = {'Algebraic (growing window)', 'Sliding window (FIR)'}
        algebraic = strncmp(kind{1}, 'Algebraic', 9);
        mb = mfc_fhat_riachy2_mimo_block('estimator', kind{1}, 'n', n_m, ...
                 'Ts', Ts, 'alpha', alpha_m, 'Kd', Kd_m, ...
                 'est_filter_window', 10, 'est_hold_time', 0.1, ...
                 'window_samples', 40, 'output_Y', true);

        kern = mfc_siso.window_kernel(2, 40, Ts);
        rs   = struct('int_km1', zeros(n_m,1), 'y_km1', zeros(n_m,1));
        est  = struct('z_km1', zeros(n_m,1), 'z_km2', zeros(n_m,1), ...
                      'num_filt_km1', zeros(n_m,1), 'num_filt_km2', zeros(n_m,1), ...
                      'den_filt_km1', 0, 'den_filt_km2', 0, ...
                      'y_buf', zeros(kern.n_intervals+1, n_m), ...
                      'u_buf', zeros(kern.n_intervals+1, n_m));
        worst = 0;  worstY = 0;  Fb = zeros(n_m, 200);
        for k = 1:200
            tk = (k-1)*Ts;
            yk = [0.3*sin(5*tk) + 0.2; -0.15*sin(3*tk)];
            uk = [cos(2*tk); 0.5*cos(7*tk)];
            [Fk, Yk] = step(mb, yk, uk, tk);
            Fb(:, k) = Fk;

            [Yr, rs] = mfc_riachy_transform(yk, Kd_m, Ts, rs);
            if algebraic
                [Fr, est] = mfc_fhat_algebraic_second_order( ...
                    Yr, uk, alpha_m, tk, Ts, 10, 0.1, 0, 0, est);
            else
                [Fr, est] = mfc_fhat_sliding_window(Yr, uk, alpha_m, tk, kern, est);
            end
            worst  = max(worst,  max(abs(Fk - Fr)));
            worstY = max(worstY, max(abs(Yk - Yr)));
        end
        check(sprintf('riachy mimo (%s): F_hat and Y are n-by-1', kind{1}), ...
              isequal(size(Fb(:, end)), [n_m 1]));
        check(sprintf('riachy mimo (%s): equals transform + estimator (max|dF|=%.3g)', ...
              kind{1}, worst), worst < 1e-12 && worstY < 1e-12);
        check(sprintf('riachy mimo (%s): F_hat is not identically zero', kind{1}), ...
              any(Fb(:) ~= 0));

        reset(mb);
        again = zeros(n_m, 200);
        for k = 1:200
            tk = (k-1)*Ts;
            again(:, k) = step(mb, [0.3*sin(5*tk) + 0.2; -0.15*sin(3*tk)], ...
                                   [cos(2*tk); 0.5*cos(7*tk)], tk);
        end
        check(sprintf('riachy mimo (%s): reset() restores the initial state', kind{1}), ...
              isequal(Fb, again));
    end

    % Diagonal Kd and alpha: every channel must reduce to the SISO block.
    % (Not bit-exact for the sliding window -- the vector buffers make the
    % tap sum a matrix product, which BLAS accumulates in a different
    % order; see tests/test_riachy.m section 4b.)
    for kind = {'Algebraic (growing window)', 'Sliding window (FIR)'}
        md = mfc_fhat_riachy2_mimo_block('estimator', kind{1}, 'n', 2, 'Ts', Ts, ...
                 'alpha', alpha_r*eye(2), 'Kd', Kd_r*eye(2), 'window_samples', 40);
        sd = mfc_fhat_riachy2_block('estimator', kind{1}, 'Ts', Ts, ...
                 'alpha', alpha_r, 'Kd', Kd_r, 'window_samples', 40);
        gapd = 0;
        for k = 1:200
            tk = (k-1)*Ts;
            yk = 0.3*sin(5*tk) + 0.2;  uk = cos(2*tk);
            gapd = max(gapd, max(abs(step(md, [yk; yk], [uk; uk], tk) - step(sd, yk, uk, tk))));
        end
        check(sprintf('riachy mimo (%s): diagonal gains reduce to the SISO block (%.3g)', ...
              kind{1}, gapd), gapd < 1e-11);
    end

    % The mask checks: alpha and Kd must both be n-by-n.
    bad_alpha = false;
    try
        mbad = mfc_fhat_riachy2_mimo_block('n', 3, 'alpha', eye(2), 'Kd', eye(3));
        step(mbad, zeros(3,1), zeros(3,1), 0);
    catch
        bad_alpha = true;
    end
    bad_Kd = false;
    try
        mbad = mfc_fhat_riachy2_mimo_block('n', 2, 'alpha', eye(2), 'Kd', eye(3));
        step(mbad, zeros(2,1), zeros(2,1), 0);
    catch
        bad_Kd = true;
    end
    check('riachy mimo: a wrong-size alpha is rejected', bad_alpha);
    check('riachy mimo: a wrong-size Kd is rejected',    bad_Kd);

    % ---------------------------------------------------------------------
    % 9) The rest of the n-channel estimator surface. Each block is claimed
    %    to be a thin wrapper over one kernel call; that claim is what is
    %    tested here, together with the port widths and the first input
    %    port NAME -- which is the library's way of telling coupled from
    %    decoupled, so it is load-bearing, not cosmetic.
    % ---------------------------------------------------------------------
    n_m     = 2;
    alpha_v = [2 0.3; -0.1 1.5];
    Kp_v    = [25 3; -2 30];
    Kd_v    = [10 1; 0 12];
    yv = @(tk) [0.3*sin(5*tk) + 0.2; -0.15*sin(3*tk)];
    uv = @(tk) [cos(2*tk); 0.5*cos(7*tk)];

    % {block, first port name, reference kernel call}
    MB = { ...
      mfc_fhat_alg1_decoupled_mimo_block('n', n_m, 'Ts', Ts, 'alpha', alpha_v, ...
          'est_filter_window', 10, 'est_hold_time', 0.1), 'y', ...
      @(z, u, tk, st) mfc_fhat_algebraic_first_order(z, u, alpha_v, tk, Ts, 10, 0.1, 0, st); ...
      mfc_fhat_alg1_coupled_mimo_block('n', n_m, 'Ts', Ts, 'alpha', alpha_v, ...
          'Kp', Kp_v, 'est_filter_window', 10, 'est_hold_time', 0.1), 'err', ...
      @(z, u, tk, st) mfc_fhat_algebraic_first_order(z, u, alpha_v, tk, Ts, 10, 0.1, -Kp_v, st); ...
      mfc_fhat_alg2_coupled_mimo_block('n', n_m, 'Ts', Ts, 'alpha', alpha_v, ...
          'Kp', Kp_v, 'Kd', Kd_v, 'est_filter_window', 10, 'est_hold_time', 0.1), 'err', ...
      @(z, u, tk, st) mfc_fhat_algebraic_second_order(z, u, alpha_v, tk, Ts, 10, 0.1, -Kd_v, -Kp_v, st)};

    for b = 1:size(MB, 1)
        blk  = MB{b, 1};
        cls  = class(blk);
        est  = struct('z_km1', zeros(n_m,1), 'z_km2', zeros(n_m,1), ...
                      'num_filt_km1', zeros(n_m,1), 'num_filt_km2', zeros(n_m,1), ...
                      'den_filt_km1', 0, 'den_filt_km2', 0);
        worst = 0;  Fb = zeros(n_m, 150);
        for k = 1:150
            tk = (k-1)*Ts;
            Fb(:, k) = step(blk, yv(tk), uv(tk), tk);
            [Fr, est] = MB{b, 3}(yv(tk), uv(tk), tk, est);
            worst = max(worst, max(abs(Fb(:, k) - Fr)));
        end
        in_names = cellstr(blk.getInputNames());
        check(sprintf('%s: first input port is ''%s''', cls, MB{b, 2}), ...
              strcmp(in_names{1}, MB{b, 2}));
        check(sprintf('%s: F_hat is n-by-1', cls), isequal(size(Fb(:, end)), [n_m 1]));
        check(sprintf('%s: equals its kernel (max|dF|=%.3g)', cls, worst), worst < 1e-12);
        check(sprintf('%s: F_hat is not identically zero', cls), any(Fb(:) ~= 0));

        reset(blk);
        again = zeros(n_m, 150);
        for k = 1:150
            tk = (k-1)*Ts;
            again(:, k) = step(blk, yv(tk), uv(tk), tk);
        end
        check(sprintf('%s: reset() restores the initial state', cls), isequal(Fb, again));
    end

    % The sliding-window MIMO block, both orders, against its own kernel.
    for ord = {'First order (dot_y = F + alpha*u)', 'Second order (ddot_y = F + alpha*u)'}
        wb = mfc_fhat_window_mimo_block('n', n_m, 'Ts', Ts, 'alpha', alpha_v, ...
                 'model_order', ord{1}, 'window_samples', 20);
        if strncmp(ord{1}, 'First', 5), o = 1; else, o = 2; end
        kern = mfc_siso.window_kernel(o, 20, Ts);
        est  = struct('y_buf', zeros(21, n_m), 'u_buf', zeros(21, n_m));
        worst = 0;  Fb = zeros(n_m, 150);
        for k = 1:150
            tk = (k-1)*Ts;
            Fb(:, k) = step(wb, yv(tk), uv(tk), tk);
            [Fr, est] = mfc_fhat_sliding_window(yv(tk), uv(tk), alpha_v, tk, kern, est);
            worst = max(worst, max(abs(Fb(:, k) - Fr)));
        end
        check(sprintf('window mimo (%d): equals mfc_fhat_sliding_window (max|dF|=%.3g)', o, worst), ...
              worst < 1e-12);
        check(sprintf('window mimo (%d): held at zeros(n,1) while the window fills', o), ...
              isequal(Fb(:, 1:20), zeros(n_m, 20)));
    end

    % The n-channel smoother: channel-wise, and equal to mfc_siso.ref_traj.
    sm = mfc_smoother_mimo_block('n', n_m, 'Ts', Ts, 'window', 10, 'output_derivatives', true);
    x1 = 0;  x2 = 0;  worst = 0;
    for k = 1:100
        tk = (k-1)*Ts;
        xk = yv(tk);
        [xf, dxf, ddxf] = step(sm, xk);
        [xr, dxr, ddxr] = mfc_siso.ref_traj(xk, x1, x2, Ts, 10, true);
        x2 = x1;  x1 = xr;
        worst = max([worst, max(abs(xf - xr)), max(abs(dxf - dxr)), max(abs(ddxf - ddxr))]);
    end
    check(sprintf('smoother mimo: equals mfc_siso.ref_traj (max|d|=%.3g)', worst), worst < 1e-12);
    check('smoother mimo: outputs are n-by-1', isequal(size(xf), [n_m 1]));

    % ---------------------------------------------------------------------
    % 10) mfc_mimo_core. Two claims, both in its help:
    %     (a) n = 1 takes the scalar code path, so it must reproduce
    %         mfc_siso_core BIT FOR BIT -- not to a tolerance;
    %     (b) at n = 2 it is the assembled form of the same stages the
    %         individual blocks call, so a hand-composed loop must match it.
    % ---------------------------------------------------------------------
    for ord = {'First order (dot_y = F + alpha*u)', 'Second order (ddot_y = F + alpha*u)'}
        for esti = {'Algebraic (growing window)', 'Sliding window (FIR taps)'}
            args = {'model_order', ord{1}, 'estimator_type', esti{1}, ...
                    'controller_structure', 'Decoupled (measurement-driven estimator, explicit iP/iPD)', ...
                    'Ts', Ts, 'Kp', 25, 'Kd', 10, 'Ki', 0.5, ...
                    'est_filter_window', 10, 'est_hold_time', 0.1};
            c1 = mfc_siso_core(args{:}, 'alpha', 1);
            cm = mfc_mimo_core(args{:}, 'n', 1, 'alpha', 1);
            same = true;
            for k = 1:200
                tk = (k-1)*Ts;
                [u1, F1] = step(c1, 1, 0.3*sin(5*tk), tk);
                [um, Fm] = step(cm, 1, 0.3*sin(5*tk), tk);
                same = same && isequal(u1, um) && isequal(F1, Fm);
            end
            check(sprintf('mimo core (n=1, %s / %s): identical to mfc_siso_core', ...
                  ord{1}(1:5), esti{1}(1:5)), same);
        end
    end

    % (b) n = 2, decoupled 2nd order algebraic, composed by hand out of the
    % smoother, the estimator block, mfc_siso.feedback (standing in for the
    % stock PID, so the comparison stays bit-exact) and the command block.
    core = mfc_mimo_core('n', n_m, 'Ts', Ts, 'alpha', alpha_v, ...
               'model_order', 'Second order (ddot_y = F + alpha*u)', ...
               'controller_structure', 'Decoupled (measurement-driven estimator, explicit iP/iPD)', ...
               'estimator_type', 'Algebraic (growing window)', ...
               'Kp', 25, 'Kd', 10, 'Ki', 0, 'ref_filter_window', 10, ...
               'est_filter_window', 10, 'est_hold_time', 0.1);
    smo = mfc_smoother_mimo_block('n', n_m, 'Ts', Ts, 'window', 10, 'output_derivatives', true);
    fes = mfc_fhat_alg2_decoupled_mimo_block('n', n_m, 'Ts', Ts, 'alpha', alpha_v, ...
               'est_filter_window', 10, 'est_hold_time', 0.1);
    cmd = mfc_command_mimo_block('n', n_m, 'alpha', alpha_v);

    sp_ref = [1; 0.5];
    err_km1 = zeros(n_m, 1);  int_err = zeros(n_m, 1);  u_prev = zeros(n_m, 1);
    y = zeros(n_m, 1);  dy = zeros(n_m, 1);
    worst = 0;
    for k = 1:200
        tk = (k-1)*Ts;

        u_core = step(core, sp_ref, y, tk);

        [sp_filt, ~, ddot_sp] = step(smo, sp_ref);
        F_hat = step(fes, y, u_prev, tk);
        err   = y - sp_filt;
        [fb, int_err] = mfc_siso.feedback(err, err_km1, int_err, Ts, 25, 10, 0, false);
        err_km1 = err;
        u_comp  = step(cmd, F_hat, ddot_sp, fb);
        u_prev  = u_comp;

        worst = max(worst, max(abs(u_core - u_comp)));

        % one shared plant, driven by the core's command
        ddy = -[9.81; 9.81] - 0.2*dy + [1 0.35; -0.2 1.2]*u_core;
        dy  = dy + Ts*ddy;
        y   = y + Ts*dy;
    end
    check(sprintf('mimo core: equals the hand-composed loop (max|du|=%.3g)', worst), ...
          worst < 1e-12);

    % ---------------------------------------------------------------------
    % 11) mfc_fhat_decoupled_dev_block. Its whole claim is that each mask
    %     setting IS the specific block it replaces, so that is what is
    %     tested: all eight settings (2 orders x 2 estimators x {n=1, n=2})
    %     against the block a finished model would ship instead. Bit-exact
    %     is demanded, not a tolerance -- the dev block calls the same
    %     kernel with the same arguments, so anything else is a bug.
    % ---------------------------------------------------------------------
    ORD2 = {'First order (dot_y = F + alpha*u)', 'Second order (ddot_y = F + alpha*u)'};
    a_s  = 2;                       % scalar alpha
    a_v  = [2 0.3; -0.1 1.5];       % matrix alpha

    for oi = 1:2
        for kind = {'Algebraic (growing window)', 'Sliding window (FIR taps)'}
            for width = [1 2]
                if width == 1, al = a_s; else, al = a_v; end
                algebraic = strncmp(kind{1}, 'Algebraic', 9);

                dev = mfc_fhat_decoupled_dev_block('n', width, 'Ts', Ts, 'alpha', al, ...
                          'model_order', ORD2{oi}, 'estimator', kind{1}, ...
                          'est_filter_window', 10, 'est_hold_time', 0.1, ...
                          'window_samples', 20);

                % the block a finished model would use instead
                if algebraic && oi == 1 && width == 1
                    ref_blk = mfc_fhat_alg1_decoupled_block('Ts', Ts, 'alpha', al, ...
                                  'est_filter_window', 10, 'est_hold_time', 0.1);
                elseif algebraic && oi == 2 && width == 1
                    ref_blk = mfc_fhat_alg2_decoupled_block('Ts', Ts, 'alpha', al, ...
                                  'est_filter_window', 10, 'est_hold_time', 0.1);
                elseif algebraic && oi == 1
                    ref_blk = mfc_fhat_alg1_decoupled_mimo_block('n', width, 'Ts', Ts, ...
                                  'alpha', al, 'est_filter_window', 10, 'est_hold_time', 0.1);
                elseif algebraic
                    ref_blk = mfc_fhat_alg2_decoupled_mimo_block('n', width, 'Ts', Ts, ...
                                  'alpha', al, 'est_filter_window', 10, 'est_hold_time', 0.1);
                elseif width == 1
                    ref_blk = mfc_fhat_window_block('Ts', Ts, 'alpha', al, ...
                                  'model_order', ORD2{oi}, 'window_samples', 20);
                else
                    ref_blk = mfc_fhat_window_mimo_block('n', width, 'Ts', Ts, 'alpha', al, ...
                                  'model_order', ORD2{oi}, 'window_samples', 20);
                end

                same = true;  nonzero = false;
                for k = 1:200
                    tk = (k-1)*Ts;
                    if width == 1
                        yk = 0.3*sin(5*tk) + 0.2;   uk = cos(2*tk);
                    else
                        yk = yv(tk);                uk = uv(tk);
                    end
                    Fd = step(dev, yk, uk, tk);
                    Fr = step(ref_blk, yk, uk, tk);
                    same    = same && isequal(Fd, Fr);
                    nonzero = nonzero || any(Fd ~= 0);
                end
                check(sprintf('dev bench (%s / %s / n=%d): identical to %s', ...
                      ORD2{oi}(1:5), kind{1}(1:5), width, class(ref_blk)), same);
                check(sprintf('dev bench (%s / %s / n=%d): F_hat is not identically zero', ...
                      ORD2{oi}(1:5), kind{1}(1:5), width), nonzero);
            end
        end
    end

    % Ports must not move when the dropdowns do -- that is the whole point.
    d1 = mfc_fhat_decoupled_dev_block('n', 2, 'Ts', Ts, 'alpha', a_v, ...
             'model_order', ORD2{1}, 'estimator', 'Sliding window (FIR taps)');
    d2 = mfc_fhat_decoupled_dev_block('n', 2, 'Ts', Ts, 'alpha', a_v, ...
             'model_order', ORD2{2}, 'estimator', 'Algebraic (growing window)');
    check('dev bench: ports are y,u_prev,t -> F_hat whatever the variant', ...
          isequal(cellstr(d1.getInputNames()),  {'y'; 'u_prev'; 't'}) && ...
          isequal(cellstr(d2.getInputNames()),  {'y'; 'u_prev'; 't'}) && ...
          isequal(cellstr(d1.getOutputNames()), {'F_hat'}) && ...
          isequal(cellstr(d2.getOutputNames()), {'F_hat'}));

    bad_dev = false;
    try
        dbad = mfc_fhat_decoupled_dev_block('n', 3, 'alpha', eye(2));
        step(dbad, zeros(3,1), zeros(3,1), 0);
    catch
        bad_dev = true;
    end
    check('dev bench: a wrong-size alpha is rejected', bad_dev);

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
