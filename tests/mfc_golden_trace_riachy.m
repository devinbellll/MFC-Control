function [traces, names, t] = mfc_golden_trace_riachy()
%MFC_GOLDEN_TRACE_RIACHY Reference closed-loop traces for Riachy's trick.
%
%   [traces, names, t] = MFC_GOLDEN_TRACE_RIACHY()
%
%   Riachy's trick was the one estimator family with no no-regression
%   contract: tests/test_riachy.m checks that it BEHAVES (the transform is a
%   trapezoid, the estimate carries Kd*dot_y, the loop tracks), but nothing
%   pinned the actual numbers. This file is that contract, for the same
%   reason the other golden traces exist: a firmware port, or a change to
%   the transform or the taps, has to reproduce these samples exactly.
%
%   Four variants, each closing the loop the way mfc_fhat_riachy2_block's
%   help says to -- F_hat from the Riachy estimator, feedback with D = 0,
%   and ff = ddot_sp + Kd*dot_sp:
%
%     'riachy_win'       SISO, sliding-window estimator on Y. The
%                        recommended pairing: finite memory, so the
%                        unbounded int y never accumulates weight.
%     'riachy_alg'       SISO, algebraic estimator on Y. Kept precisely
%                        because it is the awkward one -- the growing t^2
%                        weights see the ramping Y -- so a change in that
%                        behaviour shows up here rather than in someone's
%                        model.
%     'riachy_mimo_win'  n = 2, matrix Kd and matrix alpha, a plant whose
%                        input gain B is a genuine matrix and alpha = B.
%                        Kd is NOT diagonal, so the trace pins cross-channel
%                        derivative folding -- the thing the NxN port buys
%                        over two SISO blocks.
%     'riachy_mimo_alg'  the same, algebraic.
%
%   Runs in Octave as well as MATLAB: plain functions only, never
%   matlab.System.
%
%   Outputs
%     traces : {1 x 4} cell, each [N x (5n+1)] =
%              [u(1:n), F_hat(1:n), Y(1:n), sp_filt(1:n), err(1:n), valid]
%              Y is logged instead of u_raw: there is no output filter or
%              saturation here (u_raw == u), and Y is the one signal unique
%              to this estimator.
%     names  : {1 x 4} cell of variant labels (also the csv basenames)
%     t      : [N x 1] time vector [s]
%
%   See also GOLDEN_CAPTURE_RIACHY, TEST_GOLDEN_RIACHY, TEST_RIACHY,
%   MFC_RIACHY_TRANSFORM, MFC_GOLDEN_TRACE_MIMO.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'functions'));

% --- shared tuning (the drone plant of MFC_GOLDEN_TRACE) -----------------
Ts        = 0.01;
WFilter   = 10;
FFilter   = 10;                    % algebraic num/den smoother, window variants
FFilter_a = 5;                     % algebraic variants: HALF that, for the same
                                   % reason the window is short -- Fk moves fast,
                                   % and 10 samples of smoothing on a growing
                                   % window lags it into instability here
hold_time = 0.1;
Ki        = 0;
Twin      = 10;                    % window intervals: Tw = 0.1 s -- short on
                                   % purpose: Fk carries Kd*dot_y and moves
                                   % fast, so a long window lags it into
                                   % instability (0.4 s diverges on this plant)

tau  = 0.05;
p    = (1/tau) / 5;                % 4 rad/s closed-loop target
Kp   = p^2;   Kd = 2*p;

t   = (0:Ts:6)';
ref = double(t >= 0.2);

siso_plant = struct('g', 9.81, 'd', 0.2, 'B', 1, 'tau', tau);

T1 = run_riachy_loop(1, 1, Kd, Kp*1, Ki*1, siso_plant, Ts, t, ref, ...
                     false, WFilter, FFilter, hold_time, Twin);
T2 = run_riachy_loop(1, 1, Kd, Kp*1, Ki*1, siso_plant, Ts, t, ref, ...
                     true,  WFilter, FFilter_a, hold_time, Twin);

% --- n = 2: matrix B, alpha = B, and a NON-diagonal Kd -------------------
% Identical per-channel dynamics, all the coupling in B and in Kd, so the
% trace isolates the matrix folding rather than a tuning mismatch.
n       = 2;
B_cross = [1.00 0.35; -0.20 1.20];
Kd_m    = [Kd 1.5; -0.8 Kd];
mimo_plant = struct('g', [9.81; 9.81], 'd', [0.2; 0.2], ...
                    'B', B_cross, 'tau', [tau; tau]);

T3 = run_riachy_loop(n, B_cross, Kd_m, Kp*eye(n), Ki*eye(n), mimo_plant, Ts, t, ref, ...
                     false, WFilter, FFilter, hold_time, Twin);
T4 = run_riachy_loop(n, B_cross, Kd_m, Kp*eye(n), Ki*eye(n), mimo_plant, Ts, t, ref, ...
                     true,  WFilter, FFilter_a, hold_time, Twin);

traces = {T1, T2, T3, T4};
names  = {'riachy_win', 'riachy_alg', 'riachy_mimo_win', 'riachy_mimo_alg'};
end


function T = run_riachy_loop(n, alpha, Kd, Kp, Ki, plant, Ts, t, ref, ...
                              algebraic, WFilter, FFilter, hold_time, Twin)
%RUN_RIACHY_LOOP One Riachy loop, wired exactly as the block help states.
%
%   smoother -> mfc_riachy_transform -> 2nd-order estimator on Y ->
%   feedback with Kd = 0 (P/PI only) -> command. n = 1 takes the scalar
%   command (a division); n > 1 takes mfc_siso.command_mimo (a solve).

kern = mfc_siso.window_kernel(2, Twin, Ts);
rs   = struct('int_km1', zeros(n,1), 'y_km1', zeros(n,1));
est  = struct('z_km1', zeros(n,1), 'z_km2', zeros(n,1), ...
              'num_filt_km1', zeros(n,1), 'num_filt_km2', zeros(n,1), ...
              'den_filt_km1', 0, 'den_filt_km2', 0, ...
              'y_buf', zeros(kern.n_intervals+1, n), ...
              'u_buf', zeros(kern.n_intervals+1, n));

sp1 = zeros(n,1);  sp2 = zeros(n,1);
err_km1 = zeros(n,1);  int_err = zeros(n,1);  u_prev = zeros(n,1);
y = zeros(n,1);  dy = zeros(n,1);  u_act = zeros(n,1);

N = numel(t);
T = zeros(N, 5*n + 1);
for k = 1:N
    [sp_filt, dot_sp, ddot_sp] = mfc_siso.ref_traj( ...
        ref(k)*ones(n,1), sp1, sp2, Ts, WFilter, true);
    err = y - sp_filt;

    % Riachy: the auxiliary output, then a standard 2nd-order estimator
    [Y, rs] = mfc_riachy_transform(y, Kd, Ts, rs);
    if algebraic
        [F_hat, est, dbg] = mfc_fhat_algebraic_second_order( ...
            Y, u_prev, alpha, t(k), Ts, FFilter, hold_time, 0, 0, est);
    else
        [F_hat, est, dbg] = mfc_fhat_sliding_window(Y, u_prev, alpha, t(k), kern, est);
    end

    % Feedback with NO derivative term: Kd lives in Y and in ff
    [fb, int_err] = mfc_siso.feedback(err, err_km1, int_err, Ts, Kp, 0*Kp, Ki, false);
    err_km1 = err;

    ff = ddot_sp + Kd*dot_sp;
    if n > 1
        u = mfc_siso.command_mimo(F_hat, ff, fb, alpha);
    else
        u = mfc_siso.command(F_hat, ff, fb, alpha);
    end
    u_prev = u;

    T(k, :) = [u.', F_hat.', Y.', sp_filt.', err.', double(dbg.valid)];

    sp2 = sp1;  sp1 = sp_filt;

    u_act = u_act + (u - u_act) ./ plant.tau * Ts;
    ddy   = -plant.g - plant.d .* dy + plant.B * u_act;
    dy    = dy + ddy * Ts;
    y     = y + dy * Ts;
end
end
