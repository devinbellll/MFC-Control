function [traces, names, t] = mfc_golden_trace_mimo()
%MFC_GOLDEN_TRACE_MIMO Reference closed-loop traces for the 2x2 MIMO pair.
%
%   [traces, names, t] = MFC_GOLDEN_TRACE_MIMO()
%
%   Drives mfc_fhat_alg2_decoupled_mimo_block's math (2nd order, decoupled,
%   matrix alpha) and mfc_command_mimo_block's math through a closed loop,
%   for two variants that each pin a property the firmware port depends on:
%
%     'mimo_diag'  alpha = eye(2), two IDENTICAL channels driven by the same
%                  reference and the same plant as tests/golden/2nd_decoupled_alg.csv.
%                  With a diagonal alpha the matrix solve in mfc_command_mimo
%                  and the alpha*d2u term in the estimator must reduce to
%                  the plain scalar path, so channel 1 and channel 2 must
%                  each reproduce that SISO trace -- a direct cross-check
%                  against an existing golden trace, not only a capture of
%                  itself. See TEST_GOLDEN_MIMO.
%
%     'mimo_cross' alpha has off-diagonal terms and the two channels drive
%                  DIFFERENT plants, so F_hat is genuinely different per
%                  channel (num_raw is per-element) while the estimator's
%                  den_raw = t^2 stays scalar and shared across the vector
%                  -- that asymmetry is the whole reason one integration
%                  window can serve a vector channel, and it is the thing
%                  most likely to be got wrong in a C port. The trace has a
%                  single 'valid' column (not one per channel), which is
%                  itself evidence the denominator is shared.
%
%   Runs in Octave as well as MATLAB: it calls the plain functions directly
%   and never touches matlab.System, mirroring MFC_GOLDEN_TRACE and
%   tests/test_composed.m generalized to n-by-1 vector signals.
%
%   Outputs
%     traces : {1 x 2} cell, each [N x (5n+1)] =
%              [u(1:n), F_hat(1:n), sp_filt(1:n), err(1:n), u_raw(1:n), valid]
%     names  : {1 x 2} cell of variant labels (also the csv basenames)
%     t      : [N x 1] time vector [s]
%
%   See also GOLDEN_CAPTURE_MIMO, TEST_GOLDEN_MIMO, MFC_GOLDEN_TRACE,
%   MFC_FHAT_ALGEBRAIC_SECOND_ORDER, MFC_SISO.COMMAND_MIMO.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'functions'));

% --- shared tuning (identical to mfc_golden_trace's 2nd-order variants) --
Ts        = 0.01;
WFilter   = 10;
FFilter   = 10;
hold_time = 0.1;
Ki        = 0;

tau_design = 0.05;
p   = (1/tau_design) / 5;          % 4 rad/s closed-loop target
Kp2 = p^2;  Kd2 = 2*p;

t   = (0:Ts:6)';
ref = double(t >= 0.2);
n   = 2;

% --- Variant 1: diagonal alpha, two identical channels --------------------
plant_diag      = struct('g', [9.81; 9.81], 'd', [0.2; 0.2], ...
                          'b', [1; 1], 'tau', [0.05; 0.05]);
alpha_diag      = eye(n);
T1 = run_mimo_variant(n, alpha_diag, plant_diag, Ts, Kp2, Kd2, Ki, ...
                       WFilter, FFilter, hold_time, t, ref);

% --- Variant 2: off-diagonal alpha, two different plants ------------------
plant_cross     = struct('g', [9.81; 9.81], 'd', [0.2; 0.35], ...
                          'b', [1; 1.4], 'tau', [0.05; 0.08]);
alpha_cross     = [1.00 0.35; -0.20 1.20];
T2 = run_mimo_variant(n, alpha_cross, plant_cross, Ts, Kp2, Kd2, Ki, ...
                       WFilter, FFilter, hold_time, t, ref);

traces = {T1, T2};
names  = {'mimo_diag', 'mimo_cross'};
end


function T = run_mimo_variant(n, alpha, plant, Ts, Kp, Kd, Ki, WFilter, ...
                               FFilter, hold_time, t, ref)
%RUN_MIMO_VARIANT Decoupled 2nd-order MIMO loop, wired by hand.
%
%   smoother -> mfc_fhat_algebraic_second_order (vector z, matrix alpha) ->
%   explicit per-channel PID (mfc_siso.feedback, decoupled) ->
%   mfc_siso.command_mimo. Same wiring tests/test_composed.m uses for the
%   scalar case, generalized to n-by-1 vectors; mfc_command_mimo_block has
%   no output filter or saturation by default, so u = u_raw directly.

N = numel(t);
sp1 = zeros(n, 1);  sp2 = zeros(n, 1);
est = struct('z_km1', zeros(n, 1), 'z_km2', zeros(n, 1), ...
             'num_filt_km1', zeros(n, 1), 'num_filt_km2', zeros(n, 1), ...
             'den_filt_km1', 0, 'den_filt_km2', 0);   % den is scalar: shared t^2
err_km1 = zeros(n, 1);  int_err = zeros(n, 1);
u_km1   = zeros(n, 1);

y = zeros(n, 1);  dy = zeros(n, 1);  u_act = zeros(n, 1);
T = zeros(N, 5*n + 1);

for k = 1:N
    u_prev   = u_km1;                              % the loop's unit delay
    setpoint = ref(k) * ones(n, 1);                 % same reference, every channel

    % 1) input smoother (vectorized: plain arithmetic, no per-channel loop needed)
    [sp_filt, ~, ddot_sp] = mfc_siso.ref_traj(setpoint, sp1, sp2, Ts, WFilter, true);
    err = y - sp_filt;

    % 2) F-hat estimator: measurement-driven (decoupled), matrix alpha
    [F_hat, est, dbg] = mfc_fhat_algebraic_second_order( ...
        y, u_prev, alpha, t(k), Ts, FFilter, hold_time, 0, 0, est);

    % 3) explicit per-channel feedback (decoupled: fb = Kd*dot_err + Kp*err + Ki*int_err)
    int_err_prev = int_err;
    [fb, int_err] = mfc_siso.feedback(err, err_km1, int_err_prev, Ts, Kp, Kd, Ki, false);
    err_km1 = err;

    % 4) command: the matrix solve
    ff    = ddot_sp;
    u_raw = mfc_siso.command_mimo(F_hat, ff, fb, alpha);
    u     = u_raw;                                  % no EMA filter, no saturation
    u_km1 = u;

    T(k, :) = [u.', F_hat.', sp_filt.', err.', u_raw.', double(dbg.valid)];

    sp2 = sp1;  sp1 = sp_filt;                      % advance the smoother state

    % plant: per-channel double integrator + gravity + drag + actuator lag
    % (physically decoupled -- any cross-channel coupling in this trace
    % comes only from alpha, inside the estimator and the command law)
    u_act = u_act + (u - u_act) ./ plant.tau * Ts;
    ddy   = -plant.g - plant.d .* dy + plant.b .* u_act;
    dy    = dy + ddy * Ts;
    y     = y + dy * Ts;
end
end
