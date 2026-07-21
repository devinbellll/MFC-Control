function test_composed()
%TEST_COMPOSED Assert a hand-assembled stage pipeline equals mfc_siso.step.
%
%   octave --no-gui -q tests/test_composed.m     (also runs in MATLAB)
%
%   Rebuilds the controller from the individual pipeline stages --
%   mfc_siso.ref_traj, an mfc_fhat_* estimator, mfc_siso.feedback,
%   mfc_siso.command, mfc_siso.limit -- wired together by hand, and checks
%   it against the golden traces of the all-in-one mfc_siso.step.
%
%   This is the proof that the decomposition is FAITHFUL: the stage blocks
%   in blocks/ are thin wrappers over exactly these calls, so if this
%   pipeline matches, a Simulink model assembled from those blocks computes
%   what mfc_siso_core computes. It deliberately duplicates the wiring
%   rather than calling mfc_siso.step -- duplicating it is the whole point.
%
%   Scope: this validates the STAGE MATH, in Octave, with no Simulink. It
%   does not validate Simulink port wiring or sample-time propagation; run
%   examples/val_mfc_composed.m in MATLAB for that.
%
%   NOTE: saturation is off in the golden configuration, so the one-sample
%   anti-windup difference between a composed loop and mfc_siso_core (see
%   Knowledge/block-library-signal-flow.md) is not exercised here and the
%   match is exact. With saturation on, expect the composed loop to differ
%   by one sample of integrator freeze.
%
%   See also MFC_GOLDEN_TRACE, TEST_GOLDEN, MFC_SISO.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'functions'));
golden = fullfile(here, 'golden');

% Same plant and tuning as mfc_golden_trace (kept in step with it by hand;
% the golden csv files are the shared contract).
g = 9.81; d = 0.2; alpha_plant = 1; tau = 0.05;
Ts = 0.01; alpha = alpha_plant; WFilter = 10; FFilter = 10; Ki = 0;
hold_time = 0.1;
p = (1/tau)/5;
Kp2 = p^2; Kd2 = 2*p; Kp1 = p; Kd1 = p;

V = { ...
    2, true,  true,  '2nd_coupled_alg',   Kp2, Kd2; ...
    2, false, true,  '2nd_decoupled_alg', Kp2, Kd2; ...
    2, false, false, '2nd_decoupled_win', Kp2, Kd2; ...
    1, true,  true,  '1st_coupled_alg',   Kp1, Kd1; ...
    1, false, true,  '1st_decoupled_alg', Kp1, Kd1; ...
    1, false, false, '1st_decoupled_win', Kp1, Kd1};

t = (0:Ts:6)';  N = numel(t);  ref = double(t >= 0.2);
n_fail = 0;

for c = 1:size(V, 1)
    order = V{c, 1};  coupled = V{c, 2};  algebraic = V{c, 3};
    name  = V{c, 4};  Kp = V{c, 5};  Kd = V{c, 6};

    % --- stage-local state, exactly as each block would own it ----------
    sp1 = 0; sp2 = 0;                                   % smoother block
    est = struct('z_km1', 0, 'z_km2', 0, ...            % estimator block
                 'num_filt_km1', 0, 'num_filt_km2', 0, ...
                 'den_filt_km1', 0, 'den_filt_km2', 0);
    if ~algebraic
        kernel = mfc_siso.window_kernel(order, FFilter, Ts);
        est.y_buf = zeros(kernel.n_intervals + 1, 1);
        est.u_buf = zeros(kernel.n_intervals + 1, 1);
    end
    err_km1 = 0; int_err = 0;                           % feedback block
    u_km1 = 0;                                          % command filter block

    y = 0; dy = 0; u_act = 0;
    T = zeros(N, 6);

    for k = 1:N
        u_prev = u_km1;                                 % the unit delay in the loop

        % 1) input smoother
        [sp_filt, dot_sp, ddot_sp] = mfc_siso.ref_traj(ref(k), sp1, sp2, Ts, WFilter, true);
        err = y - sp_filt;

        % 2) F-hat estimator
        if algebraic
            if coupled, z = err; else, z = y; end
            if order == 2
                if coupled, a_f = -Kd; b_f = -Kp; else, a_f = 0; b_f = 0; end
                [F_hat, est, dbg] = mfc_fhat_algebraic_second_order( ...
                    z, u_prev, alpha, t(k), Ts, FFilter, hold_time, a_f, b_f, est);
            else
                if coupled, b_f = -Kp; else, b_f = 0; end
                [F_hat, est, dbg] = mfc_fhat_algebraic_first_order( ...
                    z, u_prev, alpha, t(k), Ts, FFilter, hold_time, b_f, est);
            end
        else
            [F_hat, est, dbg] = mfc_fhat_sliding_window(y, u_prev, alpha, t(k), kernel, est);
        end

        % 3) feedback law
        int_err_prev = int_err;
        [fb, int_err] = mfc_siso.feedback(err, err_km1, int_err_prev, Ts, Kp, Kd, Ki, coupled);

        % 4) command
        if order == 2, ff = ddot_sp; else, ff = dot_sp; end
        u_raw = mfc_siso.command(F_hat, ff, fb, alpha);

        % 5) command filter + saturation (off here)
        [u, frozen] = mfc_siso.limit(u_raw, u_prev, 1, false, -600, 600);
        if frozen, int_err = int_err_prev; end

        T(k, :) = [u, F_hat, sp_filt, err, u_raw, double(dbg.valid)];

        % advance the stage-local states
        sp2 = sp1;  sp1 = sp_filt;
        err_km1 = err;
        u_km1 = u;

        % plant
        u_act = u_act + (u - u_act)/tau * Ts;
        dy    = dy    + (-g - d*dy + alpha_plant*u_act)*Ts;
        y     = y     + dy*Ts;
    end

    f = fullfile(golden, [name '.csv']);
    refT = dlmread(f, ',', 1, 0);
    if isequaln(refT, T)
        fprintf('  PASS  %-20s  composed pipeline == mfc_siso.step\n', name);
    else
        n_fail = n_fail + 1;
        fprintf('  FAIL  %-20s  max|d| = %.6g\n', name, max(max(abs(refT - T))));
    end
end

fprintf('\ntest_composed: %d/%d variants match\n', size(V, 1) - n_fail, size(V, 1));
if n_fail > 0
    error('test_composed: %d variant(s) do not decompose faithfully', n_fail);
end
end
