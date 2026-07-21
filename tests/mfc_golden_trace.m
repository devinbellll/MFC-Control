function [traces, names, t] = mfc_golden_trace()
%MFC_GOLDEN_TRACE Reference closed-loop traces for every supported MFC variant.
%
%   [traces, names, t] = MFC_GOLDEN_TRACE()
%
%   Drives the controller core over all six supported variants on the same
%   plant used by examples/val_mfc.m (drone z-axis: double integrator +
%   gravity + linear drag + first-order actuator lag) and returns the raw
%   per-sample outputs. Runs in Octave as well as MATLAB: it calls the plain
%   controller functions directly and never touches matlab.System.
%
%   This is the no-regression contract for the decomposition refactor.
%   GOLDEN_CAPTURE writes these traces to tests/golden/*.csv; TEST_GOLDEN
%   re-runs this function and demands an EXACT match against those files.
%
%   Outputs
%     traces : {1 x 6} cell, each [N x 6] = [u, F_hat, sp_filt, err, u_raw, valid]
%     names  : {1 x 6} cell of short variant labels (also the csv basenames)
%     t      : [N x 1] time vector [s]
%
%   NOTE: the '1st_coupled_alg' variant is EXPECTED to diverge on this plant
%   (a true double integrator cannot be stabilized from a 1st-order model
%   with no derivative room to fold -- see val_mfc.m). Its divergence is
%   itself part of the contract: it must diverge the same way afterwards.
%
%   See also GOLDEN_CAPTURE, TEST_GOLDEN.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'functions'));

% --- plant (identical to examples/val_mfc.m) -----------------------------
g           = 9.81;    % gravity [m/s^2], the constant disturbance F must absorb
d           = 0.2;     % linear aerodynamic drag [1/s]
alpha_plant = 1;       % thrust-accel effectiveness
tau         = 0.05;    % actuator lag [s]

% --- shared tuning -------------------------------------------------------
Ts      = 0.01;
alpha   = alpha_plant;
WFilter = 10;
FFilter = 10;
Ki      = 0;

p   = (1/tau) / 5;          % 4 rad/s closed-loop target
Kp2 = p^2;   Kd2 = 2*p;
Kp1 = p;     Kd1 = p;

% --- variant grid: {order, structure, estimator, name, Kp, Kd} -----------
V = { ...
    2, 'coupled',   'algebraic',      '2nd_coupled_alg',   Kp2, Kd2; ...
    2, 'decoupled', 'algebraic',      '2nd_decoupled_alg', Kp2, Kd2; ...
    2, 'decoupled', 'sliding_window', '2nd_decoupled_win', Kp2, Kd2; ...
    1, 'coupled',   'algebraic',      '1st_coupled_alg',   Kp1, Kd1; ...
    1, 'decoupled', 'algebraic',      '1st_decoupled_alg', Kp1, Kd1; ...
    1, 'decoupled', 'sliding_window', '1st_decoupled_win', Kp1, Kd1};
nc = size(V, 1);

t   = (0:Ts:6)';
N   = numel(t);
ref = double(t >= 0.2);

traces = cell(1, nc);
names  = V(:, 4)';

for c = 1:nc
    cfg = mfc_siso.config( ...
        'model_order',       V{c, 1}, ...
        'structure',         V{c, 2}, ...
        'estimator',         V{c, 3}, ...
        'Ts', Ts, 'alpha', alpha, 'Kp', V{c, 5}, 'Kd', V{c, 6}, 'Ki', Ki, ...
        'ref_filter_window', WFilter, 'est_filter_window', FFilter);
    state = mfc_siso.init(cfg);

    y = 0; dy = 0; u_act = 0;
    T = zeros(N, 6);
    for k = 1:N
        [out, state] = mfc_siso.step(ref(k), y, t(k), state.u_km1, alpha, cfg, state);
        T(k, :) = [out.u, out.F_hat, out.sp_filt, out.err, out.u_raw, double(out.est_valid)];

        u_act = u_act + (out.u - u_act)/tau * Ts;                % actuator lag
        dy    = dy    + (-g - d*dy + alpha_plant*u_act)*Ts;      % gravity + drag + thrust
        y     = y     + dy*Ts;
    end
    traces{c} = T;
end
end
