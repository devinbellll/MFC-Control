clear

%% check_Darko_long.m
%
% Verify the longitudinal-only Darko variant (sysEqt_long.m +
% Darko_long.slx) against the full 6-DOF phi-theory model:
%
%   1. trim residual : at a longTrim point, sysEqt_long derivative ~ 0
%   2. lateral invariance : lateral derivative channels are exactly zero
%      for arbitrary (even corrupted) states and inputs
%   3. longitudinal equivalence : on pure-longitudinal states with
%      symmetric inputs, sysEqt_long matches sysEqt channel by channel
%   4. sim smoke test : Darko_long.slx runs from trim; lateral states
%      stay identically zero
%
% State channels: [vN vE vD p q r q0 qx qy qz]; longitudinal = [1 3 5 7 9],
% lateral = [2 4 6 8 10]. Input of sysEqt_long: u = [w_prop; delta].

drone = generateDarko();
lon = [1 3 5 7 9];
lat = [2 4 6 8 10];
rng(42);

%% --- 1. trim residual ---
theta_t = 90*pi/180;
[xts, uts, N] = longTrim(theta_t, drone);
assert(N >= 1, 'longTrim found no trim point');
xt = xts(:,1);
ut_long = [abs(uts(1,1)); uts(3,1)];

res = sysEqt_long(xt, ut_long, zeros(3,1), drone);
err1 = norm(res);
fprintf('1. trim residual |dxdt| at theta = %g deg   : %.3g\n', theta_t*180/pi, err1);

%% --- 2. lateral invariance for arbitrary states/inputs ---
err2 = 0;
for k = 1:20
    x = randn(10,1);                  % deliberately non-longitudinal state
    x(7:10) = x(7:10)/norm(x(7:10));
    u = [800*rand; 0.5*randn];
    w = randn(3,1);
    dxdt = sysEqt_long(x, u, w, drone);
    err2 = max(err2, max(abs(dxdt(lat))));
end
fprintf('2. max lateral derivative (random x, u, w)  : %.3g\n', err2);

%% --- 3. longitudinal equivalence with full sysEqt ---
err3 = 0; err3lat = 0;
for k = 1:20
    x = zeros(10,1);
    x(1) = 10*randn; x(3) = 5*randn;          % vN, vD
    x(5) = 2*randn;                           % pitch rate
    the  = pi*randn;
    x(7:10) = [cos(the/2); 0; sin(the/2); 0]; % pure-pitch quaternion
    wp = 800*rand; d = 0.5*randn;
    w  = [randn; 0; randn];
    dx_full = sysEqt(x, [wp; -wp; d; d], w, drone);
    dx_long = sysEqt_long(x, [wp; d], w, drone);
    err3    = max(err3, max(abs(dx_full(lon) - dx_long(lon))));
    err3lat = max(err3lat, max(abs(dx_full(lat))));
end
fprintf('3. max |sysEqt - sysEqt_long| (long channels): %.3g\n', err3);
fprintf('   (full model lateral derivative on same pts: %.3g)\n', err3lat);

%% --- MATLAB-side verdict ---
tol = 1e-6;
if err1 < tol && err2 == 0 && err3 < 1e-9
    disp('check_Darko_long: function checks PASSED.');
else
    warning('check_Darko_long: function checks FAILED (see values above).');
end

%% --- 4. sim smoke test: Darko_long.slx from trim, open loop ---
Ts = 1/500;
delta_mass = 0; delta_inertia = 0; delta_phi_coef = 0;
time_var_param = inf;
xj0 = xt;
uj0_long = ut_long;

out = sim('Darko_long.slx', 'StopTime', '2');

lat_max = max(abs(out.states_long(:, lat)), [], 1);
drift   = max(abs(out.states_long(end, lon)' - xt(lon)));
fprintf('4. sim: max lateral states = %.3g, 2 s open-loop drift (long) = %.3g\n', ...
    max(lat_max), drift);

if max(lat_max) == 0
    disp('check_Darko_long: sim smoke test PASSED -- motion stayed in the longitudinal plane.');
else
    warning('check_Darko_long: lateral states left the longitudinal plane (max %.3g).', max(lat_max));
end
