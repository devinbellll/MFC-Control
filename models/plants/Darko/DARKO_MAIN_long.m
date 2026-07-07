%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   DARKO Flight Simulator - LONGITUDINAL   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%          Model : Phi-Theory           %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   Plant: Darko_long.slx (x-z plane only)   %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Longitudinal-only variant of DARKO_MAIN.m. The plant Darko_long.slx is a
% copy of Darko.slx whose dynamics block calls sysEqt_long instead of
% sysEqt: motion is constrained to the north/down (pitch) plane and the
% input is reduced to u = [w_prop; delta] (symmetric props / elevons).
%
% Differences wrt Darko.slx worth knowing:
%   - the always-on +60% parameter inflation hard-coded in Darko.slx's
%     dynamics block is replaced by the delta_* logic below: perturbation
%     switches on when t > time_var_param (s). Defaults keep it OFF.
%   - initial conditions come from the (fast, analytic) longTrim instead
%     of the lateral fsolve trim grid of compute_drone_initial_points.
%   - the U and Wind inports are fed by constants (uj0_long and zero
%     wind), so the model holds trim open loop; replace the "u trim long"
%     constant with your controller when testing.
clear all;
close all;
clc;

%% Simulation parameters
t_simu = 10;
Ts     = 1/500;    % sample time of the discrete-rate (measurement) blocks

%% Delta parameters analysis (0 => nominal plant)
delta_mass      = 0;
delta_inertia   = 0;
delta_phi_coef  = 0;
time_var_param  = inf; % perturbation switch-on time (s); inf = never

%% Longitudinal trim point
% theta_init = 90 deg => hover flight; lower angles => forward flight
theta_init = 5*pi/180;
trim_sol   = 1;        % longTrim can return up to 2 solutions; pick one

drone = generateDarko();
[xts, uts, N] = longTrim(theta_init, drone);
if N < 1
    error('longTrim found no trim point for theta = %g deg.', theta_init*180/pi);
end
trim_sol = min(trim_sol, N);
xj0      = xts(:, trim_sol);                       % [vl; wb; quat], 10x1
uj0_long = [abs(uts(1, trim_sol)); uts(3, trim_sol)]; % [w_prop; delta]

fprintf('Trim @ theta = %g deg: V = %.3f m/s, w_prop = %.1f rad/s, delta = %.2f deg\n', ...
    theta_init*180/pi, xj0(1), uj0_long(1), uj0_long(2)*180/pi);

%% Run simulation Darko_long
disp('Running the simulation ...')
out = sim('Darko_long.slx', 'StopTime', num2str(t_simu));

%% Plot longitudinal results
disp('Plotting results ...')
% states_long columns: [vN vE vD p q r q0 qx qy qz]
t     = out.tout;
vN    = out.states_long(:,1);
vD    = out.states_long(:,3);
q_rate = out.states_long(:,5);
theta = 2*atan2(out.states_long(:,9), out.states_long(:,7));  % pure-pitch quaternion

% lateral channels must stay identically zero
lat_max = max(abs(out.states_long(:, [2 4 6 8 10])), [], 1);
fprintf('max lateral states [vE p r qx qz] = [%.3g %.3g %.3g %.3g %.3g]\n', lat_max);

% x-z trajectory by integrating the NED velocity
pos_n =  cumtrapz(t, vN);
alt   = -cumtrapz(t, vD);

figure(1)
subplot(2,2,1)
plot(pos_n, alt); grid on;
xlabel('North [m]'); ylabel('Altitude [m]'); title('x-z trajectory')
subplot(2,2,2)
plot(t, theta*180/pi); grid on;
xlabel('t [s]'); ylabel('\theta [deg]'); title('Pitch angle')
subplot(2,2,3)
plot(t, vN, t, -vD); grid on;
xlabel('t [s]'); legend('v_N', '-v_D'); ylabel('[m/s]'); title('NED velocity')
subplot(2,2,4)
plot(t, q_rate*180/pi); grid on;
xlabel('t [s]'); ylabel('q [deg/s]'); title('Pitch rate')
