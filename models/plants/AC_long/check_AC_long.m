clear

%% check_AC_long.m
%
% Verify that the nonlinear AC_long.mdl linearizes at trim back to the
% V0-dependent matrices A_V0, B_V0 that build_AC_long.m starts from.
%
% The model's trim/equilibrium is (V = V0, gamma = 0, alpha = 0, q = 0)
% with zero inputs (eta and deltaT are deviations from trim), so the
% integrator initial conditions ARE the operating point.
%
% Note on state coordinates: A_V0/B_V0 are written for the normalized
% state x = [dV/V0; gamma; alpha; q], while linmod works in physical
% coordinates [dV; gamma; alpha; q]. The similarity transform
% T = diag([1/V0 1 1 1]) maps between them.

build_AC_long;   % defines V0, g, X_*, Z_*, M_* and A_V0, B_V0

%%
x0 = [200; 0; 0; 0];
mdl = 'AC_long_man';
load_system(mdl);

% Operating point in Simulink's internal state ordering
[~, x0, xstr] = feval(mdl, [], [], [], 0);
if ischar(xstr), xstr = cellstr(xstr); end

% Permutation from Simulink's state order to [V; gamma; alpha; q]
want = {'Int V', 'Int gamma', 'Int alpha', 'Int q'};
p = zeros(1, 4);
for k = 1:4
    p(k) = find(endsWith(strtrim(xstr), want{k}));
end

u0 = [0; 0];
[Alinmod, Blinmod, Clinmod, Dlinmod] = linmod(mdl, x0, u0);
Alinmod = Alinmod(p, p);   Blinmod = Blinmod(p, :);   Clinmod = Clinmod(:, p);

% A_V0(1,2) was transcribed rounded (-0.049 vs exact -g/V0 = -0.04905);
% compare against the theoretically consistent value, as build_AC_long.m
% does by hardcoding g = 9.81.
A_ref = A;

errA = max(abs(Alinmod(:) - A_ref(:)));
errB = max(abs(Blinmod(:) - B(:)));
errC = max(abs(Clinmod(:) - reshape(eye(4), [], 1)));
errD = max(abs(Dlinmod(:)));

fprintf('max |Alinmod - A_ref| = %.3g\n', errA);
fprintf('max |Blinmod - B_V0|  = %.3g\n', errB);
fprintf('max |C - I|, |D|    = %.3g, %.3g\n', errC, errD);

tol = 1e-6;
if max([errA, errB, errC, errD]) < tol
    disp('check_AC_long: PASSED -- nonlinear model linearizes to A_V0, B_V0 at trim.');
else
    warning('check_AC_long: linearization mismatch exceeds %.g.', tol);
end

Alinmod

A_ref

%% --- Time-domain check: 10 s free response to an alpha perturbation ---
% Simulate the Simulink nonlinear model from trim with an initial alpha
% offset (zero inputs) and compare against the original linear model
% xdot = A_V0*x in the normalized state x = [dV/V0; gamma; alpha; q].
alpha_pert = 10*pi/180;   % 2 deg initial alpha perturbation
t = (0:0.01:3)';        % common, evenly spaced time grid

set_param([mdl '/Int alpha'], 'InitialCondition', num2str(alpha_pert, 17));
restoreIC = onCleanup(@() set_param([mdl '/Int alpha'], 'InitialCondition', 'x0(3)'));

out  = sim(mdl, 'StopTime', '3', 'SaveTime', 'on', 'SaveOutput', 'on', ...
           'SaveFormat', 'Array', 'ReturnWorkspaceOutputs', 'on', ...
           'OutputOption', 'SpecifiedOutputTimes', 'OutputTimes', 't');
y_nl = out.yout;                     % columns: [V gamma alpha q], on grid t

clear restoreIC                      % put the integrator IC back to 0
set_param(mdl, 'Dirty', 'off');

% Linear free response from the same perturbation via ss/initial, mapped
% back to physical outputs (V = V0 + V0*x1).
sys_lin = ss(A, B, eye(4), 0);
x0_lin  = [0; 0; alpha_pert; 0];
x_lin   = initial(sys_lin, x0_lin, t);
y_lin   = [V0 + x_lin(:, 1), x_lin(:, 2:4)];

names = {'V [m/s]', 'gamma [rad]', 'alpha [rad]', 'q [rad/s]'};
peak  = max(abs(y_lin - [V0, 0, 0, 0]), [], 1);   % perturbation amplitude
err   = max(abs(y_nl - y_lin), [], 1);
for k = 1:4
    fprintf('%-12s max|nl - lin| = %.3g  (%.2f%% of peak %.3g)\n', ...
        names{k}, err(k), 100*err(k)/max(peak(k), eps), peak(k));
end

figure('Name', 'AC_long: nonlinear (Simulink) vs linear A(V0), alpha perturbation');
for k = 1:4
    subplot(4, 1, k);
    plot(t, y_nl(:, k), 'b-', t, y_lin(:, k), 'r--');
    ylabel(names{k}); grid on;
    if k == 1
        title(sprintf('Free response, alpha(0) = %g deg', alpha_pert*180/pi));
        legend('nonlinear (Simulink)', 'linear A(V0)', 'Location', 'best');
    end
end
xlabel('t [s]');

% For a 2 deg perturbation the nonlinear/linear divergence over 10 s is
% well below 1% of each channel's peak (largest on V, from the quadratic
% aero terms); anything bigger points at a wiring/parameter problem.
tol_td = 0.01;
if all(err <= tol_td * max(peak, eps))
    disp('check_AC_long: time-domain PASSED -- nonlinear matches linear within 1% of peak.');
else
    warning('check_AC_long: time-domain mismatch exceeds %g%% of peak.', 100*tol_td);
end

