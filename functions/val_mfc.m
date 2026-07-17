%%% Validate the MFC SISO controller family: run every supported variant of
%%% mfc_siso_core, 1-to-1, on identical plants (each variant drives its OWN
%%% copy of the same plant) and overlay the IDEAL closed-loop response.
%%%
%%% Variants (model order x structure x estimator; coupled+sliding_window is
%%% undefined and rejected by the block):
%%%   1. 2nd order, coupled,   algebraic
%%%   2. 2nd order, decoupled, algebraic
%%%   3. 2nd order, decoupled, sliding window
%%%   4. 1st order, coupled,   algebraic
%%%   5. 1st order, decoupled, algebraic
%%%   6. 1st order, decoupled, sliding window
%%%
%%% Plant: drone z-axis (altitude), double integrator + constant gravity
%%% disturbance + linear aerodynamic drag + first-order actuator (ESC/motor)
%%% lag on the thrust command:
%%%
%%%   dot_z      = v
%%%   dot_v      = -g - d*v + alpha_plant*u_act   (F = -g - d*v is unmodeled)
%%%   tau*dot_u_act = u_cmd - u_act                (first-order input delay)
%%%
%%% The MFC ultra-local model only ever sees "ddot_z = F + alpha*u_cmd" (or
%%% the 1st-order variant on dot_z); it has NO explicit knowledge of gravity,
%%% drag, or the actuator lag -- all are unmodeled dynamics that F_hat must
%%% either absorb (gravity and drag, both algebraic in the plant states, are
%%% well within reach of the estimator) or that the gain margin must tolerate
%%% (the lag adds phase delay MFC cannot estimate away).
%%%
%%% Ideal 2nd-order model (no lag, perfect F_hat): Kp2/(s^2+Kd2*s+Kp2)
%%% Ideal 1st-order model (no lag, perfect F_hat): Kp1/(s+Kp1)
%%% NOTE: the 1st-order variants control a plant with 2 extra unmodeled poles
%%% (velocity integrator + actuator lag) through a 1st-order ultra-local
%%% model, so they only approach the 1st-order ideal, and worse so than the
%%% 2nd-order variants approach the 2nd-order ideal.
%%%
%%% NOTE on "1st order, coupled": this plant is a true double integrator, so
%%% stabilizing z from a 1st-order model fundamentally requires derivative
%%% (velocity) feedback. The coupled structure folds Kp into the estimator
%%% but has NO room to fold a derivative term at 1st order (see
%%% mfc_fhat_algebraic_first_order.m) -- so it cannot add that damping at
%%% any gain, and is expected to diverge on this plant. "1st order,
%%% decoupled" has no such restriction (Kd is explicit at 1st order too) and
%%% is tuned below with Kd1 = p to supply the missing damping directly.

clear; close all;

%%% plant: drone z-axis (SI-ish units)
g           = 9.81;   % gravity -- constant disturbance the estimator must reject
d           = 0.2;     % linear aerodynamic drag [1/s] on vertical velocity
alpha_plant = 1;       % thrust-accel effectiveness (m/s^2 per unit command)
tau         = 0.05;    % actuator (ESC/motor) lag time constant [s] -> bandwidth 1/tau = 20 rad/s

%%% shared MFC tuning
Ts      = 0.01;
alpha   = alpha_plant;   % matched gain (controller assumes the plant's true thrust effectiveness)
WFilter = 10;             % reference trajectory filter [samples]
FFilter = 10;             % estimator memory / window [samples]
Ki      = 0;              % F_hat already rejects the constant gravity disturbance; no integral needed

% Closed-loop bandwidth set from the actuator bandwidth: keep the requested
% closed-loop pole(s) at roughly 1/5 of the actuator bandwidth (1/tau) so the
% unmodeled actuator phase lag doesn't eat the stability margin.
actuator_bw = 1/tau;             % 20 rad/s
p           = actuator_bw / 5;   % 4 rad/s closed-loop bandwidth target

Kp2 = p^2;   Kd2 = 2*p;   % 2nd order: double pole at -p  (s^2 + Kd2*s + Kp2)
Kp1 = p;     Kd1 = p;     % 1st order: Kp1 = p matches the ideal 1-pole model; Kd1 = p supplies
                          % the derivative (velocity) damping this double-integrator plant needs
                          % (only takes effect when decoupled -- coupled ignores Kd at 1st order)

%%% variant grid
orders     = {'Second order (ddot_y = F + alpha*u)', ...
              'First order (dot_y = F + alpha*u)'};
structs    = {'Coupled (error-driven estimator, poles folded)', ...
              'Decoupled (measurement-driven estimator, explicit iP/iPD)'};
estimators = {'Algebraic (growing window)', ...
              'Sliding window (Simpson quadrature)'};

variants = { ...  % {order, structure, estimator, legend label, Kp, Kd}
    orders{1}, structs{1}, estimators{1}, '2nd coupled alg',   Kp2, Kd2; ...
    orders{1}, structs{2}, estimators{1}, '2nd decoupled alg', Kp2, Kd2; ...
    orders{1}, structs{2}, estimators{2}, '2nd decoupled win', Kp2, Kd2; ...
    orders{2}, structs{1}, estimators{1}, '1st coupled alg',   Kp1, Kd1; ...
    orders{2}, structs{2}, estimators{1}, '1st decoupled alg', Kp1, Kd1; ...
    orders{2}, structs{2}, estimators{2}, '1st decoupled win', Kp1, Kd1};
nc = size(variants, 1);

ctrl = cell(1, nc);
for c = 1:nc
    ctrl{c} = mfc_siso_core( ...
        'model_order',          variants{c, 1}, ...
        'controller_structure', variants{c, 2}, ...
        'estimator_type',       variants{c, 3}, ...
        'Ts', Ts, 'alpha', alpha, 'Kp', variants{c, 5}, 'Kd', variants{c, 6}, 'Ki', Ki, ...
        'ref_filter_window', WFilter, 'est_filter_window', FFilter);
end
names = variants(:, 4)';

%%% sim setup
t     = (0:Ts:6)';
ref   = double(t >= 0.2);              % step altitude command to 1 m
y     = zeros(1, nc);  dy = zeros(1, nc);   % altitude, vertical velocity per controller
u_act = zeros(1, nc);                        % lagged (actually delivered) thrust command
yi2   = 0;  dyi2 = 0;                        % ideal 2nd-order reference model
yi1   = 0;                                   % ideal 1st-order reference model

%%% logs
Y   = zeros(numel(t), nc);                % step response of each variant
Yi2 = zeros(numel(t), 1);  Yi1 = Yi2;     % ideal responses
U   = zeros(numel(t), nc);  F = U;  Yr = U;  E = U;

%%% run all variants in lockstep
for k = 1:numel(t)
    for c = 1:nc
        Y(k, c) = y(c);
        [U(k, c), F(k, c), Yr(k, c), E(k, c)] = step(ctrl{c}, ref(k), y(c), t(k));
        u_act(c) = u_act(c) + (U(k, c) - u_act(c))/tau * Ts;            % actuator lag
        dy(c)    = dy(c)    + (-g - d*dy(c) + alpha_plant*u_act(c))*Ts; % gravity + drag + thrust
        y(c)     = y(c)     + dy(c)*Ts;
    end
    Yi2(k) = yi2;                          % ideal 2nd order: Kp2/(s^2+Kd2*s+Kp2)
    dyi2   = dyi2 + (-Kd2*dyi2 - Kp2*yi2 + Kp2*ref(k))*Ts;
    yi2    = yi2  + dyi2*Ts;
    Yi1(k) = yi1;                          % ideal 1st order: Kp1/(s+Kp1)
    yi1    = yi1  + Kp1*(ref(k) - yi1)*Ts;
end

%%% step response: command, ideals, and each variant
figure;
plot(t, ref, 'k:', t, Yi2, 'k--', t, Yi1, 'k-.', t, Y, 'LineWidth', 1.2);
legend([{'command', 'ideal 2nd', 'ideal 1st'}, names], 'Location', 'southeast');
xlabel('time (s)'); ylabel('altitude z (m)'); grid on;
title(sprintf('MFC SISO variants: drone altitude step response (p=%.1f rad/s, tau=%.3fs)', p, tau));
ylim([-0.5, 2]);   % "1st order, coupled" is expected to diverge on this plant (see NOTE above);
                    % bound the axis so the other variants stay readable.

%%% compare each controller output across variants
out = {U, F, Yr, E};  olabel = {'U (thrust cmd)', 'F\_hat', 'sp\_filt', 'err'};
figure;
for i = 1:4
    subplot(4, 1, i);
    plot(t, out{i}, 'LineWidth', 1.2);
    ylabel(olabel{i}); grid on;
    if i == 1, legend(names); title('controller output comparison'); end
end
xlabel('time (s)');
