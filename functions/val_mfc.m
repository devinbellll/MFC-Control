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
%%% Plant (2nd-order truth):  ddot_y = a1*dot_y + a0*y + b*u
%%% Ideal 2nd-order model:    Kp / (s^2 + Kd*s + Kp)   (double pole -5)
%%% Ideal 1st-order model:    Kp / (s + Kp)
%%% NOTE: the 1st-order variants control a 2nd-order plant through a 1st-order
%%% ultra-local model, so they only approach the 1st-order ideal.

clear; close all;

%%% plant
a0 = 1;  a1 = 0.5;  b = 1;

%%% shared MFC tuning (Kd = 2p, Kp = p^2 -> double pole at -p, here p = 5)
Ts      = 0.01;
alpha   = b;              % matched gain
WFilter = 10;             % reference trajectory filter [samples]
FFilter = 10;             % estimator memory / window [samples]
Kp      = 25;
Kd      = 10;
Ki      = 0;

%%% variant grid
orders     = {'Second order (ddot_y = F + alpha*u)', ...
              'First order (dot_y = F + alpha*u)'};
structs    = {'Coupled (error-driven estimator, poles folded)', ...
              'Decoupled (measurement-driven estimator, explicit iP/iPD)'};
estimators = {'Algebraic (growing window)', ...
              'Sliding window (Simpson quadrature)'};

variants = { ...  % {order, structure, estimator, legend label}
    orders{1}, structs{1}, estimators{1}, '2nd coupled alg'; ...
    orders{1}, structs{2}, estimators{1}, '2nd decoupled alg'; ...
    orders{1}, structs{2}, estimators{2}, '2nd decoupled win'; ...
    orders{2}, structs{1}, estimators{1}, '1st coupled alg'; ...
    orders{2}, structs{2}, estimators{1}, '1st decoupled alg'; ...
    orders{2}, structs{2}, estimators{2}, '1st decoupled win'};
nc = size(variants, 1);

ctrl = cell(1, nc);
for c = 1:nc
    ctrl{c} = mfc_siso_core( ...
        'model_order',          variants{c, 1}, ...
        'controller_structure', variants{c, 2}, ...
        'estimator_type',       variants{c, 3}, ...
        'Ts', Ts, 'alpha', alpha, 'Kp', Kp, 'Kd', Kd, 'Ki', Ki, ...
        'ref_filter_window', WFilter, 'est_filter_window', FFilter);
end
names = variants(:, 4)';

%%% sim setup
t   = (0:Ts:6)';
ref = double(t >= 0.2);   % step command to 1
y   = zeros(1, nc);  dy = zeros(1, nc);   % plant state per controller
yi2 = 0;  dyi2 = 0;                       % ideal 2nd-order reference model
yi1 = 0;                                  % ideal 1st-order reference model

%%% logs
Y   = zeros(numel(t), nc);                % step response of each variant
Yi2 = zeros(numel(t), 1);  Yi1 = Yi2;     % ideal responses
U   = zeros(numel(t), nc);  F = U;  Yr = U;  E = U;

%%% run all variants in lockstep
for k = 1:numel(t)
    for c = 1:nc
        Y(k, c) = y(c);
        [U(k, c), F(k, c), Yr(k, c), E(k, c)] = step(ctrl{c}, ref(k), y(c), t(k));
        dy(c) = dy(c) + (a1*dy(c) + a0*y(c) + b*U(k, c))*Ts;
        y(c)  = y(c)  + dy(c)*Ts;
    end
    Yi2(k) = yi2;                          % ideal 2nd order: Kp/(s^2+Kd*s+Kp)
    dyi2   = dyi2 + (-Kd*dyi2 - Kp*yi2 + Kp*ref(k))*Ts;
    yi2    = yi2  + dyi2*Ts;
    Yi1(k) = yi1;                          % ideal 1st order: Kp/(s+Kp)
    yi1    = yi1  + Kp*(ref(k) - yi1)*Ts;
end

%%% step response: command, ideals, and each variant
figure;
plot(t, ref, 'k:', t, Yi2, 'k--', t, Yi1, 'k-.', t, Y, 'LineWidth', 1.2);
legend([{'command', 'ideal 2nd', 'ideal 1st'}, names], 'Location', 'southeast');
xlabel('time (s)'); ylabel('y'); grid on;
title('MFC SISO variants: step response');

%%% compare each controller output across variants
out = {U, F, Yr, E};  olabel = {'U', 'F\_hat', 'sp\_filt', 'err'};
figure;
for i = 1:4
    subplot(4, 1, i);
    plot(t, out{i}, 'LineWidth', 1.2);
    ylabel(olabel{i}); grid on;
    if i == 1, legend(names); title('controller output comparison'); end
end
xlabel('time (s)');
