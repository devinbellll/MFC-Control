%%% Build the SAME controller twice -- once as the all-in-one mfc_siso_core
%%% block, once assembled by hand from the individual stage blocks -- run
%%% both on identical plants, and overlay them.
%%%
%%% This is the template for tweaking the guts of the method. The assembled
%%% loop below is deliberately written out longhand: every stage is a
%%% separate object with its own parameters and its own state, so you can
%%% replace one of them (a different estimator, a smoother moved somewhere
%%% else, a hand-written feedback law) and immediately see the effect
%%% against the reference implementation.
%%%
%%% Plant: drone z-axis, as in val_mfc.m -- double integrator + constant
%%% gravity + linear drag + first-order actuator lag.
%%%
%%% THE TWO WIRING RULES the all-in-one block hides from you:
%%%
%%%   1. u_prev must be routed back EXPLICITLY. The estimator needs the
%%%      command that was applied over the last sample; mfc_siso_core keeps
%%%      it internally, an assembled loop needs a real unit delay (here,
%%%      the u_km1 variable).
%%%
%%%   2. the explicit feedback is YOURS. There is no feedback block: with a
%%%      DECOUPLED estimator you wire a stock Discrete PID Controller on the
%%%      error into the command block's fb input (mfc_siso.feedback stands
%%%      in for it below, so the comparison stays exact). With a COUPLED
%%%      estimator fb is Ground -- Kp and Kd are already inside F_hat.
%%%
%%% Composed loops have no anti-windup path at all; saturation is off here,
%%% and mfc_siso_core remains the reference if you need integral action
%%% against a real actuator limit.
%%%
%%% Run:  >> setup;  val_mfc_composed

clear; close all;

%%% plant
g = 9.81; d = 0.2; alpha_plant = 1; tau = 0.05;

%%% tuning (2nd order, decoupled, algebraic)
Ts = 0.01;  alpha = alpha_plant;
WFilter = 10;  FFilter = 10;
p = (1/tau)/5;  Kp = p^2;  Kd = 2*p;  Ki = 0;

%%% --- controller A: the all-in-one block ------------------------------
A = mfc_siso_core( ...
    'model_order',          'Second order (ddot_y = F + alpha*u)', ...
    'controller_structure', 'Decoupled (measurement-driven estimator, explicit iP/iPD)', ...
    'estimator_type',       'Algebraic (growing window)', ...
    'Ts', Ts, 'alpha', alpha, 'Kp', Kp, 'Kd', Kd, 'Ki', Ki, ...
    'ref_filter_window', WFilter, 'est_filter_window', FFilter);

%%% --- controller B: the same thing, assembled from stages -------------
% 1. input smoother -> sp_filt and the feedforward derivatives
B_smoother = mfc_smoother_block('Ts', Ts, 'window', WFilter, ...
                                'output_derivatives', true);
% 2. F-hat estimator. The DECOUPLED block: driven by the measurement,
%    nothing folded -- swap in mfc_fhat_alg2_coupled_block (fed the error,
%    given Kp/Kd) and drop the PID below to see the other structure.
B_est      = mfc_fhat_alg2_decoupled_block('Ts', Ts, 'alpha', alpha, ...
                                 'est_filter_window', FFilter, 'est_hold_time', 0.1);
% 3. model inversion (saturation available on its mask, off here)
B_command  = mfc_command_block('alpha', alpha);
% the explicit feedback is a stock Discrete PID in Simulink; its state is
% just err_km1 and the trapezoidal integral, kept here in the loop below.
B_err_km1 = 0;  B_int_err = 0;

%%% sim setup: each controller drives its OWN copy of the plant
t   = (0:Ts:6)';
ref = double(t >= 0.2);

yA = 0; dyA = 0; uActA = 0;
yB = 0; dyB = 0; uActB = 0;
u_km1 = 0;                    % <-- wiring rule 1: the explicit unit delay

UA = zeros(size(t)); FA = UA; YA = UA;
UB = zeros(size(t)); FB = UB; YB = UB;

for k = 1:numel(t)
    YA(k) = yA;  YB(k) = yB;

    % --- A: all-in-one ------------------------------------------------
    [UA(k), FA(k)] = step(A, ref(k), yA, t(k));

    % --- B: assembled from stages -------------------------------------
    [sp_filt, dot_sp, ddot_sp] = step(B_smoother, ref(k));
    err     = yB - sp_filt;                        % note the sign: y - sp
    FB(k)   = step(B_est, yB, u_km1, t(k));        % decoupled: driven by measurement
    [fb, B_int_err] = mfc_siso.feedback(err, B_err_km1, B_int_err, ...
                                        Ts, Kp, Kd, Ki, false);   % the PID
    B_err_km1 = err;
    ff      = ddot_sp;                             % 2nd order => ddot_sp
    UB(k)   = step(B_command, FB(k), ff, fb);
    u_km1   = UB(k);

    % --- plants --------------------------------------------------------
    uActA = uActA + (UA(k) - uActA)/tau * Ts;
    dyA   = dyA   + (-g - d*dyA + alpha_plant*uActA)*Ts;
    yA    = yA    + dyA*Ts;

    uActB = uActB + (UB(k) - uActB)/tau * Ts;
    dyB   = dyB   + (-g - d*dyB + alpha_plant*uActB)*Ts;
    yB    = yB    + dyB*Ts;
end

%%% --- report ----------------------------------------------------------
worst = max(abs(UA - UB));
fprintf('val_mfc_composed: max |u_allinone - u_composed| = %.3g\n', worst);
if worst < 1e-12
    fprintf('  the decomposition is exact.\n');
else
    fprintf('  MISMATCH -- check the wiring rules in the header.\n');
end

figure;
subplot(3,1,1);
plot(t, ref, 'k:', t, YA, 'b-', t, YB, 'r--', 'LineWidth', 1.2);
legend('command', 'all-in-one', 'composed from stages', 'Location', 'southeast');
ylabel('altitude z (m)'); grid on; ylim([-0.5 2]);
title('mfc\_siso\_core vs the same controller assembled from stage blocks');

subplot(3,1,2);
plot(t, UA, 'b-', t, UB, 'r--', 'LineWidth', 1.2);
ylabel('u (thrust cmd)'); grid on;

subplot(3,1,3);
plot(t, UA - UB, 'k-', 'LineWidth', 1.2);
ylabel('u difference'); xlabel('time (s)'); grid on;
