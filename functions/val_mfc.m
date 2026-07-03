%%% Compare the three MFC SISO controllers, 1-to-1, on identical plants.
%%% Each controller drives its OWN copy of the same plant, same params.
%%% Also overlay the IDEAL 2nd-order response set by Kp,Kd: Kp/(s^2+Kd*s+Kp).
%%% Plant (ultra-local truth):  ddot_y = a1*dot_y + a0*y + b*u

clear; close all;

%%% plant
a0 = 1;  a1 = 0.5;  b = 1;

%%% shared MFC params  (Kd=2p, Kp=p^2 -> double pole at -p, here p=5)
%%% FFilter is the SAME knob in all three (estimator memory, FFilter*Ts).
Ts      = 0.01;
alpha   = b;              % matched gain
WFilter = 10;
FFilter = 10;
Kp      = 25;
Kd      = 10;
Ki      = 0;
p       = sqrt(Kp);       % folded uses the pole directly (double pole at -p)

%%% controllers (each gets the same shared knobs)
ctrl  = { ...
    mfc_siso_non_algebraic('Ts',Ts,'FFilter',FFilter,'alpha',alpha, ...
                           'WFilter',WFilter,'Kp',Kp,'Kd',Kd,'Ki',Ki, ...
                           'use_first_order',false), ...
    mfc_siso_decoupled('Ts',Ts,'FFilter',FFilter,'alpha',alpha, ...
                       'WFilter',WFilter,'Kp',Kp,'Kd',Kd,'Ki',Ki), ...
    mfc_siso_core('Ts',Ts,'int_window',FFilter,'alpha',alpha, ...
                 'time_trajec',WFilter,'kp',p) };
names = {'non\_algebraic','decoupled','core'};
nc    = numel(ctrl);

%%% sim setup
t   = (0:Ts:6)';
ref = double(t >= 0.2);   % step command to 1
y   = zeros(1,nc);  dy = zeros(1,nc);   % plant state per controller
yi  = 0;            dyi = 0;            % ideal reference-model state

%%% logs
Y  = zeros(numel(t), nc);    % step response of each controller
Yi = zeros(numel(t), 1);     % ideal response
U  = zeros(numel(t), nc);  F = U;  Yr = U;  E = U;   % the 4 outputs

%%% run all in lockstep
for k = 1:numel(t)
    for c = 1:nc
        Y(k,c) = y(c);                          % measured output = response sample
        [U(k,c), F(k,c), Yr(k,c), E(k,c)] = step(ctrl{c}, ref(k), y(c), t(k));
        dy(c) = dy(c) + (a1*dy(c) + a0*y(c) + b*U(k,c))*Ts;
        y(c)  = y(c)  + dy(c)*Ts;
    end
    Yi(k) = yi;                                 % ideal 2nd-order ref model
    dyi   = dyi + (-Kd*dyi - Kp*yi + Kp*ref(k))*Ts;
    yi    = yi  + dyi*Ts;
end

%%% step response: command, ideal, and each controller
figure;
plot(t, ref, 'k:', t, Yi, 'k--', t, Y, 'LineWidth', 1.2);
legend([{'command','ideal'}, names]); xlabel('time (s)'); ylabel('y'); grid on;
title('step response');

%%% compare each of the 4 outputs across controllers
out = {U, F, Yr, E};  olabel = {'U','F','yref\_filter','err'};
figure;
for i = 1:4
    subplot(4,1,i);
    plot(t, out{i}, 'LineWidth', 1.2);
    ylabel(olabel{i}); grid on;
    if i == 1, legend(names); title('controller output comparison'); end
end
xlabel('time (s)');
