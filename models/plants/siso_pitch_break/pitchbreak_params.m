%% pitchbreak_params.m  -- run before the Simulink model

clear P

%% Plant
P.c     = 0.8;                  % quadratic pitch damping     [1/rad]
P.k     = 60;                   % linear stiffness = wn^2     [1/s^2]
P.y_sad = 0.35;                 % pitch break angle           [rad]
P.gam   = P.k/P.y_sad^2;        % softening coefficient       [1/(rad^2 s^2)]
P.b     = -70;                  % control effectiveness       [1/s^2]

%% Actuator
P.tau  = 0.030;                 % first-order lag             [s]
P.rate = deg2rad(300);          % rate limit                  [rad/s]
P.umax = deg2rad(25);           % position limit              [rad]

%% Sampling and measurement noise
P.Ts       = Ts;
P.sigma_y  = deg2rad(0.3);
P.noise_pw = P.sigma_y^2 * P.Ts;    % Band-Limited White Noise "Noise power"
P.seed     = 23341;

%% Filters
P.wf = 60;                      % measurement pre-filter, 2nd-order Butterworth
P.filt_num = P.wf^2;
P.filt_den = [1, sqrt(2)*P.wf, P.wf^2];
P.N  = 30;                      % derivative filter corner

%% Controller
P.Kp = 60;
P.Ki = 200;
P.Kd = 15;
P.Kb = P.Ki/P.Kd;               % back-calculation anti-windup gain

%% Derived, for reference
P.wn_ol    = sqrt(P.k);                     % 7.75 rad/s
P.wn_cl    = sqrt(P.k + P.Kp);              % 11.0
P.w_act    = 1/P.tau;                       % 33.3
P.y_sad_cl = sqrt((P.k + P.Kp)/P.gam);      % 0.495 rad = 28.4 deg

%% Trajectory
P.traj.order = 7;               % 5 = quintic, 7 = zero jerk at endpoints
P.traj.y0    = 0;
P.traj.y1    = deg2rad(10);
P.traj.t0    = 1.0;             % hold at y0 before the move
P.traj.T     = 2.0;             % duration of the move
P.traj.tend  = 8.0;

t = (0:P.Ts:P.traj.tend).';
[r, rd, rdd] = rest_to_rest(t, P.traj);

% exact inverse of the plant along the reference
u_ff = (rdd + P.c*rd.*abs(rd) + P.k*r - P.gam*r.^3) / P.b;

ref = timeseries([r rd rdd u_ff], t);

report_margins(P, r, rd, u_ff);


%% ---------------------------------------------------------------------------
function [r, rd, rdd] = rest_to_rest(t, tr)
switch tr.order
    case 5, a = [0 0 0 10 -15 6];
    case 7, a = [0 0 0 0 35 -84 70 -20];
    otherwise, error('traj.order must be 5 or 7');
end
p   = a(end:-1:1);                      % descending, for polyval
pd  = polyder(p);
pdd = polyder(pd);

s      = min(max((t - tr.t0)/tr.T, 0), 1);
active = (t >= tr.t0) & (t <= tr.t0 + tr.T);

dy  = tr.y1 - tr.y0;
r   = tr.y0 + dy*polyval(p, s);
rd  = zeros(size(t));
rdd = zeros(size(t));
rd(active)  = dy*polyval(pd,  s(active)) / tr.T;
rdd(active) = dy*polyval(pdd, s(active)) / tr.T^2;
end

function report_margins(P, r, rd, u_ff)
fprintf('\n  open-loop wn        %6.2f rad/s\n', P.wn_ol);
fprintf('  closed-loop wn      %6.2f rad/s\n', P.wn_cl);
fprintf('  closed-loop saddle  %6.2f deg\n', rad2deg(P.y_sad_cl));
fprintf('  peak |r|            %6.2f deg   (%.0f%% of saddle)\n', ...
        rad2deg(max(abs(r))), 100*max(abs(r))/P.y_sad_cl);
fprintf('  peak |rdot|         %6.2f deg/s\n', rad2deg(max(abs(rd))));
fprintf('  peak |u_ff|         %6.2f deg\n\n', rad2deg(max(abs(u_ff))));

if max(abs(r)) > 0.8*P.y_sad_cl
    warning('reference approaches the closed-loop saddle');
end
end
