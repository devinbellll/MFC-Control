function dxdt = sysEqt( x, u, w, drone )
%% SYSEQT tiltbody equations of motion
%
%  this function computes the derivative of the system state given any
%    state and inputs: dx/dt = sysEqt(x,u,w)
%
%  INPUTS:
%    state def. : x = (vl wb q)     in R(10x1)
%    input def. : u = (w1 w2 d1 d2) in R( 4x1)
%    disturbance: w = (wn we wd)    in R( 3x1)
%    drone specs: drone             struct
%    + required drone specs: 
%    |---- G (gravity, m/s^2)      in R( 1x1)
%    |---- MASS (kg)               in R( 1x1)
%    |---- INERTIA (kg m^2)        in R( 3x3)
%    |---- P_P1_CG (m)             in R( 1x1)
%    |      (position of prop1 wrt cg)
%    |---- P_P2_CG (m)             in R( 1x1)
%    |      (position of prop2 wrt cg)
%    |---- P_A1_CG (m)             in R( 1x1)
%    |      (position of aero wrench1 wrt cg)
%    |---- P_A2_CG (m)             in R( 1x1)
%    |      (position of aero wrench2 wrt cg)
%    |---- INERTIA_PROP_X (kg m^2) in R( 1x1)
%    |---- INERTIA_PROP_N (kg m^2) in R( 1x1)
%
%  OUTPUTS:
%    dxdt (state derivative)        in R(10x1) 
%
%  vl: vehicle velocity in NED axis (m/s) [3x1 Real]
%  wb: vehicle angular velocity in body axis (rad/s) [3x1 Real]
%  q:  quaternion attitude (according to MATLAB convention) [4x1 Real]
%  w1: left propeller angular speed (rad/s) [Scalar Real]
%  w2: right propeller angular speed (rad/s) [Scalar Real]
%  d1: left elevon (rad) [Scalar Real]
%  d2: right elevon (rad) [Scalar Real]
%
%  NOTE1: notice that w1>0 while w2<0 (normally) due to counter-rotating
%    propellers;
%  NOTE2: elevon sign convention is positive pictch-up deflections.
%  
%  refer to [1] for further information.
% 
%  REFERENCES
%    [1] Lustosa L.R., Defay F., Moschetta J.-M., "The Phi-theory 
%    approach to flight control design of tail-sitter vehicles"
%    @ http://lustosa-leandro.github.io 

%% physical constants
% gravity
G = drone.G;
% vehicle mass
MASS = drone.MASS;
% vehicle inertia matrix
INERTIA = drone.INERTIA; % [3x3 real symmetric] in body axis
% position of propellers wrt center of gravity
P_P1_CG = drone.P_P1_CG; % [3x1 real] in body axis
P_P2_CG = drone.P_P2_CG; % [3x1 real] in body axis
% position of aerodynamic wrenches wrt center of gravity
P_A1_CG = drone.P_A1_CG; % [3x1 real] in body axis
P_A2_CG = drone.P_A2_CG; % [3x1 real] in body axis
% propeller blade inertia
INERTIA_PROP_X = drone.INERTIA_PROP_X; % [Scalar Real > 0]
INERTIA_PROP_N = drone.INERTIA_PROP_N; % [Scalar Real > 0]

%% state demultiplexing
wb = x(4:6);
q  = x(7:10);

%% input demultiplexing
w1 = u(1);
w2 = u(2);
d1 = u(3);
d2 = u(4);

%% linear velocity derivative
D = q2dcm(q');
% thrust computation
[T1,N1] = thrust(x, w1, d1, w, drone);
[T2,N2] = thrust(x, w2, d2, w, drone);
% aeordynamics computation (propwash included)
[A1,M1] = aero(x, T1, d1, w, drone);
[A2,M2] = aero(x, T2, d2, w, drone);
% all forces except gravity
Fb = T1 + T2 + A1 + A2;
% kinematics
dvdt = [0;0;G] + 1/MASS*D'*Fb;

%% angular velocity derivative
% Aerodynamic moments
Ma = M1 + M2 + cross(P_A1_CG, A1) + cross(P_A2_CG, A2);
% gyroscopic and reaction moments
P = wb(1); Q = wb(2); R = wb(3);
Jpx = INERTIA_PROP_X; Jpn = INERTIA_PROP_N;
tau1 = N1 - (P+w1)*(Jpx-Jpn)*[0;R;-Q];
tau2 = N2 - (P+w2)*(Jpx-Jpn)*[0;R;-Q];
Mg = tau1 + tau2 + cross(P_P1_CG, T1) + cross(P_P2_CG, T2);
% sum of moments
Mb = Ma + Mg;
% kinematics
dwdt = -INERTIA\cross(wb, INERTIA*wb) + INERTIA\Mb;

%% quaternion derivative
dqdt = zeros(4,1);
dqdt(1)   = -1/2*dot(wb, q(2:4));
dqdt(2:4) = 1/2*( q(1)*wb - cross( wb, q(2:4) ) );

%% derivative multiplexing
dxdt = [ dvdt; dwdt; dqdt ];

end

