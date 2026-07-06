function [ Ac, Bc, Cc, Dc ] = jacob( x, u, drone )
%% JACOB jacobian of equations of motion in euler angles differentials
%
%  this functions linearizes the nonlinear equations of motion with respect
%    to a given operation point. the local linearized equations are
%    described in function of local euler angle misalignment instead of
%    local linear quaternions (which are less meaningful and diffcult to
%    work with.
%
%  INPUTS:
%    state def. : x = (vb wb q)     in R(10x1)
%    input def. : u = (w1 w2 d1 d2) in R( 4x1)
%    drone specs: drone             struct
%    + required drone specs: 
%
%  OUTPUTS:
%    Ac                             in R(9x9)
%    Bc                             in R(9x4)
%    Cc                             in R(9x9)
%    Dc                             in R(9x4)
%
%  vb: vehicle velocity in body axis (m/s) [3x1 Real]
%  wb: vehicle angular velocity in body axis (rad/s) [3x1 Real]
%  q:  quaternion attitude (according to MATLAB convention) [4x1 Real]
%  w1: left propeller angular speed (rad/s) [Scalar Real]
%  w2: right propeller angular speed (rad/s) [Scalar Real]
%  d1: left elevon (rad) [Scalar Real]
%  d2: right elevon (rad) [Scalar Real]
%
%  NOTE1: notice that w1<0 while w2>0 (normally) due to counter-rotating
%    propellers;
%  NOTE2: elevon sign convention is positive pictch-up deflections.
%  
%  refer to [1] for further information.
% 
%  REFERENCES
%    [1] Lustosa L.R., Defay F., Moschetta J.-M., "The Phi-theory 
%    approach to flight control design of tail-sitter vehicles"
%    @ http://lustosa-leandro.github.io 

%% jacobians init
Ac = zeros(9,9);
Bc = zeros(9,4);
Cc = eye(9);
Dc = zeros(9,4);

%% input multiplexing
vb = x(1:3, 1);
wb = x(4:6, 1);
q  = x(7:10,1);
w1 = u(1,1);
w2 = u(2,1);
d1 = u(3,1);
d2 = u(4,1);

%% drone physical values
m = drone.MASS;
rho = drone.RHO;
c  = drone.CHORD;
ws = drone.WINGSPAN;
Phi_fv = drone.PHI(1:3,1:3);
Phi_mv = drone.PHI(4:6,1:3);
Phi_mw = drone.PHI(4:6,4:6);
I = eye(3);
xif = drone.ELEVON_FEFFICIENCY;
xim = drone.ELEVON_MEFFICIENCY;
J = drone.INERTIA;
Jpx = drone.INERTIA_PROP_X;
Jpn = drone.INERTIA_PROP_N;
Sp = pi*(drone.PROP_RADIUS)^2;
kf = drone.PROP_KP;
km = drone.PROP_KM;

%% quantities of interest
D  = quat2dcm(q');
w_cross = [0 -wb(3) wb(2); wb(3) 0 -wb(1); -wb(2) wb(1) 0];
v_cross = [0 -vb(3) vb(2); vb(3) 0 -vb(1); -vb(2) vb(1) 0];
gb = D*[0;0;drone.G];
g_cross =   [0 -gb(3) gb(2); gb(3) 0 -gb(1); -gb(2) gb(1) 0];
xif_cross = [0 -xif(3) xif(2); xif(3) 0 -xif(1); -xif(2) xif(1) 0];
xim_cross = [0 -xim(3) xim(2); xim(3) 0 -xim(1); -xim(2) xim(1) 0];
a1 = drone.P_A1_CG;
a2 = drone.P_A2_CG;
p1 = drone.P_P1_CG;
p2 = drone.P_P2_CG;
a1_cross = [0 -a1(3) a1(2); a1(3) 0 -a1(1); -a1(2) a1(1) 0];
a2_cross = [0 -a2(3) a2(2); a2(3) 0 -a2(1); -a2(2) a2(1) 0];
p1_cross = [0 -p1(3) p1(2); p1(3) 0 -p1(1); -p1(2) p1(1) 0];
p2_cross = [0 -p2(3) p2(2); p2(3) 0 -p2(1); -p2(2) p2(1) 0];
Jw = J*wb;
Jw_cross = [0 -Jw(3) Jw(2); Jw(3) 0 -Jw(1); -Jw(2) Jw(1) 0];
S = c*ws;
v = norm(vb);
if v > 0
    vhat = vb/v;
else
    vhat = [1 0 0]';
end
B = diag([ws c ws]);

%% jacobian formula (from [1]) for A

% partial vdot partial v
Ac(1:3,1:3) = Ac(1:3,1:3) - w_cross;
Ac(1:3,1:3) = Ac(1:3,1:3) - 1/m*1/4*rho*S*Phi_fv*(2*I - d1*xif_cross - d2*xif_cross)*(v*I+vhat*vb'); 
Ac(1:3,1:3) = Ac(1:3,1:3) - 1/m*1/4*rho*S*Phi_mv*B*(2*I-d1*xif_cross-d2*xif_cross)*wb*vhat';

% partial vdot partial w
Ac(1:3,4:6) = Ac(1:3,4:6) + v_cross;
Ac(1:3,4:6) = Ac(1:3,4:6) - 1/4*1/m*rho*S*Phi_mv*(2*I-d1*xif_cross-d2*xif_cross)*v*B;

% partial vdot partial psi
Ac(1:3,7:9) = Ac(1:3,7:9) + g_cross;

% partial wdot partial v
Ac(4:6,1:3) = Ac(4:6,1:3) - 1/2*rho*S*B*Phi_mv*(v*I+vhat*vb');
Ac(4:6,1:3) = Ac(4:6,1:3) - 1/2*rho*S*B*Phi_mw*B*wb*vhat';
Ac(4:6,1:3) = Ac(4:6,1:3) + 1/4*rho*S*(d1*a1_cross+d2*a2_cross)*Phi_fv*xif_cross*(v*I+vhat*vb');
Ac(4:6,1:3) = Ac(4:6,1:3) + 1/4*rho*S*(d1*a1_cross+d2*a2_cross)*Phi_mv*B*xif_cross*wb*vhat';
Ac(4:6,1:3) = Ac(4:6,1:3) + 1/4*rho*S*B*Phi_mv*xim_cross*(d1+d2)*(v*I+vhat*vb');
Ac(4:6,1:3) = Ac(4:6,1:3) + 1/4*rho*S*B*Phi_mw*B*xim_cross*(d1+d2)*wb*vhat';
Ac(4:6,1:3) = J\Ac(4:6,1:3);

% partial wdot partial w
Ac(4:6,4:6) = Ac(4:6,4:6) + Jw_cross;
Ac(4:6,4:6) = Ac(4:6,4:6) - w_cross*J;
Ac(4:6,4:6) = Ac(4:6,4:6) - 1/2*rho*S*B*Phi_mw*v*B;
Ac(4:6,4:6) = Ac(4:6,4:6) + 1/4*rho*S*(d1*a1_cross+d2*a2_cross)*Phi_mv*B*xif_cross*v;
Ac(4:6,4:6) = Ac(4:6,4:6) + 1/4*rho*S*B*Phi_mw*B*xim_cross*(d1+d2)*v;
Ac(4:6,4:6) = Ac(4:6,4:6) + (Jpx-Jpn)*[0 0 0; -wb(3) 0 -(wb(1)+w1); wb(2) (wb(1)+w1) 0]; 
Ac(4:6,4:6) = Ac(4:6,4:6) + (Jpx-Jpn)*[0 0 0; -wb(3) 0 -(wb(1)+w2); wb(2) (wb(1)+w2) 0];
Ac(4:6,4:6) = J\Ac(4:6,4:6);

% partial wdot partial psi
% this is zero!

% partial psidot partial v
% this is zero!

% partial psidot partial w
Ac(7:9,4:6) = eye(3);

% partial psidot partial psi
Ac(7:9,7:9) = -w_cross;

%% jacobian formula (from [1]) for B

% partial vdot partial w1
Bc(1:3,1) = Bc(1:3,1) + 1/m*(I - S/4/Sp*Phi_fv*(I-d1*xif_cross))*[2*kf*w1;0;0];

% partial vdot partial w2
Bc(1:3,2) = Bc(1:3,2) + 1/m*(I - S/4/Sp*Phi_fv*(I-d2*xif_cross))*[2*kf*w2;0;0];

% partial vdot partial d1
Bc(1:3,3) = Bc(1:3,3) + 1/m*1/4*rho*S*Phi_fv*xif_cross*v*vb;
Bc(1:3,3) = Bc(1:3,3) + 1/m*S/4/Sp*Phi_fv*xif_cross*kf*w1^2*[1;0;0];
Bc(1:3,3) = Bc(1:3,3) + 1/m*1/4*rho*S*Phi_mv*B*xif_cross*v*wb;

% partial vdot partial d2
Bc(1:3,4) = Bc(1:3,4) + 1/m*1/4*rho*S*Phi_fv*xif_cross*v*vb;
Bc(1:3,4) = Bc(1:3,4) + 1/m*S/4/Sp*Phi_fv*xif_cross*kf*w2^2*[1;0;0];
Bc(1:3,4) = Bc(1:3,4) + 1/m*1/4*rho*S*Phi_mv*B*xif_cross*v*wb;

% partial wdot partial w1
Bc(4:6,1) = Bc(4:6,1) - S/4/Sp*B*Phi_mv*[2*kf*w1;0;0];
Bc(4:6,1) = Bc(4:6,1) + p1_cross*[2*kf*w1;0;0];
Bc(4:6,1) = Bc(4:6,1) - S/4/Sp*a1_cross*Phi_fv*[2*kf*w1;0;0];
Bc(4:6,1) = Bc(4:6,1) + S/4/Sp*a1_cross*Phi_fv*xif_cross*d1*[2*kf*w1;0;0];
Bc(4:6,1) = Bc(4:6,1) + S/4/Sp*B*Phi_mv*xim_cross*d1*[2*kf*w1;0;0];
Bc(4:6,1) = Bc(4:6,1) + sign(-w1)*2*km*w1*[1;0;0];
Bc(4:6,1) = Bc(4:6,1) - (Jpx-Jpn)*[0;wb(3);-wb(2)];
Bc(4:6,1) = J\Bc(4:6,1);

% partial wdot partial w2
Bc(4:6,2) = Bc(4:6,2) - S/4/Sp*B*Phi_mv*[2*kf*w2;0;0];
Bc(4:6,2) = Bc(4:6,2) + p2_cross*[2*kf*w2;0;0];
Bc(4:6,2) = Bc(4:6,2) - S/4/Sp*a2_cross*Phi_fv*[2*kf*w2;0;0];
Bc(4:6,2) = Bc(4:6,2) + S/4/Sp*a2_cross*Phi_fv*xif_cross*d2*[2*kf*w2;0;0];
Bc(4:6,2) = Bc(4:6,2) + S/4/Sp*B*Phi_mv*xim_cross*d2*[2*kf*w2;0;0];
Bc(4:6,2) = Bc(4:6,2) + sign(-w2)*2*km*w2*[1;0;0];
Bc(4:6,2) = Bc(4:6,2) - (Jpx-Jpn)*[0;wb(3);-wb(2)];
Bc(4:6,2) = J\Bc(4:6,2);

% partial wdot partial d1
Bc(4:6,3) = Bc(4:6,3) + S/4/Sp*a1_cross*Phi_fv*xif_cross*[kf*w1^2;0;0];
Bc(4:6,3) = Bc(4:6,3) + 1/4*rho*S*a1_cross*Phi_fv*xif_cross*v*vb;
Bc(4:6,3) = Bc(4:6,3) + 1/4*rho*S*a1_cross*Phi_mv*B*xif_cross*v*wb;
Bc(4:6,3) = Bc(4:6,3) + 1/4*rho*S*B*Phi_mv*xim_cross*v*vb;
Bc(4:6,3) = Bc(4:6,3) + 1/4*rho*S*B*Phi_mw*B*xim_cross*v*wb;
Bc(4:6,3) = Bc(4:6,3) + S/4/Sp*B*Phi_mv*xim_cross*[kf*w1^2;0;0];
Bc(4:6,3) = J\Bc(4:6,3);

% partial wdot partial d2
Bc(4:6,4) = Bc(4:6,4) + S/4/Sp*a2_cross*Phi_fv*xif_cross*[kf*w2^2;0;0];
Bc(4:6,4) = Bc(4:6,4) + 1/4*rho*S*a2_cross*Phi_fv*xif_cross*v*vb;
Bc(4:6,4) = Bc(4:6,4) + 1/4*rho*S*a2_cross*Phi_mv*B*xif_cross*v*wb;
Bc(4:6,4) = Bc(4:6,4) + 1/4*rho*S*B*Phi_mv*xim_cross*v*vb;
Bc(4:6,4) = Bc(4:6,4) + 1/4*rho*S*B*Phi_mw*B*xim_cross*v*wb;
Bc(4:6,4) = Bc(4:6,4) + S/4/Sp*B*Phi_mv*xim_cross*[kf*w2^2;0;0];
Bc(4:6,4) = J\Bc(4:6,4);

end

