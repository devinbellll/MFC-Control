function [Fb, Mb] = aero( x, T, del, w, drone )
%% AERO aerodynamic forces and moments computation (in a wing section!)
%
%  this function computes the aerodynamic forces and moments in body axis
%    in view of Phi-theory in a wing section. it includes propwash effects 
%    due to thrust.
%
%  INPUTS:
%    state def. : x = (vl wb q)     in R(10x1)
%    thrust     : T (N)             in R( 3x1)
%    elevon def.: del (rad)         in R( 1x1)
%    motor speed: w (rad/s)         in R( 1x1)
%    drone specs: drone             struct
%    + required drone specs: 
%    |---- PHI                      in R( 6x6)
%    |---- RHO (kg/m^3)             in R( 1x1)
%    |---- WET_SURFACE (m^2)        in R( 1x1)
%    |---- DRY_SURFACE (m^2)        in R( 1x1)
%    |---- PHI_n                    in R( 1x1)
%    |---- CHORD (m)                in R( 1x1)
%    |---- WINGSPAN (m)             in R( 1x1)
%    |       (of the wing section, half of full drone)
%    |---- PROP_RADIUS (m)          in R( 1x1)
%    |---- ELEVON_MEFFICIENCY       in R( 3x1)
%    |---- ELEVON_FEFFICIENCY       in R( 3x1)
%
%  OUTPUTS:
%    Aero force : Fb (body axis)    in R( 3x1) 
%    Aero moment: Mb (body axis)    in R( 3x1) 
%
%  vl: vehicle velocity in NED axis (m/s) [3x1 Real]
%  wb: vehicle angular velocity in body axis (rad/s) [3x1 Real]
%  q:  quaternion attitude (according to MATLAB convention) [4x1 Real]
%
%  NOTE1: notice that w's sign depend on which section we are due to
%    counter-rotating propellers;
%  NOTE2: elevon sign convention is positive pictch-up deflections.
%  
%  refer to [1] for further information.
% 
%  REFERENCES
%    [1] Lustosa L.R., Defay F., Moschetta J.-M., "The Phi-theory 
%    approach to flight control design of tail-sitter vehicles"
%    @ http://lustosa-leandro.github.io 

%% data extraction from drone struct
PHI_fv = drone.PHI(1:3,1:3);
PHI_mv = drone.PHI(4:6,1:3);
PHI_mw = drone.PHI(4:6,4:6);
phi_n  = drone.PHI_n;
RHO    = drone.RHO;
Swet   = drone.WET_SURFACE;
Sdry   = drone.DRY_SURFACE;
chord  = drone.CHORD;
ws     = drone.WINGSPAN;
Prop_R = drone.PROP_RADIUS;
Thetam = drone.ELEVON_MEFFICIENCY;
Thetaf = drone.ELEVON_FEFFICIENCY;

%% derivative data
Sp = pi*Prop_R^2;

%% state demultiplexing
vl = x(1:3);
wb = x(4:6);
q  = x(7:10);

%% DCM computation
D = q2dcm(q');

%% freestream velocity computation in body frame
vinf = D*(vl-w);

%% computation of total wing section area
S = Swet + Sdry;

%% computation of chord matrix
B = diag([ws; chord; ws]);

%% eta computation
eta = sqrt( norm(vinf)^2 + phi_n*norm(B*wb)^2 );

%% force computation
% airfoil contribution
Fb = -1/2*RHO*S*eta*PHI_fv*vinf - 1/2*RHO*S*eta*PHI_mv*B*wb - 1/2*Swet/Sp*PHI_fv*T;
% elevon contribution
Fb = Fb + 1/2*RHO*S*eta*PHI_fv*cross(del*Thetaf,vinf) + 1/2*RHO*S*eta*PHI_mv*B*cross(del*Thetaf,wb) + 1/2*Swet/Sp*PHI_fv*cross(del*Thetaf,T);

%% moment computation
% airfoil contribution
Mb = -1/2*RHO*S*eta*B*PHI_mv*vinf - 1/2*RHO*S*eta*B*PHI_mw*B*wb - 1/2*Swet/Sp*B*PHI_mv*T;
% elevon contribution
Mb = Mb + 1/2*RHO*S*eta*B*PHI_mv*cross(del*Thetam,vinf) + 1/2*RHO*S*eta*B*PHI_mw*B*cross(del*Thetam,wb) + 1/2*Swet/Sp*B*PHI_mv*cross(del*Thetam,T);

end
