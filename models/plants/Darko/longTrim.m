function [xts, uts, N] = longTrim( the, drone )
%LONGTRIM deliver longitudinal trimpoints to a given pitch angle
%
%  this function computes the set of trimpoints for a given pitch angle
%  during longitudinal flight. Notice that the number of trim points is
%  proven [1] to be either 0, 1 or 2. Therefore, this fucntion returns a
%  varible matrix with N columns representing N possible solutions. N can
%  be 0, 1, or 2.
%
%  INPUTS:
%    pitch angle     : the   (rad)  in R(1x1)
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
%    |---- PROP_KP                  in R( 1x1)
%    |---- PROP_KM                  in R( 1x1)
%
%  OUTPUTS:
%    xts (trim states)              in R(10xN) 
%      (x = [vl wb q])
%    uts (trim inputs)              in R( 4xN) 
%      (u = [w1 w2 d1 d2])
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
%  NOTE3: trim points are compuited assuming no wind disturbances
%  
%  refer to [1] for further information.
% 
%  REFERENCES
%    [1] Lustosa L.R., Defay F., Moschetta J.-M., "The Phi-theory 
%    approach to flight control design of tail-sitter vehicles"
%    @ http://lustosa-leandro.github.io 

%% extract drone relevant specs
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
G = drone.G;
MASS = drone.MASS;
kp = drone.PROP_KP;

%% derivative data
% area computations
Sp = pi*Prop_R^2;
S  = 2*(Swet + Sdry); 
% longitudinal Phi computation
PHI_fv = [PHI_fv(1,1) PHI_fv(1,3); PHI_fv(3,1) PHI_fv(3,3)];
PHI_mv = [PHI_mv(2,1) PHI_mv(2,3)];

%% TRIM linear system construction
% TRIM matrix construction (see [1])
them_cross = [0 1; -1 0]*Thetam(2);
thef_cross = [0 1; -1 0]*Thetaf(2);
TRIM = zeros(3,4);
TRIM(1:2,1) = -1/2*RHO*S*PHI_fv*[cos(the); sin(the)];
TRIM(1:2,2) = 1/2*RHO*S*PHI_fv*thef_cross*[cos(the); sin(the)];
TRIM(1:2,3) = 1/2*S/Sp*PHI_fv*thef_cross*[1;0];
TRIM(1:2,4) = [2;0]-1/2*S/Sp*PHI_fv*[1;0];
TRIM(  3,1) = -1/2*RHO*S*chord*PHI_mv*[cos(the); sin(the)];
TRIM(  3,2) = 1/2*RHO*S*chord*PHI_mv*them_cross*[cos(the); sin(the)];
TRIM(  3,3) = 1/2*S/Sp*chord*PHI_mv*them_cross*[1;0];
TRIM(  3,4) = -1/2*S/Sp*chord*PHI_mv*[1;0];
% b vector construction (see [1])
bt = -MASS*G*[-sin(the); cos(the); 0];

%% linear system solving by explicit rref - init
sol = rref([TRIM bt]);
% kernel is given by x0 + epsilon*dx
x0 = [ sol(:,5); 0 ];
dx = [ -sol(:,4); 1];
% delT polynomial coefficients construction (ax^2+bx+c=0)
pa = dx(2)-dx(1)*dx(3);
pb = x0(2)-dx(3)*x0(1)-dx(1)*x0(3);
pc = -x0(1)*x0(3);
% polynomial discriminant
Delta = pb^2-4*pa*pc;

% check the degree of the delT polynomial
if (pa == 0 && pb == 0)
    No = 0;
    % notice that N=0 means no solution or infinite solutions
    % for us, both cases are bad!
elseif (pa == 0)
    No = 1;
else
    No = 2;
end

%% check if we are in delta-inneficient angles
% in this case, none of the below computations are valid
if ( sin(the) == 0 || (PHI_mv*[cos(the);sin(the)]) == 0 )
    xts = [];
    uts = [];
    N = 0;
    return;
end

%% linear system solving by explicit rref - computation

% a priori, we have no solution:
N = 0;
xts = [];
uts = [];

% in case the delT polynomial is second-order
if Delta >= 0 && No == 2
    % compute the two real roots of polynomial
    epss = [(-pb-sqrt(Delta))/(2*pa); (-pb+sqrt(Delta))/(2*pa)];
    % init output
    xts = zeros(10,2);
    uts = zeros(4,2);
    % deliver solutions
    for i=1:2
        % compute adequate solution to linear system
        x = x0 + epss(i)*dx;
        % check if there is a valid freesream velocity
        if x(1) >= 0
            xt = zeros(10,1);
            vt = sqrt(x(1));
            Tt = x(4);
            dt = x(3)/x(4);
            % set NED velocity
            xt(1:3,1) = [ vt; 0; 0 ];
            % set quaternion
            xt(7:10,1) = angle2quat(0, the, 0);
            % set motor rotations
            wt = sqrt(Tt/kp);
            ut(1:2,1) = wt*[-1;1];
            % set elevon deflections
            ut(3:4,1) = dt*[1;1];
            % update the list of trim points, if we have positive thrust
            if (Tt>0)
                N = N+1;
                % the use of N as index garantees that the solutions fill
                % the output vector from left to right no matter in which
                % order we find appropriate epsilons and solutions.
                xts(:,N) = xt; 
                uts(:,N) = ut;
            end
        end
    end
end

% in case the delT polynomial is first-order
if No == 1
    % init output
    xts = zeros(10,1);
    uts = zeros(4,1);
    % compute the one real root of polynomial
    epss = -pc/pb; % pb is different than 0 due to poly order
    % compute adequate solution to linear system
    x = x0 + epss*dx;
    % check if there is a valid freesream velocity
    if x(1) >= 0
        xt = zeros(10,1);
        vt = sqrt(x(1));
        Tt = x(4);
        dt = x(3)/x(4);
        % set NED velocity
        xt(1:3,1) = [ vt; 0; 0 ];
        % set quaternion
        xt(7:10,1) = angle2quat(0, the, 0);
        % set motor rotations
        wt = sqrt(Tt/kp);
        ut(1:2,1) = wt*[1;-1];
        % set elevon deflections
        ut(3:4,1) = dt*[1;1];
        % update the list of trim points, if we have positive thrust
        if (Tt>0)
            N = 1; % we found one solution!
            xts(:,1) = xt; % here it is..
            uts(:,1) = ut;
        end
    end
end

end

