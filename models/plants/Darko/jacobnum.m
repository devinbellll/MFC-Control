function [ Ac, Bc, Cc, Dc ] = jacobnum( x, u, drone )
%% JACOBNUM jacobian of equations of motion numerically
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
Ac = zeros(10,10);
Bc = zeros(10,4);
Cc = eye(10);
Dc = zeros(10,4);

%% infinitesimal definition
dxn = 1e-5;
dun = 1e-5;

%% no wind here
w = zeros(3,1);

%% A computation
for i = 1:10
	dx = zeros(10,1);
	dx(i) = dxn;
	Ac(:,i) = 1/dxn*( sysEqtb( x+dx, u, w, drone ) - sysEqtb( x, u, w, drone ) );
end

%% B computation
for i = 1:4
	du = zeros(4,1);
	du(i) = dun;
	Bc(:,i) = 1/dun*( sysEqtb( x, u+du, w, drone ) - sysEqtb( x, u, w, drone ) );
end

end

