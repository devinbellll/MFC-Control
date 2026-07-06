function cost = latCost( x, v0, w0, drone )
%LATCOST compute the difference between predicted state derivative at
%  candidate points and the state derivative at a coordinated curve.
%  
%  INPUTS:
%    search variable : x = (theta phi w1 w2 d1 d2)  in R( 6x1)
%    desired velocity: v0 (m/s)                     in R( 1x1)
%    desired velocity: w0 (rad/s)                   in R( 1x1)
%    drone                                          struct 
%    + required drone specs: 
%    |---- none in caller
%
%  OUTPUTS:
%    cost                                           in R(10x1) 
%
%  theta: vehicle pitch (rad)                [1x1 Real]
%  phi: vehicle roll (rad)                   [3x1 Real]
%  w1: left propeller angular speed (rad/s)  [Scalar Real]
%  w2: right propeller angular speed (rad/s) [Scalar Real]
%  d1: left elevon (rad)                     [Scalar Real]
%  d2: right elevon (rad)                    [Scalar Real]
%
%  refer to [1] for further information.
% 
%  REFERENCES
%    [1] Lustosa L.R., Defay F., Moschetta J.-M., "The Phi-theory 
%    approach to flight control design of tail-sitter vehicles"
%    @ http://lustosa-leandro.github.io 

%% demultiplex stuff
the = x(1);
phi = x(2);
w1 = x(3);
w2 = x(4);
d1 = x(5);
d2 = x(6);

%% solve for attitude
q = angle2quat(0,the,phi)';
D = quat2dcm(q');

%% given velocities
vl = [v0;0;0];
wl = [0;0;w0];
wb = D*wl;

%% state construction
x = [vl; wb; q ];
u = [w1; w2; d1; d2];
w = zeros(3,1);

%% cost computation
dxdt = sysEqt( x, u, w, drone );
cost = dxdt(1:6) - [ cross(wl,vl); zeros(3,1) ];
cost = diag([1 1 1 1 1 1])*cost;

end

