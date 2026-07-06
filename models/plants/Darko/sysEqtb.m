function dxdt = sysEqtb( x, u, w, drone )
%% SYSEQTB tiltbody equations of motion in body axis
%
%  this function computes the derivative of the system state given any
%    state and inputs: dx/dt = sysEqt(x,u,w)
%
%  INPUTS:
%    state def. : x = (vb wb q)     in R(10x1)
%    input def. : u = (w1 w2 d1 d2) in R( 4x1)
%    disturbance: w = (wn we wd)    in R( 3x1)
%    drone specs: drone             struct
%    + required drone specs: 
%
%  OUTPUTS:
%    dxdt (state derivative)        in R(10x1) 
%
%  vb: vehicle velocity in body axis (m/s) [3x1 Real]
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

%% coordinate transformation
vb = x(1:3, 1);
wb = x(4:6, 1);
q  = x(7:10,1)/norm(x(7:10,1));
D = quat2dcm(q');
vl = D'*vb;
xl = [vl; wb; q];

%% computation of original derivative
dxdt0 = sysEqt( xl, u, w, drone );

%% computation of transformed derivative
dxdt = dxdt0;
dxdt(1:3,1) = D*dxdt(1:3,1) - cross(wb,vb);

end

