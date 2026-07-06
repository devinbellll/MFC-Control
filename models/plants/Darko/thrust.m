function [T,N] = thrust( x, rot, del, w, drone )
%% THRUST thrust forces and moments modeling
%
%  this function computes the propeller thrust by means of simple
%    prop formulas commonly applied in quadrotor literature.
%
%  INPUTS:
%    state def. : x = (vl wb q)     in R(10x1)
%    motor rot. : rot               in R( 1x1)
%    elevon def.: del               in R( 1x1)
%    disturbance: w = (wn we wd)    in R( 3x1)
%    drone specs: drone             struct
%    + required drone specs: 
%    |---- PROP_KP                  in R( 1x1)
%    |---- PROP_KM                  in R( 1x1)
%
%  OUTPUTS:
%    Prop force : T (N) (body axis)  in R( 3x1) 
%    Prop moment: N (Nm) (body axis) in R( 3x1) 
%
%  vl: vehicle velocity in NED axis (m/s) [3x1 Real]
%  wb: vehicle angular velocity in body axis (rad/s) [3x1 Real]
%  q:  quaternion attitude (according to MATLAB convention) [4x1 Real]
%
%  NOTE1: notice that rot's sign depend on which section we are due to
%    counter-rotating propellers;
%  NOTE2: elevon sign convention is positive pictch-up deflections.
%  NOTE3: gyroscopic effects are implemented in the caller function!
%  
%  refer to [1] for further information.
% 
%  REFERENCES
%    [1] Lustosa L.R., Defay F., Moschetta J.-M., "The Phi-theory 
%    approach to flight control design of tail-sitter vehicles"
%    @ http://lustosa-leandro.github.io 

%% extract drone specs
kp = drone.PROP_KP;
km = drone.PROP_KM;

%% prop forces computation
% notice that negative thrust is not implemented
T = kp*rot^2*[1; 0; 0];

%% prop moments computation
N = sign(-rot)*km*rot^2*[1; 0; 0];

end

