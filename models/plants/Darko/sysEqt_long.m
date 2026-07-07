function dxdt = sysEqt_long( x, u, w, drone )
%% SYSEQT_LONG tiltbody equations of motion, longitudinal plane only
%
%  this function computes the derivative of the system state constrained
%    to the longitudinal (north/down, pitch) plane: dx/dt = sysEqt_long(x,u,w)
%
%  it is a wrapper around SYSEQT that
%    1) projects the incoming state onto the longitudinal plane
%       (vE = 0, p = r = 0, pure-pitch quaternion),
%    2) expands the symmetric 2x1 input to the full 4x1 input
%       (counter-rotating props at equal speed, equal elevons),
%    3) zeroes the east wind component,
%    4) hard-zeroes the lateral derivative channels so lateral states
%       stay frozen under integration.
%
%  INPUTS:
%    state def. : x = (vl wb q)     in R(10x1)  (same as SYSEQT)
%    input def. : u = (wp d)        in R( 2x1)
%    disturbance: w = (wn we wd)    in R( 3x1)  (we is ignored)
%    drone specs: drone             struct      (same as SYSEQT)
%
%  OUTPUTS:
%    dxdt (state derivative)        in R(10x1)
%      with dxdt([2 4 6 8 10]) = 0  (vE, p, r, qx, qz channels)
%
%  wp: propeller angular speed magnitude (rad/s) [Scalar Real]
%      (expanded internally to w1 = wp, w2 = -wp)
%  d:  symmetric elevon deflection (rad) [Scalar Real]
%      (expanded internally to d1 = d2 = d)
%
%  NOTE1: with equal counter-rotating prop speeds the prop reaction and
%    gyroscopic moments cancel, so the prop sign convention does not
%    affect the longitudinal dynamics.
%  NOTE2: elevon sign convention is positive pitch-up deflections.
%
%  see also SYSEQT, LONGTRIM.

%% state projection onto longitudinal plane
xl = x;
xl(2) = 0;              % vE (east velocity)
xl(4) = 0;              % p  (roll rate)
xl(6) = 0;              % r  (yaw rate)
% pure-pitch quaternion: keep [q0; 0; qy; 0], renormalize
ql = [x(7); 0; x(9); 0];
xl(7:10) = ql/norm(ql);

%% input expansion to (w1 w2 d1 d2)
ul = [u(1); -u(1); u(2); u(2)];

%% wind projection (no east component)
wl = [w(1); 0; w(3)];

%% full phi-theory dynamics
dxdt = sysEqt(xl, ul, wl, drone);

%% enforce longitudinal constraint on the derivative
dxdt(2)  = 0;   % dvE
dxdt(4)  = 0;   % dp
dxdt(6)  = 0;   % dr
dxdt(8)  = 0;   % dqx
dxdt(10) = 0;   % dqz

end
