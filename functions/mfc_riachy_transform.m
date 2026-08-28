function [Y, state] = mfc_riachy_transform(y, Kd, Ts, state)
%MFC_RIACHY_TRANSFORM Riachy's auxiliary output Y = y + Kd*int y.
%
%   [Y, state] = MFC_RIACHY_TRANSFORM(y, Kd, Ts, state)
%
%   Riachy's trick (Riachy et al.) removes the derivative from an iPD/iPID
%   without estimating one. Start from the second-order ultra-local model
%   and add Kd*dot_y to both sides,
%
%       ddot_y + Kd*dot_y = F + Kd*dot_y + alpha*u
%
%   then define the auxiliary output
%
%       Y(t) = y(t) + Kd * int_c^t y(sigma) dsigma ,     0 <= c < t
%
%   so that ddot_Y = ddot_y + Kd*dot_y exactly. With
%
%       Fk = F + Kd*dot_y
%
%   the model reads  ddot_Y = Fk + alpha*u : the SAME second-order
%   ultra-local model, in Y instead of y. Any second-order F estimator --
%   algebraic or sliding window -- applied to Y therefore returns Fk, and
%   Fk already contains the derivative term the iPD would otherwise need.
%   The command law loses its Kd*dot_y feedback and keeps Kd only on the
%   reference:
%
%       u = -( Fk_est - ddot_sp - Kd*dot_sp + Kp*err + Ki*int err ) / alpha
%
%   with err = y - sp_filt (this repo's sign convention: the command law
%   SUBTRACTS fb). Closed loop: ddot_err + Kd*dot_err + Kp*err + Ki*int
%   err = 0, i.e. the same poles as the ordinary decoupled iPD -- obtained
%   with no derivative of the measurement anywhere in the loop.
%
%   DISCRETIZATION: trapezoidal integration of y,
%
%       int[k] = int[k-1] + (Ts/2)*(y[k] + y[k-1]),   Y[k] = y[k] + Kd*int[k]
%
%   Trapezoidal, not forward Euler: the estimators downstream weight the
%   window by powers of t, so a half-sample bias in Y is not free.
%
%   c IS FREE, and int y is unbounded. For a non-zero steady state y, Y
%   ramps forever. That is harmless for ddot_Y and for the sliding-window
%   estimator (finite memory), but the algebraic estimator's growing t^2
%   weights see an ever-larger signal -- see Knowledge/riachy-trick.md.
%
%   VECTOR-SAFE. y may be n-by-1 with a square n-by-n Kd (the MIMO
%   ultra-local model): Y = y + Kd*int y is a matrix product, and the
%   derivation is unchanged because adding Kd*dot_y to both sides is
%   linear. See MFC_FHAT_RIACHY2_MIMO_BLOCK.
%
%   Inputs
%     y     : plant measurement at step k (scalar or n-by-1)
%     Kd    : derivative gain folded into Y (the Kd of the closed-loop
%             polynomial s^2 + Kd*s + Kp; n-by-n for a vector y)
%     Ts    : sample time [s]
%     state : struct with fields
%               .int_km1  integral of y through step k-1
%               .y_km1    previous measurement
%
%   Outputs
%     Y     : auxiliary output y + Kd*int y
%     state : updated state struct
%
%   See also MFC_FHAT_RIACHY2_BLOCK, MFC_FHAT_RIACHY2_MIMO_BLOCK, MFC_FHAT_ALGEBRAIC_SECOND_ORDER,
%   MFC_FHAT_SLIDING_WINDOW.

int_y = state.int_km1 + 0.5*Ts*(y + state.y_km1);
Y     = y + Kd*int_y;

state.int_km1 = int_y;
state.y_km1   = y;
end
