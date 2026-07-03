%% build_AC_long.m
%
% Given the V0-dependent longitudinal state-space matrices A_V0, B_V0
% (evaluated at a known trim airspeed V0), recover the V0-INDEPENDENT
% "non-scaled" derivative matrices A, B.
%
% State vector:  x = [V*; gamma; alpha; q]
% Input vector:  u = [eta; deltaT]
%
% -------------------------------------------------------------------
% LaTeX -- scaled (V0-dependent) matrices, as originally given:
%
%   A(V_0) =
%   \begin{bmatrix}
%   X_V & -\dfrac{g}{V_0} & X_\alpha & X_q \\
%   -Z_V V_0 & 0 & -Z_\alpha & -Z_q \\
%   Z_V V_0 & 0 & Z_\alpha & Z_q+1 \\
%   M_V V_0 & 0 & M_\alpha & M_q
%   \end{bmatrix}
%   ,\qquad
%   B(V_0) =
%   \begin{bmatrix}
%   \dfrac{X_\eta}{V_0} & \dfrac{X_{\delta T}}{V_0} \\
%   -Z_\eta & -Z_{\delta T} \\
%   Z_\eta & Z_{\delta T} \\
%   M_\eta & M_{\delta T}
%   \end{bmatrix}
%
% LaTeX -- non-scaled (V0-independent) matrices, i.e. what this script
% computes as A, B below:
%
%   \bar{A} =
%   \begin{bmatrix}
%   X_V & -g & X_\alpha & X_q \\
%   -Z_V & 0 & -Z_\alpha & -Z_q \\
%   Z_V & 0 & Z_\alpha & Z_q+1 \\
%   M_V & 0 & M_\alpha & M_q
%   \end{bmatrix}
%   ,\qquad
%   \bar{B} =
%   \begin{bmatrix}
%   X_\eta & X_{\delta T} \\
%   -Z_\eta & -Z_{\delta T} \\
%   Z_\eta & Z_{\delta T} \\
%   M_\eta & M_{\delta T}
%   \end{bmatrix}
%
% Relation between the two: every entry of A(V0)/B(V0) equals the
% corresponding entry of Abar/Bbar times V0, EXCEPT column 1 of A and
% row 1 of B, which instead get divided by V0. (See extraction below.)
% -------------------------------------------------------------------

%% --- INPUTS: set these to your V0-dependent matrices and trim speed ---
V0 = 200;   % <-- trim airspeed [m/s] at which A_V0, B_V0 were computed

A_V0 = [-0.01  -0.049  -0.046  -0.001;
         0.11   0       1.14    0.043;
        -0.11   0      -1.14    0.957;
         0.1    0     -15.34   -3.00 ];

B_V0 = [ 0.0002   0.014;
        -0.059    0;
         0.059    0;
         6.56     0.13 ];

%% --- Extract the V0-independent derivatives from A_V0 ---
X_V     = A_V0(1,1);
g       = 9.81;                   % (Theoretically: -A_V0(1,2) * V0) -> A(1,2) = -g/V0
X_alpha = A_V0(1,3);
X_q     = A_V0(1,4);

Z_V_a   = -A_V0(2,1) / V0;        % A(2,1) = -Z_V*V0
Z_V_b   =  A_V0(3,1) / V0;        % A(3,1) =  Z_V*V0  (redundant check)
Z_V     = mean([Z_V_a, Z_V_b]);

Z_alpha_a = -A_V0(2,3);           % A(2,3) = -Z_alpha
Z_alpha_b =  A_V0(3,3);           % A(3,3) =  Z_alpha (redundant check)
Z_alpha   = mean([Z_alpha_a, Z_alpha_b]);

Z_q_a   = -A_V0(2,4);             % A(2,4) = -Z_q
Z_q_b   =  A_V0(3,4) - 1;         % A(3,4) =  Z_q + 1 (redundant check)
Z_q     = mean([Z_q_a, Z_q_b]);

M_V     = A_V0(4,1) / V0;
M_alpha = A_V0(4,3);
M_q     = A_V0(4,4);

%% --- Extract the V0-independent derivatives from B_V0 ---
X_eta = B_V0(1,1) * V0;           % B(1,1) = X_eta/V0
X_dT  = B_V0(1,2) * V0;           % B(1,2) = X_dT/V0

Z_eta_a = -B_V0(2,1);             % B(2,1) = -Z_eta
Z_eta_b =  B_V0(3,1);             % B(3,1) =  Z_eta (redundant check)
Z_eta   = mean([Z_eta_a, Z_eta_b]);

Z_dT_a  = -B_V0(2,2);             % B(2,2) = -Z_dT
Z_dT_b  =  B_V0(3,2);             % B(3,2) =  Z_dT  (redundant check)
Z_dT    = mean([Z_dT_a, Z_dT_b]);

M_eta = B_V0(4,1);
M_dT  = B_V0(4,2);

%% --- Consistency check on redundant (sign-paired) entries ---
tol = 1e-3;
mismatch = [abs(Z_V_a - Z_V_b), abs(Z_alpha_a - Z_alpha_b), ...
            abs(Z_q_a - Z_q_b), abs(Z_eta_a - Z_eta_b), abs(Z_dT_a - Z_dT_b)];
if any(mismatch > tol)
    warning('One or more redundant entries in A_V0/B_V0 disagree by more than %.g. Check transcription.', tol);
end

%% --- Assemble the non-scaled (V0-independent) matrices ---
A = [ X_V,   -g,      X_alpha,  X_q;
     -Z_V,    0,     -Z_alpha, -Z_q;
      Z_V,    0,      Z_alpha,  Z_q+1;
      M_V,    0,      M_alpha,  M_q ];

B = [ X_eta,   X_dT;
     -Z_eta,  -Z_dT;
      Z_eta,   Z_dT;
      M_eta,   M_dT ];

disp('Non-scaled (V0-independent) A matrix:');
disp(A);
disp('Non-scaled (V0-independent) B matrix:');
disp(B);

%% --- Sanity check: re-scale A, B back up to V0 and compare to input ---
A_check = A; B_check = B;
A_check(1,2) = A(1,2)/V0;   % -g/V0
A_check(2,1) = A(2,1)*V0;   % -Z_V*V0
A_check(3,1) = A(3,1)*V0;   %  Z_V*V0
A_check(4,1) = A(4,1)*V0;   %  M_V*V0
B_check(1,:) = B(1,:)/V0;   % X_eta/V0, X_dT/V0

fprintf('Max reconstruction error: A -> %.3g, B -> %.3g\n', ...
    max(abs(A_check(:)-A_V0(:))), max(abs(B_check(:)-B_V0(:))));