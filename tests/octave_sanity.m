% octave_sanity.m  --  MATLAB-free sanity check for the sliding-window F estimators.
%
%   Run:  octave --no-gui -q tests/octave_sanity.m
%         (also runs in MATLAB: >> run tests/octave_sanity.m)
%
%   Octave cannot instantiate matlab.System, so this file re-implements the
%   ESTIMATOR MATH that lives in functions/mfc_siso_non_algebraic.m (the
%   use_first_order=true/false branches of the estimator block of stepImpl)
%   and checks it against plants with an analytically known F. It is the
%   CI-friendly guard on the kernels / signs. For a test of the actual System
%   objects in MATLAB, use tests/test_estimators.m instead.
%
%   IMPORTANT: if you change the estimator formula in the source, mirror the
%   change in est1/est2 below. The asserts (constant-F recovery, sinusoid
%   tracking, and POSITIVE-sign correlation) will catch the historical bugs:
%     - 2nd order prefactor sign (-60/T^5 returns -F)
%     - 1st order u-term sign  (- u_kernel returns F + extra*u)

1;

% ---- estimator math, must match functions/*_F_estimator.m stepImpl ----
function w = simpweights(N)   % composite Simpson (N even): 1 4 2 .. 4 1
  w = ones(N+1,1); w(2:2:end-1) = 4; w(3:2:end-1) = 2;
end

function [F, st] = est1(ym, u, alpha, T, Ts, st)   % mfc_siso_non_algebraic, use_first_order=true (Eq.11)
  if isempty(st)
    N=round(T/Ts); N=N+mod(N,2); st.Tw=N*Ts;       % even intervals (Simpson)
    st.sig=(0:N)'*Ts; st.yk=st.Tw-2*st.sig; st.w=simpweights(N);
    st.y=zeros(N+1,1); st.u=st.y;
  end
  st.y=[st.y(2:end);ym]; st.u=[st.u(2:end);u];
  uk = alpha.*st.sig.*(st.Tw-st.sig);
  F  = (-6/st.Tw^3)*(Ts/3)*sum(st.w.*(st.yk.*st.y + uk.*st.u));   % +uk.*u, Simpson
end

function [F, st] = est2(ym, u, alpha, T, Ts, st)   % mfc_siso_non_algebraic, use_first_order=false (Eq.16)
  if isempty(st)
    N=round(T/Ts); N=N+mod(N,2); st.Tw=N*Ts;       % even intervals (Simpson)
    st.sig=(0:N)'*Ts; st.yk=st.Tw^2-6*st.Tw*st.sig+6*st.sig.^2; st.w=simpweights(N);
    st.y=zeros(N+1,1); st.u=st.y;
  end
  st.y=[st.y(2:end);ym]; st.u=[st.u(2:end);u];
  uk = (alpha/2).*st.sig.^2.*(st.Tw-st.sig).^2;
  F  = (60/st.Tw^5)*(Ts/3)*sum(st.w.*(st.yk.*st.y - uk.*st.u));   % +60/T^5, Simpson
end

% ---- tiny test harness ----
global N_PASS N_FAIL; N_PASS=0; N_FAIL=0;
function check(name, cond)
  global N_PASS N_FAIL;
  if cond, N_PASS++; printf('  PASS  %s\n', name);
  else      N_FAIL++; printf('  FAIL  %s\n', name); end
end

% Fine Ts so the trapezoidal quadrature is accurate and the precision
% tolerances below are meaningful (the 2nd-order estimate of a 2nd derivative
% needs a well-resolved window). Sign guards below catch the historical bugs
% regardless of Ts.
Ts=1e-4; alpha=2;

% ===== TEST 1: constant true F recovered (sign + scale) =====
% plant fed consistently: ydot/yddot = F0 + alpha*u0  ->  F = F0
printf('Constant-F recovery (true F = 3):\n');
T=0.1; F0=3; u0=1.5; t=0:Ts:1.5; u=u0*ones(size(t));
y1=(F0+alpha*u0)*t; y2=0.5*(F0+alpha*u0)*t.^2;
s1=[];s2=[];F1=0;F2=0;
for k=1:numel(t)
  [F1,s1]=est1(y1(k),u(k),alpha,T,Ts,s1);
  [F2,s2]=est2(y2(k),u(k),alpha,T,Ts,s2);
end
check(sprintf('1st order recovers F=3 (got %.3f)',F1), abs(F1-3) < 0.05);
check(sprintf('2nd order recovers F=3 (got %.3f)',F2), abs(F2-3) < 0.25);  % residual = trapz leakage on unbounded ramp
check('2nd order has CORRECT sign (not -F)', F2 > 0);

% ===== TEST 2: slow-varying F tracked with positive correlation =====
printf('Slow-sinusoid tracking:\n');
w0=0.5; t=0:Ts:3;
y=sin(w0*t); yd=w0*cos(w0*t); ydd=-w0^2*sin(w0*t); u=cos(w0*t);
F1t=yd-alpha*u; F2t=ydd-alpha*u;
s1=[];s2=[];F1=zeros(size(t));F2=zeros(size(t));
for k=1:numel(t)
  [F1(k),s1]=est1(y(k),u(k),alpha,T,Ts,s1);
  [F2(k),s2]=est2(y(k),u(k),alpha,T,Ts,s2);
end
m=t>0.5;
c1=corr(F1(m)',F1t(m)'); c2=corr(F2(m)',F2t(m)');
check(sprintf('1st order corr=%.4f > 0.99',c1), c1 > 0.99);
check(sprintf('2nd order corr=%.4f > 0.99',c2), c2 > 0.99);
check('2nd order corr POSITIVE (sign guard)', c2 > 0);
check(sprintf('1st order max|err|=%.4g < 0.05',max(abs(F1(m)-F1t(m)))), max(abs(F1(m)-F1t(m)))<0.05);
check(sprintf('2nd order max|err|=%.4g < 0.25',max(abs(F2(m)-F2t(m)))), max(abs(F2(m)-F2t(m)))<0.25);

% ===== TEST 3: COARSE Ts guard (Simpson vs trapezoidal regression) =====
% At Ts=0.01, T=0.1 the trapezoidal rule gave ~60x error on the 2nd-order
% estimate (F approx +17 instead of -1). Simpson must recover F here.
printf('Coarse-Ts robustness (Ts=0.01, T=0.1, constant true F=3):\n');
Tsc=0.01; Tc=0.1; tc=0:Tsc:2; uc=u0*ones(size(tc));
y1c=(F0+alpha*u0)*tc; y2c=0.5*(F0+alpha*u0)*tc.^2;
s1=[];s2=[];G1=0;G2=0;
for k=1:numel(tc)
  [G1,s1]=est1(y1c(k),uc(k),alpha,Tc,Tsc,s1);
  [G2,s2]=est2(y2c(k),uc(k),alpha,Tc,Tsc,s2);
end
check(sprintf('1st order @Ts=0.01 recovers F=3 (got %.3f)',G1), abs(G1-3) < 0.05);
check(sprintf('2nd order @Ts=0.01 recovers F=3 (got %.3f)',G2), abs(G2-3) < 0.1);

% ===== summary =====
printf('\n%d passed, %d failed\n', N_PASS, N_FAIL);
if N_FAIL > 0
  error('octave_sanity: %d test(s) failed', N_FAIL);
end
