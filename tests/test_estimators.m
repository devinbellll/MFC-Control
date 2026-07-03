function test_estimators()
% test_estimators  Closed-loop regression tests for the MFC SISO System objects.
%
%   Run in MATLAB:  >> cd <repo root>;  test_estimators
%   (Requires the actual matlab.System objects, so this does NOT run in Octave;
%    use tests/octave_sanity.m for a MATLAB-free check of the estimator math.)
%
%   Exercises the REAL objects in functions/ via step():
%     - mfc_siso_non_algebraic(use_first_order=false)  (sliding-window iPD,  Eq.16)
%     - mfc_siso_non_algebraic(use_first_order=true)   (sliding-window iP,   Eq.11)
%     - mfc_siso_core(use_first_order=false)           (algebraic, validated point-of-record)
%     - mfc_siso_core(use_first_order=true)            (algebraic iP, folded control law)
%     - mfc_siso_decoupled                              (algebraic iPD)
%
%   Checks: setpoint tracking, steady-state F accuracy/sign, and a sign-flip
%   regression on an OPEN-LOOP-UNSTABLE plant that only stays bounded when F is
%   estimated with the correct sign (the -60/T^5 bug diverges here).

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..', 'functions'));

    global N_PASS N_FAIL; N_PASS = 0; N_FAIL = 0;
    Ts = 1e-3;

    % ---------------------------------------------------------------------
    % 1) mfc_siso_non_algebraic (use_first_order=false): track a step on a
    %    stable 2nd-order plant, and check the F output equals the true
    %    lumped dynamics at steady state.
    %    Plant: yddot = -a1*yd - a0*y + b*u,  alpha = b (matched).
    % ---------------------------------------------------------------------
    a0 = 1; a1 = 1.4; b = 2; alpha = 2; ref = 1;
    c = mfc_siso_non_algebraic('Ts',Ts,'FFilter',100,'alpha',alpha, ...
                               'WFilter',10,'Kp',25,'Kd',10,'Ki',20, ...
                               'use_first_order',false);
    [t,y,~,F] = sim2(c, a0,a1,b, Ts, 6, @(tt) ref*(tt>=0.2));
    fin = t > 5;
    check('2nd: step tracking (|y-ref|<0.02)',  abs(mean(y(fin))-ref) < 0.02);
    check('2nd: bounded response',              max(abs(y)) < 5);
    % steady state: ydd=yd=0 -> u_ss=a0*ref/b ; F_true = -alpha*u_ss
    Ftrue = -alpha * (a0*ref/b);
    check(sprintf('2nd: steady F=%.3f vs true %.3f',mean(F(fin)),Ftrue), ...
          abs(mean(F(fin)) - Ftrue) < 0.1);
    check('2nd: F has correct sign', mean(F(fin)) < 0);

    % ---------------------------------------------------------------------
    % 2) SIGN-FLIP REGRESSION: open-loop-unstable plant, weak PD so that correct
    %    F-cancellation is REQUIRED for stability. yddot = +4*y + 2*u.
    %    Correct (+60/T^5): bounded.  Bug (-60/T^5): diverges.
    % ---------------------------------------------------------------------
    cu = mfc_siso_non_algebraic('Ts',Ts,'FFilter',100,'alpha',2, ...
                                'WFilter',10,'Kp',1,'Kd',2,'Ki',0, ...
                                'use_first_order',false);
    [tu,yu] = sim2(cu, -4,0,2, Ts, 6, @(tt) 1.0*(tt>=0.2));
    check(sprintf('2nd: unstable plant stays bounded (max|y|=%.3g)',max(abs(yu))), ...
          max(abs(yu)) < 10);

    % ---------------------------------------------------------------------
    % 3) mfc_siso_non_algebraic (use_first_order=true): track a step on a
    %    stable 1st-order plant.
    %    Plant: ydot = -a0*y + b*u,  alpha = b.  Steady F = -alpha*u_ss.
    % ---------------------------------------------------------------------
    a0 = 1; b = 2; alpha = 2; ref = 1;
    c1 = mfc_siso_non_algebraic('Ts',Ts,'FFilter',100,'alpha',alpha, ...
                                'WFilter',10,'Kp',5,'Ki',5, ...
                                'use_first_order',true);
    [t1,y1,~,F1] = sim1(c1, a0,b, Ts, 6, @(tt) ref*(tt>=0.2));
    fin1 = t1 > 5;
    check('1st (non-algebraic): step tracking (|y-ref|<0.02)', abs(mean(y1(fin1))-ref) < 0.02);
    check('1st (non-algebraic): bounded response',             max(abs(y1)) < 5);
    F1true = -alpha * (a0*ref/b);
    check(sprintf('1st (non-algebraic): steady F=%.3f vs true %.3f',mean(F1(fin1)),F1true), ...
          abs(mean(F1(fin1)) - F1true) < 0.1);

    % ---------------------------------------------------------------------
    % 4) Regression: mfc_siso_core (folded, 2nd order) & mfc_siso_decoupled
    %    still track (validated estimators).
    % ---------------------------------------------------------------------
    cf = mfc_siso_core('Ts',Ts,'alpha',2,'int_window',5,'time_trajec',10, ...
                       'kp',5);
    [tf,yf] = sim2(cf, 1,1.4,2, Ts, 6, @(tt) 1.0*(tt>=0.2));
    check('core (2nd order): step tracking',   abs(mean(yf(tf>5))-1) < 0.05);
    check('core (2nd order): bounded',         max(abs(yf)) < 5);

    cd = mfc_siso_decoupled('Ts',Ts,'alpha',2,'FFilter',10,'WFilter',10, ...
                            'Kp',25,'Kd',10,'Ki',20);
    [td,yd2] = sim2(cd, 1,1.4,2, Ts, 6, @(tt) 1.0*(tt>=0.2));
    check('decoupled: step tracking', abs(mean(yd2(td>5))-1) < 0.05);
    check('decoupled: bounded',       max(abs(yd2)) < 5);

    % ---------------------------------------------------------------------
    % 5) mfc_siso_core (use_first_order=true): NEW algebraic 1st-order
    %    estimator (F identified raw, kp applied as explicit feedback).
    %    Plant: ydot = -a0*y + b*u,  alpha = b.  Steady F = -alpha*u_ss.
    % ---------------------------------------------------------------------
    a0 = 1; b = 2; alpha = 2; ref = 1;
    c2 = mfc_siso_core('Ts',Ts,'alpha',alpha,'int_window',5,'time_trajec',10, ...
                       'kp',5,'use_first_order',true);
    [t2,y2,~,F2] = sim1(c2, a0,b, Ts, 6, @(tt) ref*(tt>=0.2));
    fin2 = t2 > 5;
    check('1st (core, algebraic): step tracking (|y-ref|<0.02)', abs(mean(y2(fin2))-ref) < 0.02);
    check('1st (core, algebraic): bounded response',             max(abs(y2)) < 5);
    F2true = -alpha * (a0*ref/b);
    check(sprintf('1st (core, algebraic): steady F=%.3f vs true %.3f',mean(F2(fin2)),F2true), ...
          abs(mean(F2(fin2)) - F2true) < 0.1);

    fprintf('\n%d passed, %d failed\n', N_PASS, N_FAIL);
    if N_FAIL > 0
        error('test_estimators: %d test(s) failed', N_FAIL);
    end
end

% ===== helpers =====
function check(name, cond)
    global N_PASS N_FAIL;
    if cond
        N_PASS = N_PASS + 1; fprintf('  PASS  %s\n', name);
    else
        N_FAIL = N_FAIL + 1; fprintf('  FAIL  %s\n', name);
    end
end

function [t,y,U,F] = sim2(ctrl, a0,a1,b, Ts, Tend, ref)
% Closed loop on yddot = -a1*yd - a0*y + b*u (ZOH on u, forward Euler).
    t = 0:Ts:Tend;  y = 0; yd = 0;
    Y = zeros(size(t)); U = Y; F = Y;
    for k = 1:numel(t)
        [u, f] = step(ctrl, ref(t(k)), y, t(k));  % measure y(t_k), compute u_k
        Y(k) = y; U(k) = u; F(k) = f;
        ydd = -a1*yd - a0*y + b*u;                % apply u_k over [t_k, t_k+Ts]
        yd  = yd + ydd*Ts;
        y   = y  + yd*Ts;
    end
    y = Y;
end

function [t,y,U,F] = sim1(ctrl, a0,b, Ts, Tend, ref)
% Closed loop on ydot = -a0*y + b*u (ZOH on u, forward Euler).
    t = 0:Ts:Tend;  y = 0;
    Y = zeros(size(t)); U = Y; F = Y;
    for k = 1:numel(t)
        [u, f] = step(ctrl, ref(t(k)), y, t(k));
        Y(k) = y; U(k) = u; F(k) = f;
        yd = -a0*y + b*u;
        y  = y + yd*Ts;
    end
    y = Y;
end
