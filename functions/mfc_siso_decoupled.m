classdef mfc_siso_decoupled < matlab.System
    % mfc_siso_decoupled  Second-order Model-Free Control (MFC) SISO controller.
    %
    %   DECOUPLED variant (canonical iPD; second-order analogue of the
    %   algebraic first-order F estimator, e.g. mfc_siso_core with
    %   use_first_order = true). Two things are kept strictly separate:
    %
    %     (a) F is estimated from the PURE plant model using ONLY (Ym, U,
    %         alpha, t). The error never enters the estimator, so F is the TRUE
    %         plant lumped dynamics ddot_y - alpha*u.
    %     (b) the stabilizing feedback is an EXPLICIT PD(I) law on the error:
    %
    %       u = ( -F + ddot_yref ) / alpha  -  ( Kd*dot_e + Kp*e + Ki*∫e ) / alpha
    %
    %   The feedback acts on the current error (no estimator lag), so gains map
    %   cleanly onto the closed-loop poles. Trade-off vs. mfc_siso_core (the
    %   folded, use_first_order=false variant): the PD term needs an explicit
    %   error derivative (finite difference) which amplifies measurement noise.
    %
    %   EQUIVALENCE: Kd = 2*p, Kp = p^2, Ki = 0 gives ddot_e + 2p*dot_e +
    %   p^2*e = 0, i.e. (s + p)^2 -- the same double pole as
    %   mfc_siso_core(kp = p, use_kd = false).
    %
    %   Inputs : Yref, Ym, t     Outputs: U, F, yref_filter, error
    %   (t is the simulation time, fed from a clock -- as in mfc_siso_run.)

    properties
        alpha    = 1
        FFilter  = 10
        WFilter  = 10
        Kp       = 25     % = p^2  for double pole at -p
        Kd       = 10     % = 2*p  for double pole at -p
        Ki       = 0
        Umax     = 600
        Umin     = -600
    end

    properties (Nontunable)
        Ts       = 0.01   % sample time; also fixes the block's discrete rate
    end

    properties (Nontunable, Logical)
        use_control_sat = false
        use_ref_filter  = true
    end

    properties (DiscreteState)
        yref_km1
        yref_km2
        Ym_km1
        Ym_km2
        Fnum_km1
        Fnum_km2
        Fden_km1
        Fden_km2
        err_km1
        Ierr
        U_km1
    end

    methods
        function obj = mfc_siso_decoupled(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)
        function [U, F, yref_filter, err] = stepImpl(obj, Yref, Ym, t)
            Ts = obj.Ts;  WF = obj.WFilter;  FF = obj.FFilter;
            al = obj.alpha;

            % 1) Reference pre-filter (smooth, twice-differentiable trajectory).
            yref_filter = (Yref + (2*WF^2 + 2*WF)*obj.yref_km1 + (-WF^2)*obj.yref_km2) ...
                          / (WF^2 + 2*WF + 1);
            if ~obj.use_ref_filter
                yref_filter = Yref;
            end
            ddot_yref = (yref_filter - 2*obj.yref_km1 + obj.yref_km2) / Ts^2;

            err = Ym - yref_filter;

            % 2) PURE plant F estimate from (Ym, U, alpha, t) ONLY -- no error term.
            sdy   = -(t*Ym - (t-Ts)*obj.Ym_km1) / Ts;
            s2d2y =  (t^2*Ym - 2*(t-Ts)^2*obj.Ym_km1 + (t-2*Ts)^2*obj.Ym_km2) / Ts^2;
            d2u   =  t^2*obj.U_km1;

            num = 2*Ym + 4*sdy + s2d2y - al*d2u;
            den = t^2;

            Fnum = (num + (2*FF^2 + 2*FF)*obj.Fnum_km1 + (-FF^2)*obj.Fnum_km2) ...
                   / (FF^2 + 2*FF + 1);
            Fden = (den + (2*FF^2 + 2*FF)*obj.Fden_km1 + (-FF^2)*obj.Fden_km2) ...
                   / (FF^2 + 2*FF + 1);

            F = 0;
            if (Fden ~= 0) && (t > 0.1)
                F = Fnum / Fden;
            end

            % 3) Explicit iPD(I) feedback on the error (no estimator lag).
            dot_err = (err - obj.err_km1) / Ts;                 % noise-sensitive term
            Ierr    = obj.Ierr + (err + obj.err_km1)/2 * Ts;    % trapezoidal integral

            U = (-F + ddot_yref)/al - (obj.Kd*dot_err + obj.Kp*err + obj.Ki*Ierr)/al;

            if obj.use_control_sat
                U_sat = max(obj.Umin, min(obj.Umax, U));
                if U_sat ~= U          % saturated -> freeze integrator (anti-windup)
                    Ierr = obj.Ierr;
                end
                U = U_sat;
            end

            obj.yref_km2 = obj.yref_km1;  obj.yref_km1 = yref_filter;
            obj.Ym_km2   = obj.Ym_km1;    obj.Ym_km1   = Ym;
            obj.Fnum_km2 = obj.Fnum_km1;  obj.Fnum_km1 = Fnum;
            obj.Fden_km2 = obj.Fden_km1;  obj.Fden_km1 = Fden;
            obj.err_km1  = err;
            obj.Ierr     = Ierr;
            obj.U_km1    = U;
        end

        function resetImpl(obj)
            obj.yref_km1 = 0;  obj.yref_km2 = 0;
            obj.Ym_km1   = 0;  obj.Ym_km2   = 0;
            obj.Fnum_km1 = 0;  obj.Fnum_km2 = 0;
            obj.Fden_km1 = 0;  obj.Fden_km2 = 0;
            obj.err_km1  = 0;
            obj.Ierr     = 0;
            obj.U_km1    = 0;
        end

        function [sz, dt, cp] = getDiscreteStateSpecificationImpl(~, ~)
            % All discrete states are scalar, real, double.
            sz = [1 1];
            dt = 'double';
            cp = false;
        end

        function sts = getSampleTimeImpl(obj)
            % Run at a fixed discrete rate Ts (do NOT inherit), so the
            % per-sample recursions (ref filter, F filter) advance at the same
            % rate as the original delay-based subsystem.
            sts = createSampleTime(obj, 'Type', 'Discrete', 'SampleTime', obj.Ts);
        end

        % All four outputs (U, F, yref_filter, error) are scalar real double.
        function [o1, o2, o3, o4] = getOutputSizeImpl(~)
            o1 = [1 1];  o2 = [1 1];  o3 = [1 1];  o4 = [1 1];
        end
        function [o1, o2, o3, o4] = getOutputDataTypeImpl(~)
            o1 = 'double';  o2 = 'double';  o3 = 'double';  o4 = 'double';
        end
        function [o1, o2, o3, o4] = isOutputComplexImpl(~)
            o1 = false;  o2 = false;  o3 = false;  o4 = false;
        end
        function [o1, o2, o3, o4] = isOutputFixedSizeImpl(~)
            o1 = true;  o2 = true;  o3 = true;  o4 = true;
        end
    end
end
