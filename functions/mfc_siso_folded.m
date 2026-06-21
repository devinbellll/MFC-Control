classdef mfc_siso_folded < matlab.System
    % mfc_siso_folded  Second-order Model-Free Control (MFC) SISO controller.
    %
    %   FOLDED variant: the auxiliary stabilizing feedback (the desired
    %   closed-loop "control poles") is folded INTO the algebraic estimate of
    %   F, so the command carries only F and the reference acceleration:
    %
    %       u = ( -F + ddot_yref ) / alpha
    %
    %   The pole terms a*dot_e + b*e enter the estimator numerator, so the
    %   PD-like action is realized through the same noise-robust algebraic
    %   integrals as F (no explicit error differentiation). Price: the feedback
    %   inherits the estimator's lag, coupling tuning to FFilter. Compare with
    %   mfc_siso_decoupled. Direct port of functions/mfc_siso_run.m.
    %
    %   Ultra-local model:  ddot_y = F + alpha*u,   error e = Ym - yref_filter.
    %   With a = -2*p, b = -p^2 the ideal closed loop is (s + p)^2, i.e. a
    %   DOUBLE POLE AT -pole_controller. Pass pole_controller > 0.
    %
    %   Inputs : Yref, Ym, t     Outputs: U, F, yref_filter, error
    %   (t is the simulation time, fed from a clock -- as in mfc_siso_run.)

    properties
        alpha           = 10
        FFilter         = 5
        WFilter         = 250
        pole_controller = 1      % positive magnitude; double pole at -pole_controller
        Umax            = 600
        Umin            = -600
    end

    properties (Nontunable)
        Ts              = 0.001  % sample time; also fixes the block's discrete rate
    end

    properties (Nontunable, Logical)
        use_control_sat = false
        use_ref_filter  = true
    end

    properties (DiscreteState)
        yref_km1
        yref_km2
        Fnum_km1
        Fnum_km2
        Fden_km1
        Fden_km2
        err_km1
        err_km2
        U_km1
    end

    methods
        function obj = mfc_siso_folded(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)
        function [U, F, yref_filter, err] = stepImpl(obj, Yref, Ym, t)
            Ts = obj.Ts;  WF = obj.WFilter;  FF = obj.FFilter;
            al = obj.alpha;  p = obj.pole_controller;

            % 1) Reference pre-filter (smooth, twice-differentiable trajectory).
            yref_filter = (Yref + (2*WF^2 + 2*WF)*obj.yref_km1 + (-WF^2)*obj.yref_km2) ...
                          / (WF^2 + 2*WF + 1);
            if ~obj.use_ref_filter
                yref_filter = Yref;
            end
            ddot_yref = (yref_filter - 2*obj.yref_km1 + obj.yref_km2) / Ts^2;

            % 2) Control poles FOLDED into the estimator: (s+p)^2 -> a=-2p, b=-p^2.
            a = -2*p;
            b = -p^2;
            err = Ym - yref_filter;

            sde   = -(t*err - (t-Ts)*obj.err_km1) / Ts;
            s2d2e =  (t^2*err - 2*(t-Ts)^2*obj.err_km1 + (t-2*Ts)^2*obj.err_km2) / Ts^2;
            sd2e  =  (t^2*err - (t-Ts)^2*obj.err_km1) / Ts;
            de    = -t*err;
            d2e   =  t^2*err;
            d2u   =  t^2*obj.U_km1;

            num = 2*err + 4*sde + s2d2e - a*(2*de + sd2e) - b*d2e - al*d2u;
            den = t^2;

            Fnum = (num + (2*FF^2 + 2*FF)*obj.Fnum_km1 + (-FF^2)*obj.Fnum_km2) ...
                   / (FF^2 + 2*FF + 1);
            Fden = (den + (2*FF^2 + 2*FF)*obj.Fden_km1 + (-FF^2)*obj.Fden_km2) ...
                   / (FF^2 + 2*FF + 1);

            F = 0;
            if (Fden ~= 0) && (t > 0.1)
                F = Fnum / Fden;
            end

            U = -F/al + ddot_yref/al;
            if obj.use_control_sat
                U = max(obj.Umin, min(obj.Umax, U));
            end

            obj.yref_km2 = obj.yref_km1;  
            obj.yref_km1 = yref_filter;
            obj.err_km2  = obj.err_km1;   
            obj.err_km1  = err;
            obj.Fnum_km2 = obj.Fnum_km1;  
            obj.Fnum_km1 = Fnum;
            obj.Fden_km2 = obj.Fden_km1;  
            obj.Fden_km1 = Fden;
            obj.U_km1    = U;
        end

        function resetImpl(obj)
            obj.yref_km1 = 0;  obj.yref_km2 = 0;
            obj.Fnum_km1 = 0;  obj.Fnum_km2 = 0;
            obj.Fden_km1 = 0;  obj.Fden_km2 = 0;
            obj.err_km1  = 0;  obj.err_km2  = 0;
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
