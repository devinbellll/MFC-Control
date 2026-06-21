classdef first_order_F_estimator < matlab.System
    % first_order_F_estimator  First-order Model-Free Control (MFC) SISO
    %   controller built on the SLIDING-WINDOW weighted-integral F estimator
    %   (Eq. 11), wrapped as a full DECOUPLED iP(I) controller.
    %
    %   First-order analogue of second_order_F_estimator / mfc_siso_decoupled
    %   (and the System-object form of CSM_NG_alphaV). Two things are kept
    %   strictly separate:
    %
    %     (a) F is estimated from the PURE plant using ONLY (Ym, U, alpha) over
    %         a window of length T. The error never enters the estimator, so F
    %         is the TRUE plant lumped dynamics dot_y - alpha*u.
    %     (b) the stabilizing feedback is an EXPLICIT iP(I) law on the error:
    %
    %       u = ( -F + dot_yref ) / alpha  -  ( Kp*e + Ki*∫e ) / alpha
    %
    %   A first-order ultra-local model dot_y = F + alpha*u admits a single
    %   closed-loop pole, so there is NO derivative term (no Kd): with Ki = 0,
    %   dot_e + Kp*e = 0, i.e. a single pole at -Kp.
    %
    %   Estimator (Eq. 11), with sigma in [0, T] the window-local time:
    %       F = -(6/T^3) * ∫_0^T [ (T - 2*σ)*y_m(σ)
    %                              + alpha*σ*(T-σ) * u(σ) ] dσ
    %   evaluated by the trapezoidal rule over the (N+1)-sample window.
    %
    %   Ultra-local model: dot_y = F + alpha*u,   error e = Ym - yref_filter.
    %
    %   Inputs : Yref, Ym, t     Outputs: U, F, yref_filter, error
    %   (t is the simulation time, fed from a clock -- as in mfc_siso_run.
    %    It is used only to hold F at 0 until the window has filled, t > T.)

    properties
        alpha    = 1
        WFilter  = 10     % reference pre-filter: tau = WFilter*Ts (first order)
        Kp       = 5      % single closed-loop pole at -Kp (with Ki = 0)
        Ki       = 0
        Umax     = 600
        Umin     = -600
    end

    properties (Nontunable)
        Ts       = 0.01   % sample time; also fixes the block's discrete rate
        T        = 0.1    % estimator sliding-window length (s); sets buffer size
    end

    properties (Nontunable, Logical)
        use_control_sat = false
        use_ref_filter  = true
    end

    properties (DiscreteState)
        yref_km1
        y_buf            % [(N+1)x1] window of measurements, newest last
        u_buf            % [(N+1)x1] window of applied inputs, newest last
        err_km1
        Ierr
        U_km1
    end

    properties (Access = private)
        % Pre-computed window constants (set once in setupImpl).
        N                % number of intervals: window holds N+1 samples
        sigma            % [(N+1)x1] window-local time nodes, 0 .. T
        ym_kernel        % [(N+1)x1] constant measurement kernel (T - 2*sigma)
        w                % [(N+1)x1] trapezoidal weights
    end

    methods
        function obj = first_order_F_estimator(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)
        function setupImpl(obj)
            % One-time window quadrature constants. The y_m kernel and nodes
            % depend only on the (nontunable) window T/Ts; the u kernel depends
            % on the tunable alpha and is rebuilt each step.
            obj.N     = round(obj.T / obj.Ts);
            obj.sigma = (0:obj.N).' * obj.Ts;             % [(N+1)x1], 0 .. T
            obj.ym_kernel = obj.T - 2*obj.sigma;
            wv = ones(obj.N+1, 1);  wv(1) = 0.5;  wv(end) = 0.5;
            obj.w = wv;
        end

        function [U, F, yref_filter, err] = stepImpl(obj, Yref, Ym, t)
            Ts = obj.Ts;  WF = obj.WFilter;  al = obj.alpha;  T = obj.T;

            % 1) First-order reference pre-filter (tau = WF*Ts), unity DC gain.
            yref_filter = (Yref + WF*obj.yref_km1) / (WF + 1);
            if ~obj.use_ref_filter
                yref_filter = Yref;
            end
            dot_yref = (yref_filter - obj.yref_km1) / Ts;

            err = Ym - yref_filter;

            % 2) PURE plant F via sliding-window quadrature (Eq. 11), (Ym, U)
            %    only -- no error term. u uses the previously applied U_km1
            %    (the current U is not known yet), as in the algebraic variants.
            ybuf = [obj.y_buf(2:end); Ym];
            ubuf = [obj.u_buf(2:end); obj.U_km1];

            u_kernel  = al .* obj.sigma .* (T - obj.sigma);
            integrand = obj.ym_kernel .* ybuf + u_kernel .* ubuf;

            F = 0;
            if t > T                              % hold until the window fills
                F = (-6 / T^3) * Ts * sum(obj.w .* integrand);
            end

            % 3) Explicit iP(I) feedback on the error (no estimator lag).
            Ierr = obj.Ierr + (err + obj.err_km1)/2 * Ts;      % trapezoidal integral

            U = (-F + dot_yref)/al - (obj.Kp*err + obj.Ki*Ierr)/al;

            if obj.use_control_sat
                U_sat = max(obj.Umin, min(obj.Umax, U));
                if U_sat ~= U          % saturated -> freeze integrator (anti-windup)
                    Ierr = obj.Ierr;
                end
                U = U_sat;
            end

            % 4) Advance discrete state.
            obj.yref_km1 = yref_filter;
            obj.y_buf    = ybuf;          obj.u_buf    = ubuf;
            obj.err_km1  = err;
            obj.Ierr     = Ierr;
            obj.U_km1    = U;
        end

        function resetImpl(obj)
            n = round(obj.T / obj.Ts);
            obj.yref_km1 = 0;
            obj.y_buf    = zeros(n+1, 1);
            obj.u_buf    = zeros(n+1, 1);
            obj.err_km1  = 0;
            obj.Ierr     = 0;
            obj.U_km1    = 0;
        end

        function [sz, dt, cp] = getDiscreteStateSpecificationImpl(obj, name)
            % Window buffers are [(N+1)x1]; everything else is scalar. All real
            % double.
            switch name
                case {'y_buf', 'u_buf'}
                    sz = [round(obj.T / obj.Ts) + 1, 1];
                otherwise
                    sz = [1 1];
            end
            dt = 'double';
            cp = false;
        end

        function sts = getSampleTimeImpl(obj)
            % Run at a fixed discrete rate Ts (do NOT inherit), so the
            % per-sample recursions (ref filter, sliding window) advance at the
            % same rate as the original delay-based subsystem.
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
