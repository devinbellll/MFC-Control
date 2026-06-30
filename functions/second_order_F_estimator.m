classdef second_order_F_estimator < matlab.System
    % second_order_F_estimator  Second-order Model-Free Control (MFC) SISO
    %   controller built on the SLIDING-WINDOW weighted-integral F estimator
    %   (Eq. 16), wrapped as a full DECOUPLED iPD controller.
    %
    %   This is the decoupled architecture of mfc_siso_decoupled, but F is
    %   estimated with the finite, bounded sliding-window quadrature instead of
    %   the absolute-time algebraic recursion. Two things are kept strictly
    %   separate:
    %
    %     (a) F is estimated from the PURE plant using ONLY (Ym, U, alpha) over
    %         a window of FFilter samples (length T = FFilter*Ts). The error
    %         never enters the estimator, so F is the TRUE plant lumped dynamics
    %         ddot_y - alpha*u.
    %     (b) the stabilizing feedback is an EXPLICIT iPD(I) law on the error:
    %
    %       u = ( -F + ddot_yref ) / alpha  -  ( Kd*dot_e + Kp*e + Ki*∫e ) / alpha
    %
    %   Estimator (Eq. 16), with sigma in [0, T] the window-local time:
    %       F = (60/T^5) * ∫_0^T [ (T^2 - 6*T*σ + 6*σ^2)*y_m(σ)
    %                              - (alpha/2)*σ^2*(T-σ)^2 * u(σ) ] dσ
    %   evaluated by composite SIMPSON quadrature over the window. (Trapezoidal
    %   is NOT usable here: the 60/T^5 prefactor amplifies its O(Ts^2/T^4)
    %   leakage to ~60x error at Ts=0.01, T=0.1. Simpson is exact through cubics,
    %   so the estimate matches mfc_siso_decoupled. The interval count is forced
    %   even and the realized window is N*Ts.)
    %
    %   Ultra-local model: ddot_y = F + alpha*u,   error e = Ym - yref_filter.
    %   With Kd = 2*p, Kp = p^2, Ki = 0 the ideal closed loop is (s + p)^2, the
    %   same double pole at -p as mfc_siso_folded(pole_controller = p).
    %
    %   Inputs : Yref, Ym, t     Outputs: U, F, yref_filter, error
    %   (t is the simulation time, fed from a clock -- as in mfc_siso_run.
    %    It is used only to hold F at 0 until the window has filled, t > T.)

    properties
        alpha    = 1
        WFilter  = 10     % reference pre-filter: tau = WFilter*Ts
        Kp       = 25     % = p^2  for double pole at -p
        Kd       = 10     % = 2*p  for double pole at -p
        Ki       = 0
        Umax     = 600
        Umin     = -600
    end

    properties (Nontunable)
        Ts       = 0.01   % sample time; also fixes the block's discrete rate
        FFilter  = 10     % estimator window length in SAMPLES (window T = FFilter*Ts);
                          % same knob/units as mfc_siso_decoupled.FFilter
    end

    properties (Nontunable, Logical)
        use_control_sat = false
        use_ref_filter  = true
    end

    properties (DiscreteState)
        yref_km1
        yref_km2
        y_buf            % [(N+1)x1] window of measurements, newest last
        u_buf            % [(N+1)x1] window of applied inputs, newest last
        err_km1
        Ierr
        U_km1
    end

    properties (Access = private)
        % Pre-computed window constants (set once in setupImpl).
        N                % number of intervals (forced EVEN for Simpson): N+1 samples
        Tw               % realized window length = N*Ts (>= T, within one Ts)
        sigma            % [(N+1)x1] window-local time nodes, 0 .. Tw
        ym_kernel        % [(N+1)x1] constant measurement kernel
        w                % [(N+1)x1] composite-Simpson weights (1 4 2 .. 4 1)
    end

    methods
        function obj = second_order_F_estimator(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)
        function setupImpl(obj)
            % One-time window quadrature constants. The y_m kernel and nodes
            % depend only on the (nontunable) window T/Ts; the u kernel depends
            % on the tunable alpha and is rebuilt each step.
            %
            % Composite SIMPSON quadrature (not trapezoidal): the 60/T^5
            % prefactor amplifies any quadrature error enormously, and the
            % trapezoidal rule leaves an O(Ts^2/T^4) leakage that makes the
            % estimate unusable at practical Ts (e.g. ~60x error at Ts=0.01,
            % T=0.1). Simpson is exact through cubics, so the constant/linear
            % leakage vanishes and the estimate matches the algebraic variant
            % (mfc_siso_decoupled). Simpson needs an EVEN interval count, so N
            % is rounded up to even and the realized window is Tw = N*Ts.
            obj.N  = obj.FFilter;                         % window length in samples
            obj.N  = obj.N + mod(obj.N, 2);               % force even
            obj.Tw = obj.N * obj.Ts;
            obj.sigma = (0:obj.N).' * obj.Ts;             % [(N+1)x1], 0 .. Tw
            obj.ym_kernel = obj.Tw^2 - 6*obj.Tw*obj.sigma + 6*obj.sigma.^2;
            wv = ones(obj.N+1, 1);  wv(2:2:end-1) = 4;  wv(3:2:end-1) = 2;
            obj.w = wv;
        end

        function [U, F, yref_filter, err] = stepImpl(obj, Yref, Ym, t)
            Ts = obj.Ts;  WF = obj.WFilter;  al = obj.alpha;

            % 1) Reference pre-filter (smooth, twice-differentiable trajectory).
            yref_filter = (Yref + (2*WF^2 + 2*WF)*obj.yref_km1 + (-WF^2)*obj.yref_km2) ...
                          / (WF^2 + 2*WF + 1);
            if ~obj.use_ref_filter
                yref_filter = Yref;
            end
            ddot_yref = (yref_filter - 2*obj.yref_km1 + obj.yref_km2) / Ts^2;

            err = Ym - yref_filter;

            % 2) PURE plant F via sliding-window quadrature (Eq. 16), (Ym, U)
            %    only -- no error term. u uses the previously applied U_km1
            %    (the current U is not known yet), as in the algebraic variants.
            Tw   = obj.Tw;                        % realized (even-interval) window
            ybuf = [obj.y_buf(2:end); Ym];
            ubuf = [obj.u_buf(2:end); obj.U_km1];

            u_kernel  = (al/2) .* obj.sigma.^2 .* (Tw - obj.sigma).^2;
            integrand = obj.ym_kernel .* ybuf - u_kernel .* ubuf;

            F = 0;
            if t > Tw                             % hold until the window fills
                F = (60 / Tw^5) * (Ts/3) * sum(obj.w .* integrand);   % Simpson
            end

            % 3) Explicit iPD(I) feedback on the error (no estimator lag).
            dot_err = (err - obj.err_km1) / Ts;                % noise-sensitive term
            Ierr    = obj.Ierr + (err + obj.err_km1)/2 * Ts;   % trapezoidal integral

            U = (-F + ddot_yref)/al - (obj.Kd*dot_err + obj.Kp*err + obj.Ki*Ierr)/al;

            if obj.use_control_sat
                U_sat = max(obj.Umin, min(obj.Umax, U));
                if U_sat ~= U          % saturated -> freeze integrator (anti-windup)
                    Ierr = obj.Ierr;
                end
                U = U_sat;
            end

            % 4) Advance discrete state.
            obj.yref_km2 = obj.yref_km1;  obj.yref_km1 = yref_filter;
            obj.y_buf    = ybuf;          obj.u_buf    = ubuf;
            obj.err_km1  = err;
            obj.Ierr     = Ierr;
            obj.U_km1    = U;
        end

        function resetImpl(obj)
            n = obj.FFilter;  n = n + mod(n, 2);   % window samples, even (Simpson)
            obj.yref_km1 = 0;  obj.yref_km2 = 0;
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
                    n = obj.FFilter;  n = n + mod(n, 2);
                    sz = [n + 1, 1];
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
