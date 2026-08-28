classdef mfc_fhat_riachy2_mimo_block < matlab.System
    % mfc_fhat_riachy2_mimo_block  F estimator, 2nd order, RIACHY'S TRICK,
    % with MATRIX gains.
    %
    %   The vector-valued twin of mfc_fhat_riachy2_block: estimates
    %
    %       Fk = F + Kd*dot_y
    %
    %   in the second-order ultra-local model rewritten on the auxiliary
    %   output Y (Riachy et al.), with y, u, F and Fk n-by-1 and both Kd and
    %   alpha square n-by-n:
    %
    %       Y = y + Kd*int y      =>      ddot_Y = Fk + alpha*u
    %
    %   The rewrite is the SISO one with matrix products: adding Kd*dot_y to
    %   both sides of ddot_y = F + alpha*u is a linear operation, so it goes
    %   through unchanged for a matrix Kd. Kd need not be diagonal -- a full
    %   Kd folds cross-channel derivative feedback into Y, and the closed
    %   loop is the matrix polynomial ddot_e + Kd*dot_e + Kp*e = 0.
    %
    %   DECOUPLED, in the same sense as the SISO block: Fk is a property of
    %   the plant, nothing about the setpoint or the error enters here. The
    %   remaining feedback is explicit and it is a PI, NOT a PID:
    %
    %       fb <- Kp*err + Ki*int err, D = 0, on err = y - sp_filt
    %       ff <- ddot_sp + Kd*dot_sp
    %
    %   so that mfc_command_mimo_block computes
    %
    %       u = alpha \ ( -Fk + ddot_sp + Kd*dot_sp - Kp*err - Ki*int err )
    %
    %   Putting D on the feedback applies Kd twice; dropping Kd*dot_sp from
    %   ff leaves the feedforward inconsistent with Y and costs tracking on a
    %   moving setpoint. Both mistakes are silent.
    %
    %   ESTIMATOR TYPE IS A PARAMETER, as in the SISO block: algebraic
    %   (growing window, needs a run-start clock, weights the unbounded
    %   int y by t^2) or sliding window (FIR, fixed memory, no internal
    %   smoothing). Both act channel-wise on Y through ONE shared window --
    %   the algebraic denominator t^2 is scalar, and the sliding window's
    %   taps are scalar quadrature weights. Prefer the sliding window here,
    %   for the reason in Knowledge/riachy-trick.md.
    %
    %   Ports
    %     In : y       plant measurement                       (n-by-1)
    %          u_prev  command that actually reached the plant over the LAST
    %                  sample, through an explicit unit delay  (n-by-1)
    %          t       clock, starting with the run            (scalar)
    %          alpha   optional live gain                      (n-by-n)
    %     Out: F_hat   the estimate of Fk = F + Kd*dot_y        (n-by-1)
    %          Y       optional, the auxiliary output (logging / debug)
    %
    %   The math is mfc_riachy_transform (matrix Kd) followed by
    %   mfc_fhat_algebraic_second_order (a_fold = b_fold = 0) or
    %   mfc_fhat_sliding_window; this class only maps parameters and
    %   Simulink state onto them.
    %
    %   See also mfc_fhat_riachy2_block, mfc_riachy_transform,
    %   mfc_fhat_alg2_decoupled_mimo_block, mfc_command_mimo_block.

    properties
        % alpha Ultra-local model input gain, square n-by-n (ignored if the live alpha input is enabled)
        alpha = eye(2)
        % Kd Derivative gain folded into Y, square n-by-n (the Kd of ddot_e + Kd*dot_e + Kp*e); keep D = 0 on the external feedback
        Kd = 10*eye(2)
    end

    properties (Nontunable)
        % n Number of channels (alpha and Kd are n-by-n, signals are n-by-1)
        n = 2
        % estimator F estimator applied to the auxiliary output Y
        estimator = 'Algebraic (growing window)'
        % Ts Sample time [s] (fixes the block's discrete rate)
        Ts = 0.01
        % est_filter_window Algebraic only: internal num/den smoother memory [samples]
        est_filter_window = 10
        % est_hold_time Algebraic only: F_hat held at zero until t exceeds this [s]
        est_hold_time = 0.1
        % window_samples Sliding window only: window length [intervals]; Tw = window_samples*Ts
        window_samples = 10
    end

    properties (Hidden, Constant)
        estimatorSet = matlab.system.StringSet({ ...
            'Algebraic (growing window)', ...
            'Sliding window (FIR)'});
    end

    properties (Nontunable, Logical)
        % use_live_alpha Add the alpha input port (overrides the alpha parameter)
        use_live_alpha = false
        % output_Y Add the Y output port (the auxiliary output, for logging)
        output_Y = false
    end

    properties (DiscreteState)
        % Riachy transform
        int_km1
        y_km1
        % algebraic estimator
        z_km1
        z_km2
        num_filt_km1
        num_filt_km2
        den_filt_km1
        den_filt_km2
        % sliding-window estimator (1-row placeholders when algebraic)
        y_buf
        u_buf
    end

    properties (Access = private)
        kernel     % quadrature kernel built by mfc_siso.window_kernel in setupImpl
    end

    methods
        function obj = mfc_fhat_riachy2_mimo_block(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = private)
        function tf = isAlgebraic(obj)
            tf = strncmp(obj.estimator, 'Algebraic', 9);
        end

        function n_buf = bufferLength(obj)
            % 1-row placeholder when algebraic: the (dead) sliding-window
            % branch is still compiled under code generation, so the
            % buffers must exist and be typed either way. Built only from
            % Nontunable properties, so it is a compile-time constant.
            if isAlgebraic(obj)
                n_buf = 1;
            else
                n_buf = obj.window_samples + 1;
            end
        end
    end

    methods (Access = protected)

        function setupImpl(obj)
            % Always built: harmless when algebraic, and keeping it
            % unconditional keeps the kernel a fixed type under codegen.
            obj.kernel = mfc_siso.window_kernel(2, obj.window_samples, obj.Ts);
        end

        function varargout = stepImpl(obj, y, u_prev, t, varargin)
            if obj.use_live_alpha
                alpha_k = varargin{1};
            else
                alpha_k = obj.alpha;
            end

            % 1. Riachy's auxiliary output: Y = y + Kd*int y (matrix Kd)
            rc = struct('int_km1', obj.int_km1, 'y_km1', obj.y_km1);
            [Y, rc] = mfc_riachy_transform(y, obj.Kd, obj.Ts, rc);
            obj.int_km1 = rc.int_km1;
            obj.y_km1   = rc.y_km1;

            % 2. A standard 2nd-order estimator, driven by Y instead of y.
            %    Nothing is folded (a_fold = b_fold = 0): the Kd term is
            %    already inside Y, and it comes back out inside F_hat.
            est = struct( ...
                'z_km1',        obj.z_km1, ...
                'z_km2',        obj.z_km2, ...
                'num_filt_km1', obj.num_filt_km1, ...
                'num_filt_km2', obj.num_filt_km2, ...
                'den_filt_km1', obj.den_filt_km1, ...
                'den_filt_km2', obj.den_filt_km2, ...
                'y_buf',        obj.y_buf, ...
                'u_buf',        obj.u_buf);

            if isAlgebraic(obj)
                [F_hat, est] = mfc_fhat_algebraic_second_order( ...
                    Y, u_prev, alpha_k, t, obj.Ts, obj.est_filter_window, ...
                    obj.est_hold_time, 0, 0, est);
            else
                [F_hat, est] = mfc_fhat_sliding_window( ...
                    Y, u_prev, alpha_k, t, obj.kernel, est);
            end

            obj.z_km1        = est.z_km1;
            obj.z_km2        = est.z_km2;
            obj.num_filt_km1 = est.num_filt_km1;
            obj.num_filt_km2 = est.num_filt_km2;
            obj.den_filt_km1 = est.den_filt_km1;
            obj.den_filt_km2 = est.den_filt_km2;
            obj.y_buf        = est.y_buf;
            obj.u_buf        = est.u_buf;

            varargout{1} = F_hat;
            if obj.output_Y
                varargout{2} = Y;
            end
        end

        function resetImpl(obj)
            obj.int_km1      = zeros(obj.n, 1);
            obj.y_km1        = zeros(obj.n, 1);
            obj.z_km1        = zeros(obj.n, 1);
            obj.z_km2        = zeros(obj.n, 1);
            obj.num_filt_km1 = zeros(obj.n, 1);
            obj.num_filt_km2 = zeros(obj.n, 1);
            obj.den_filt_km1 = 0;   % denominator is t^2: scalar, shared
            obj.den_filt_km2 = 0;
            % One COLUMN per channel, with a codegen-constant length: code
            % generation types the discrete states from these assignments.
            n_buf     = bufferLength(obj);
            obj.y_buf = zeros(n_buf, obj.n);
            obj.u_buf = zeros(n_buf, obj.n);
        end

        function validatePropertiesImpl(obj)
            if ~obj.use_live_alpha && ~isequal(size(obj.alpha), [obj.n obj.n])
                error('mfc:mimo:alphaSize', ...
                    'alpha must be %d-by-%d to match n.', obj.n, obj.n);
            end
            if ~isequal(size(obj.Kd), [obj.n obj.n])
                error('mfc:mimo:KdSize', ...
                    'Kd must be %d-by-%d to match n.', obj.n, obj.n);
            end
        end

        % ---- save/load of the locked object ----------------------------
        % kernel is PRIVATE and is built in setupImpl, so the base class
        % does not carry it through a save/load of a LOCKED object. Simulink
        % saves and clones locked System objects for fast restart and for
        % array sim(), and without these two methods kernel comes back
        % empty on every run after the first -- which surfaces as
        % "simulations completed with errors at indices [2 3]" while
        % individual 1 succeeds, i.e. it reads as a bad candidate rather than
        % a broken block. Measured 2026-08-21 driving a GA over this block.
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            if isLocked(obj)
                s.kernel = obj.kernel;
            end
        end

        function loadObjectImpl(obj, s, wasLocked)
            if wasLocked
                obj.kernel = s.kernel;
            end
            loadObjectImpl@matlab.System(obj, s, wasLocked);
        end

        function [sz, dt, cp] = getDiscreteStateSpecificationImpl(obj, name)
            switch name
                case {'y_buf', 'u_buf'}
                    sz = [bufferLength(obj), obj.n];
                case {'den_filt_km1', 'den_filt_km2'}
                    sz = [1 1];
                otherwise
                    sz = [obj.n 1];
            end
            dt = 'double';  cp = false;
        end

        function sts = getSampleTimeImpl(obj)
            % Fixed discrete rate (do NOT inherit): the integrator, the
            % backward-difference recursion and the window buffers all
            % assume they advance exactly once per Ts.
            sts = createSampleTime(obj, 'Type', 'Discrete', 'SampleTime', obj.Ts);
        end

        % ---- Ports -------------------------------------------------------
        function num = getNumInputsImpl(obj), num = 3 + obj.use_live_alpha; end
        function varargout = getInputNamesImpl(obj)
            names = {'y', 'u_prev', 't'};
            if obj.use_live_alpha, names{end+1} = 'alpha'; end
            varargout = names;
        end

        function num = getNumOutputsImpl(obj), num = 1 + obj.output_Y; end
        function varargout = getOutputNamesImpl(obj)
            names = {'F_hat'};
            if obj.output_Y, names{end+1} = 'Y'; end
            varargout = names;
        end
        function varargout = getOutputSizeImpl(obj)
            [varargout{1:getNumOutputsImpl(obj)}] = deal([obj.n 1]);
        end
        function varargout = getOutputDataTypeImpl(obj)
            [varargout{1:getNumOutputsImpl(obj)}] = deal('double');
        end
        function varargout = isOutputComplexImpl(obj)
            [varargout{1:getNumOutputsImpl(obj)}] = deal(false);
        end
        function varargout = isOutputFixedSizeImpl(obj)
            [varargout{1:getNumOutputsImpl(obj)}] = deal(true);
        end

        function icon = getIconImpl(obj)
            if isAlgebraic(obj)
                kind = 'algebraic';
            else
                kind = sprintf('window, Tw = %g s', obj.window_samples*obj.Ts);
            end
            icon = sprintf('F-hat Riachy 2nd\nY = y + Kd*int y\nmatrix Kd, alpha (%dx%d)\n%s', ...
                obj.n, obj.n, kind);
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_fhat_riachy2_mimo_block', ...
                'Title', 'MFC F-hat: Riachy trick, 2nd order, matrix gains', ...
                'Text', sprintf(['Estimates Fk = F + Kd*dot_y from the auxiliary ', ...
                    'output Y = y + Kd*int y, for which ddot_Y = Fk + alpha*u, with ', ...
                    'n-by-1 signals and square n-by-n Kd and alpha.\n\n', ...
                    'The derivative term comes back inside F_hat, so the loop needs no ', ...
                    'derivative estimate at all. Wire the remaining feedback with P = Kp, ', ...
                    'I = Ki and D = 0 on err = y - sp_filt, and feed ', ...
                    'mfc_command_mimo_block ff = ddot_sp + Kd*dot_sp. Putting D on that ', ...
                    'feedback applies Kd twice.\n\nThe estimator applied to Y is a ', ...
                    'parameter, not a separate block: algebraic (growing window, needs a ', ...
                    'run-start clock) or sliding window (FIR, fixed memory, no internal ', ...
                    'smoothing). int y is unbounded -- with a non-zero steady state the ', ...
                    'sliding window is the better-behaved choice.']));
        end
    end
end
