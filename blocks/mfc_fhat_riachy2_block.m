classdef mfc_fhat_riachy2_block < matlab.System
    % mfc_fhat_riachy2_block  F estimator, 2nd order, RIACHY'S TRICK.
    %
    %   Estimates
    %
    %       Fk = F + Kd*dot_y
    %
    %   in the second-order ultra-local model rewritten on the auxiliary
    %   output Y (Riachy et al.):
    %
    %       Y = y + Kd*int y      =>      ddot_Y = Fk + alpha*u
    %
    %   The block integrates the measurement, forms Y, and runs a standard
    %   second-order F estimator on it. Because Fk already carries Kd*dot_y,
    %   the loop it closes is an iPD (or iPID with Ki) that never estimates
    %   a derivative of the measurement -- that is the whole point of the
    %   trick. See mfc_riachy_transform for the derivation.
    %
    %   DECOUPLED. Fk is a property of the plant, not of the loop: nothing
    %   about the setpoint or the error enters here. The remaining feedback
    %   is explicit and it is a PI, NOT a PID:
    %
    %       fb <- stock Discrete PID with P = Kp, I = Ki, D = 0, driven by
    %             err = y - sp_filt
    %       ff <- ddot_sp + Kd*dot_sp    (a stock Gain and Sum on the
    %             smoother's derivative outputs)
    %
    %   so that mfc_command_block computes
    %
    %       u = ( -Fk + ddot_sp + Kd*dot_sp - Kp*err - Ki*int err ) / alpha
    %
    %   and the closed loop is ddot_err + Kd*dot_err + Kp*err = 0. The two
    %   ways to get Kd wrong are worth stating: putting D on the PID applies
    %   Kd twice, and dropping Kd*dot_sp from ff leaves the reference
    %   feedforward inconsistent with Y and costs tracking on a moving
    %   setpoint (it is exact for a settled step).
    %
    %   ESTIMATOR TYPE IS A PARAMETER HERE, unlike coupled/decoupled and
    %   model order, which are separate blocks. Both choices estimate the
    %   same Fk from the same Y with the same ports and the same wiring --
    %   only the numerics differ, so it is a tuning knob, not a structure:
    %
    %     Algebraic (growing window)  operational-calculus recursion; the
    %         window grows from t = 0, so t must be a clock that starts with
    %         the run, and est_hold_time guards the near-zero denominator.
    %         Smooths its own numerator and denominator.
    %     Sliding window (Simpson)    fixed-length quadrature over the last
    %         Tw = window_samples*Ts seconds (rounded up to even), held at
    %         zero until the window fills. Insensitive to the time origin,
    %         and with no internal smoothing at all -- follow it with an
    %         mfc_smoother_block if Fk is noisy.
    %
    %   The growing-window caveat bites harder here than on a plain
    %   estimator: int y is unbounded, so with a non-zero steady state Y
    %   ramps forever and the algebraic estimator weights that ramp by t^2.
    %   The sliding window has finite memory and does not care. See
    %   Knowledge/riachy-trick.md.
    %
    %   Ports
    %     In : y       plant measurement
    %          u_prev  command that actually reached the plant over the LAST
    %                  sample; through an explicit unit delay
    %          t       clock, starting with the run
    %          alpha   optional live gain (overrides the mask parameter)
    %     Out: F_hat   the estimate of Fk = F + Kd*dot_y
    %          Y       optional, the auxiliary output (logging / debug)
    %
    %   The math is mfc_riachy_transform followed by
    %   mfc_fhat_algebraic_second_order (a_fold = b_fold = 0) or
    %   mfc_fhat_sliding_window; this class only maps parameters and
    %   Simulink state onto them.
    %
    %   See also mfc_riachy_transform, mfc_fhat_alg2_decoupled_block,
    %   mfc_fhat_window_block, mfc_command_block, mfc_smoother_block.

    properties
        % alpha Ultra-local model input gain (ignored if the live alpha input is enabled)
        alpha = 1
        % Kd Derivative gain folded into Y (the Kd of s^2 + Kd*s + Kp); keep D = 0 on the external PID
        Kd = 10
    end

    properties (Nontunable)
        % estimator F estimator applied to the auxiliary output Y
        estimator = 'Algebraic (growing window)'
        % Ts Sample time [s] (fixes the block's discrete rate)
        Ts = 0.01
        % est_filter_window Algebraic only: internal num/den smoother memory [samples]
        est_filter_window = 10
        % est_hold_time Algebraic only: F_hat held at zero until t exceeds this [s]
        est_hold_time = 0.1
        % window_samples Sliding window only: window length [samples], rounded up to even for Simpson
        window_samples = 10
    end

    properties (Hidden, Constant)
        estimatorSet = matlab.system.StringSet({ ...
            'Algebraic (growing window)', ...
            'Sliding window (Simpson)'});
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
        % sliding-window estimator (1x1 placeholders when algebraic)
        y_buf
        u_buf
    end

    properties (Access = private)
        kernel     % quadrature kernel built by mfc_siso.window_kernel in setupImpl
    end

    methods
        function obj = mfc_fhat_riachy2_block(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = private)
        function tf = isAlgebraic(obj)
            tf = strncmp(obj.estimator, 'Algebraic', 9);
        end

        function n_buf = bufferLength(obj)
            % 1x1 placeholder when algebraic: the (dead) sliding-window
            % branch is still compiled under code generation, so the
            % buffers must exist and be typed either way. Built only from
            % Nontunable properties, so it is a compile-time constant.
            if isAlgebraic(obj)
                n_buf = 1;
            else
                n     = obj.window_samples;
                n     = n + mod(n, 2);            % even (Simpson)
                n_buf = n + 1;
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

            % 1. Riachy's auxiliary output: Y = y + Kd*int y
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
            obj.int_km1      = 0;
            obj.y_km1        = 0;
            obj.z_km1        = 0;
            obj.z_km2        = 0;
            obj.num_filt_km1 = 0;
            obj.num_filt_km2 = 0;
            obj.den_filt_km1 = 0;
            obj.den_filt_km2 = 0;
            % Full-size zeros() with a codegen-constant length: code
            % generation types the discrete states from these assignments.
            n_buf     = bufferLength(obj);
            obj.y_buf = zeros(n_buf, 1);
            obj.u_buf = zeros(n_buf, 1);
        end

        function [sz, dt, cp] = getDiscreteStateSpecificationImpl(obj, name)
            switch name
                case {'y_buf', 'u_buf'}
                    sz = [bufferLength(obj), 1];
                otherwise
                    sz = [1 1];
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
            [varargout{1:getNumOutputsImpl(obj)}] = deal([1 1]);
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
                n    = obj.window_samples + mod(obj.window_samples, 2);
                kind = sprintf('window, Tw = %g s', n*obj.Ts);
            end
            icon = sprintf('F-hat Riachy 2nd\nY = y + %g*int y\n%s', obj.Kd, kind);
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_fhat_riachy2_block', ...
                'Title', 'MFC F-hat: Riachy trick, 2nd order', ...
                'Text', sprintf(['Estimates Fk = F + Kd*dot_y from the auxiliary ', ...
                    'output Y = y + Kd*int y, for which ddot_Y = Fk + alpha*u.\n\n', ...
                    'The derivative term comes back inside F_hat, so the loop needs no ', ...
                    'derivative estimate at all. Wire the remaining feedback as a stock ', ...
                    'Discrete PID with P = Kp, I = Ki and D = 0 on err = y - sp_filt, ', ...
                    'and feed the command block ff = ddot_sp + Kd*dot_sp. Putting D on ', ...
                    'that PID applies Kd twice.\n\n', ...
                    'The estimator applied to Y is a parameter, not a separate block: ', ...
                    'algebraic (growing window, needs a run-start clock) or sliding ', ...
                    'window (Simpson, fixed memory, no internal smoothing). int y is ', ...
                    'unbounded -- with a non-zero steady state the sliding window is ', ...
                    'the better-behaved choice.']));
        end
    end
end
