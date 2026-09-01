classdef mfc_fhat_alg2_coupled_mimo_block < matlab.System
    % mfc_fhat_alg2_coupled_mimo_block  Algebraic F estimator, 2nd order,
    % COUPLED, with MATRIX gains.
    %
    %   The vector-valued twin of mfc_fhat_alg2_coupled_block: the
    %   second-order ultra-local model written about the TRACKING ERROR,
    %   with the desired closed-loop polynomial folded in,
    %
    %       ddot_e = F + alpha*u - Kd*dot_e - Kp*e,      e = y - y_sp
    %
    %   where e, u and F are n-by-1 and alpha, Kp, Kd are square n-by-n. The
    %   folded quantity is a MATRIX polynomial: the closed loop is
    %   ddot_e + Kd*dot_e + Kp*e = 0, so off-diagonal Kp/Kd terms fold
    %   cross-channel proportional and derivative action into the estimate.
    %   That is the one thing this block does that n SISO coupled blocks
    %   side by side cannot.
    %
    %   WHAT THAT BUYS YOU, unchanged from the SISO case: P and D action are
    %   already inside F_hat, so tie mfc_command_mimo_block's fb input to
    %   Ground, or to a stock discrete integrator for integral action (Ki is
    %   the only gain still applied outside). Applying Kp or Kd again in an
    %   external controller doubles them.
    %
    %   The price is also unchanged: F_hat is a closed-loop quantity, not
    %   the plant's own dynamics. Use mfc_fhat_alg2_decoupled_mimo_block
    %   when you want to look at the plant.
    %
    %   SIGN CONVENTION: e = y - y_sp (measurement minus filtered setpoint).
    %   Feeding y_sp - y inverts the folded poles and the loop diverges.
    %
    %   Ports
    %     In : err     tracking error e = y - y_sp   (n-by-1)
    %          u_prev  command that actually reached the plant over the LAST
    %                  sample, through an explicit unit delay  (n-by-1)
    %          t       clock, starting with the run (scalar)
    %          alpha   optional live gain            (n-by-n)
    %     Out: F_hat                                 (n-by-1)
    %
    %   The math is mfc_fhat_algebraic_second_order with a_fold = -Kd and
    %   b_fold = -Kp, both matrices.
    %
    %   See also mfc_fhat_alg2_coupled_block, mfc_fhat_alg1_coupled_mimo_block,
    %   mfc_fhat_alg2_decoupled_mimo_block, mfc_command_mimo_block.

    properties
        % alpha Ultra-local model input gain, square n-by-n (ignored if the live alpha input is enabled)
        alpha = eye(2)
        % Kp Proportional gain matrix, folded in as b_fold = -Kp (n-by-n). Do NOT also apply it externally.
        Kp = 25*eye(2)
        % Kd Derivative gain matrix, folded in as a_fold = -Kd (n-by-n). Do NOT also apply it externally.
        Kd = 10*eye(2)
        % est_filter_window Internal num/den smoother memory [samples]
        est_filter_window = 10
        % est_hold_time F_hat held at zero until t exceeds this [s] (guards the near-zero denominator at startup)
        est_hold_time = 0.1
    end

    properties (Nontunable)
        % n Number of channels (alpha, Kp and Kd are n-by-n, signals are n-by-1)
        n = 2
        % Ts Sample time [s] (fixes the block's discrete rate)
        Ts = 0.01
    end

    properties (Nontunable, Logical)
        % use_live_alpha Add the alpha input port (overrides the alpha parameter)
        use_live_alpha = false
    end

    properties (DiscreteState)
        z_km1
        z_km2
        num_filt_km1
        num_filt_km2
        den_filt_km1
        den_filt_km2
    end

    methods
        function obj = mfc_fhat_alg2_coupled_mimo_block(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)

        function F_hat = stepImpl(obj, err, u_prev, t, varargin)
            if obj.use_live_alpha
                alpha_k = varargin{1};
            else
                alpha_k = obj.alpha;
            end

            state = struct( ...
                'z_km1',        obj.z_km1, ...
                'z_km2',        obj.z_km2, ...
                'num_filt_km1', obj.num_filt_km1, ...
                'num_filt_km2', obj.num_filt_km2, ...
                'den_filt_km1', obj.den_filt_km1, ...
                'den_filt_km2', obj.den_filt_km2);

            % Coupled: error-driven, matrix polynomial s^2 + Kd*s + Kp folded in.
            [F_hat, state] = mfc_fhat_algebraic_second_order( ...
                err, u_prev, alpha_k, t, obj.Ts, obj.est_filter_window, ...
                obj.est_hold_time, -obj.Kd, -obj.Kp, state);

            obj.z_km1        = state.z_km1;
            obj.z_km2        = state.z_km2;
            obj.num_filt_km1 = state.num_filt_km1;
            obj.num_filt_km2 = state.num_filt_km2;
            obj.den_filt_km1 = state.den_filt_km1;
            obj.den_filt_km2 = state.den_filt_km2;
        end

        function resetImpl(obj)
            obj.z_km1        = zeros(obj.n, 1);
            obj.z_km2        = zeros(obj.n, 1);
            obj.num_filt_km1 = zeros(obj.n, 1);
            obj.num_filt_km2 = zeros(obj.n, 1);
            obj.den_filt_km1 = 0;   % denominator is t^2: scalar, shared
            obj.den_filt_km2 = 0;
        end

        function validatePropertiesImpl(obj)
            if ~obj.use_live_alpha && ~isequal(size(obj.alpha), [obj.n obj.n])
                error('mfc:mimo:alphaSize', ...
                    'alpha must be %d-by-%d to match n.', obj.n, obj.n);
            end
            if ~isequal(size(obj.Kp), [obj.n obj.n]) || ~isequal(size(obj.Kd), [obj.n obj.n])
                error('mfc:mimo:gainSize', ...
                    'Kp and Kd must be %d-by-%d to match n.', obj.n, obj.n);
            end
        end

        function [sz, dt, cp] = getDiscreteStateSpecificationImpl(obj, name)
            switch name
                case {'den_filt_km1', 'den_filt_km2'}
                    sz = [1 1];
                otherwise
                    sz = [obj.n 1];
            end
            dt = 'double';  cp = false;
        end

        function sts = getSampleTimeImpl(obj)
            sts = createSampleTime(obj, 'Type', 'Discrete', 'SampleTime', obj.Ts);
        end

        % ---- Ports -------------------------------------------------------
        function num = getNumInputsImpl(obj), num = 3 + obj.use_live_alpha; end
        function varargout = getInputNamesImpl(obj)
            names = {'err', 'u_prev', 't'};
            if obj.use_live_alpha, names{end+1} = 'alpha'; end
            varargout = names;
        end

        function num = getNumOutputsImpl(~), num = 1; end
        function varargout = getOutputNamesImpl(~),    varargout = {'F_hat'}; end
        function varargout = getOutputSizeImpl(obj),   varargout = {[obj.n 1]}; end
        function varargout = getOutputDataTypeImpl(~), varargout = {'double'}; end
        function varargout = isOutputComplexImpl(~),   varargout = {false}; end
        function varargout = isOutputFixedSizeImpl(~), varargout = {true}; end

        function icon = getIconImpl(obj)
            icon = sprintf('F-hat algebraic\n2nd order\ncoupled (err-driven)\nmatrix Kp, Kd (%dx%d)', ...
                obj.n, obj.n);
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_fhat_alg2_coupled_mimo_block', ...
                'Title', 'MFC F-hat: algebraic 2nd order, coupled, matrix gains', ...
                'Text', sprintf(['Error-driven growing-window estimator for ', ...
                    'ddot_e = F + alpha*u - Kd*dot_e - Kp*e with n-by-1 signals and ', ...
                    'square n-by-n alpha, Kp, Kd.\n\n', ...
                    'The matrix polynomial ddot_e + Kd*dot_e + Kp*e is folded into ', ...
                    'F_hat, so no explicit P or D feedback is needed: tie ', ...
                    'mfc_command_mimo_block''s fb input to Ground, or to a discrete ', ...
                    'integrator for Ki. Applying Kp or Kd again outside doubles them.\n\n', ...
                    'err must be y - y_sp, measurement minus FILTERED setpoint; the ', ...
                    'other sign inverts the folded poles.']));
        end
    end
end
