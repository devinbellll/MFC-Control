classdef mfc_fhat_alg1_coupled_mimo_block < matlab.System
    % mfc_fhat_alg1_coupled_mimo_block  Algebraic F estimator, 1st order,
    % COUPLED, with MATRIX gains.
    %
    %   The vector-valued twin of mfc_fhat_alg1_coupled_block: the
    %   first-order ultra-local model written about the TRACKING ERROR, with
    %   the desired closed-loop pole folded in,
    %
    %       dot_e = F + alpha*u - Kp*e,      e = y - y_sp
    %
    %   where e, u and F are n-by-1 and alpha, Kp are square n-by-n. The
    %   folded quantity is a MATRIX: off-diagonal Kp terms fold
    %   cross-channel proportional action into the estimate, which is what
    %   this block does that n SISO coupled blocks cannot.
    %
    %   NO DERIVATIVE ROOM AT FIRST ORDER. There is no Kd here, and that is
    %   structural, not an omission: the first-order model has no dot_e term
    %   to fold one into. Use mfc_fhat_alg2_coupled_mimo_block if the plant
    %   needs D action.
    %
    %   P action is already inside F_hat, so tie mfc_command_mimo_block's fb
    %   input to Ground, or to a stock discrete integrator for Ki; feed its
    %   ff input dot_sp (first order). Applying Kp again outside doubles it.
    %
    %   SIGN CONVENTION: e = y - y_sp (measurement minus filtered setpoint).
    %   Feeding y_sp - y inverts the folded pole and the loop diverges.
    %
    %   Ports
    %     In : err     tracking error e = y - y_sp   (n-by-1)
    %          u_prev  command that actually reached the plant over the LAST
    %                  sample, through an explicit unit delay  (n-by-1)
    %          t       clock, starting with the run (scalar)
    %          alpha   optional live gain            (n-by-n)
    %     Out: F_hat                                 (n-by-1)
    %
    %   The math is mfc_fhat_algebraic_first_order with b_fold = -Kp, a matrix.
    %
    %   See also mfc_fhat_alg1_coupled_block, mfc_fhat_alg2_coupled_mimo_block,
    %   mfc_fhat_alg1_decoupled_mimo_block, mfc_command_mimo_block.

    properties
        % alpha Ultra-local model input gain, square n-by-n (ignored if the live alpha input is enabled)
        alpha = eye(2)
        % Kp Proportional gain matrix, folded in as b_fold = -Kp (n-by-n). Do NOT also apply it externally.
        Kp = 25*eye(2)
        % est_filter_window Internal num/den smoother memory [samples]
        est_filter_window = 10
        % est_hold_time F_hat held at zero until t exceeds this [s] (guards the near-zero denominator at startup)
        est_hold_time = 0.1
    end

    properties (Nontunable)
        % n Number of channels (alpha and Kp are n-by-n, signals are n-by-1)
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
        function obj = mfc_fhat_alg1_coupled_mimo_block(varargin)
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

            % Coupled: error-driven, the matrix pole s + Kp folded in.
            [F_hat, state] = mfc_fhat_algebraic_first_order( ...
                err, u_prev, alpha_k, t, obj.Ts, obj.est_filter_window, ...
                obj.est_hold_time, -obj.Kp, state);

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
            obj.den_filt_km1 = 0;   % denominator is t: scalar, shared
            obj.den_filt_km2 = 0;
        end

        function validatePropertiesImpl(obj)
            if ~obj.use_live_alpha && ~isequal(size(obj.alpha), [obj.n obj.n])
                error('mfc:mimo:alphaSize', ...
                    'alpha must be %d-by-%d to match n.', obj.n, obj.n);
            end
            if ~isequal(size(obj.Kp), [obj.n obj.n])
                error('mfc:mimo:gainSize', ...
                    'Kp must be %d-by-%d to match n.', obj.n, obj.n);
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
            icon = sprintf('F-hat algebraic\n1st order, coupled\nmatrix Kp (%dx%d), no Kd', ...
                obj.n, obj.n);
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_fhat_alg1_coupled_mimo_block', ...
                'Title', 'MFC F-hat: algebraic 1st order, coupled, matrix gains', ...
                'Text', sprintf(['Error-driven growing-window estimator for ', ...
                    'dot_e = F + alpha*u - Kp*e with n-by-1 signals and square n-by-n ', ...
                    'alpha and Kp.\n\n', ...
                    'The matrix pole s + Kp is folded into F_hat, so tie ', ...
                    'mfc_command_mimo_block''s fb input to Ground (or a discrete ', ...
                    'integrator for Ki) and feed its ff input dot_sp. Applying Kp again ', ...
                    'outside doubles it.\n\nThere is no Kd at first order: the model has ', ...
                    'no derivative term to fold one into.']));
        end
    end
end
