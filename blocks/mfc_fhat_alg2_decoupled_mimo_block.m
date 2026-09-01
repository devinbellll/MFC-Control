classdef mfc_fhat_alg2_decoupled_mimo_block < matlab.System
    % mfc_fhat_alg2_decoupled_mimo_block  Algebraic F estimator, 2nd order,
    % DECOUPLED, with a MATRIX input gain.
    %
    %   The vector-valued twin of mfc_fhat_alg2_decoupled_block: estimates F
    %   in the second-order ultra-local model
    %
    %       ddot_y = F + alpha*u
    %
    %   where y, u and F are n-by-1 and alpha is a square n-by-n matrix, so
    %   the channels are coupled through the input gain only. Everything
    %   else -- the operational-calculus window, the shared num/den
    %   smoother, the hold time -- is identical and acts channel-wise; the
    %   denominator t^2 is scalar and shared, which is what keeps one window
    %   valid for the whole vector.
    %
    %   DECOUPLED still means what it means in the SISO case: F_hat is the
    %   true-plant lumped dynamics ddot_y - alpha*u, nothing folded, not
    %   stabilizing on its own. Add explicit feedback into the fb input of
    %   mfc_command_mimo_block.
    %
    %   Ports
    %     In : y       plant measurement            (n-by-1)
    %          u_prev  command that actually reached the plant over the LAST
    %                  sample, through an explicit unit delay  (n-by-1)
    %          t       clock, starting with the run (scalar)
    %          alpha   optional live gain            (n-by-n)
    %     Out: F_hat                                 (n-by-1)
    %
    %   The math is mfc_fhat_algebraic_second_order with a_fold = b_fold = 0,
    %   unchanged -- it is already vector-safe.
    %
    %   See also mfc_fhat_alg2_decoupled_block, mfc_command_mimo_block,
    %   mfc_fhat_algebraic_second_order.

    properties
        % alpha Ultra-local model input gain, square n-by-n (ignored if the live alpha input is enabled)
        alpha = eye(2)
        % est_filter_window Internal num/den smoother memory [samples]
        est_filter_window = 10
        % est_hold_time F_hat held at zero until t exceeds this [s] (guards the near-zero denominator at startup)
        est_hold_time = 0.1
    end

    properties (Nontunable)
        % n Number of channels (alpha is n-by-n, signals are n-by-1)
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
        function obj = mfc_fhat_alg2_decoupled_mimo_block(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)

        function F_hat = stepImpl(obj, y, u_prev, t, varargin)
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

            % Decoupled: driven by the measurement, nothing folded.
            [F_hat, state] = mfc_fhat_algebraic_second_order( ...
                y, u_prev, alpha_k, t, obj.Ts, obj.est_filter_window, ...
                obj.est_hold_time, 0, 0, state);

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
            % Fixed discrete rate (do NOT inherit): the backward-difference
            % recursion assumes it advances exactly once per Ts.
            sts = createSampleTime(obj, 'Type', 'Discrete', 'SampleTime', obj.Ts);
        end

        % ---- Ports -------------------------------------------------------
        function num = getNumInputsImpl(obj), num = 3 + obj.use_live_alpha; end
        function varargout = getInputNamesImpl(obj)
            names = {'y', 'u_prev', 't'};
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
            icon = sprintf('F-hat algebraic\n2nd order\ndecoupled (y-driven)\nmatrix alpha (%dx%d)', ...
                obj.n, obj.n);
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_fhat_alg2_decoupled_mimo_block', ...
                'Title', 'MFC F-hat: algebraic 2nd order, decoupled, matrix alpha', ...
                'Text', sprintf(['Growing-window operational-calculus estimator for ', ...
                    'ddot_y = F + alpha*u with vector y, u, F and a square n-by-n ', ...
                    'alpha, driven by the plant measurement.\n\n', ...
                    'F_hat is the true-plant lumped dynamics -- nothing is folded in, so ', ...
                    'it does not stabilize by itself. Add an explicit feedback law into ', ...
                    'the fb input of mfc_command_mimo_block.\n\nThe window grows from ', ...
                    't = 0, so t must be a clock that starts with the run.']));
        end
    end
end
