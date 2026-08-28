classdef mfc_fhat_window_mimo_block < matlab.System
    % mfc_fhat_window_mimo_block  Sliding-window (FIR) F estimator with a
    % MATRIX input gain.
    %
    %   The vector-valued twin of mfc_fhat_window_block: estimates F in
    %
    %       dot_y  = F + alpha*u        (first order)
    %       ddot_y = F + alpha*u        (second order)
    %
    %   with y, u and F n-by-1 and alpha a square n-by-n matrix, from a
    %   fixed-length weighted integral of the last Tw = window_samples*Ts
    %   seconds of the measurement and applied-input histories. The taps are
    %   scalar quadrature weights and are SHARED across channels -- one
    %   window serves the whole vector, exactly as the algebraic estimators'
    %   scalar t/t^2 denominator does. alpha multiplies the input-window sum
    %   as a matrix, so it mixes channels the same way mfc_siso.command_mimo
    %   does when the loop inverts it again.
    %
    %   DECOUPLED ONLY, and structurally so: an FIR window has nowhere to
    %   fold closed-loop poles, so there is no coupled counterpart at any
    %   width. F_hat is the true-plant lumped dynamics; add explicit
    %   feedback into mfc_command_mimo_block's fb input.
    %
    %   No internal smoothing and no poles: a sample enters, is weighted for
    %   exactly Tw seconds, and leaves completely. If F_hat is noisy, follow
    %   this block with mfc_smoother_mimo_block -- unlike the algebraic
    %   estimators, this one does not filter anything for you.
    %
    %   The estimate is held at zero (n-by-1 zeros) until the window fills.
    %
    %   Ports
    %     In : y       plant measurement            (n-by-1)
    %          u_prev  command that actually reached the plant over the LAST
    %                  sample, through an explicit unit delay  (n-by-1)
    %          t       clock; only used for the window-fill hold (scalar)
    %          alpha   optional live gain            (n-by-n)
    %     Out: F_hat                                 (n-by-1)
    %
    %   The math is mfc_siso.window_kernel + mfc_fhat_sliding_window.
    %
    %   See also mfc_fhat_window_block, mfc_fhat_alg2_decoupled_mimo_block,
    %   mfc_smoother_mimo_block, mfc_command_mimo_block, mfc_mimo_core.

    properties
        % alpha Ultra-local model input gain, square n-by-n (ignored if the live alpha input is enabled)
        alpha = eye(2)
    end

    properties (Nontunable)
        % n Number of channels (alpha is n-by-n, signals are n-by-1)
        n = 2
        % model_order Order of the ultra-local model
        model_order = 'Second order (ddot_y = F + alpha*u)'
        % Ts Sample time [s] (fixes the block's discrete rate)
        Ts = 0.01
        % window_samples Window length [intervals]; realized window Tw = window_samples*Ts. Nontunable: sizes the buffers and the taps.
        window_samples = 10
    end

    properties (Hidden, Constant)
        model_orderSet = matlab.system.StringSet({ ...
            'First order (dot_y = F + alpha*u)', ...
            'Second order (ddot_y = F + alpha*u)'});
    end

    properties (Nontunable, Logical)
        % use_live_alpha Add the alpha input port (overrides the alpha parameter)
        use_live_alpha = false
    end

    properties (DiscreteState)
        y_buf
        u_buf
    end

    properties (Access = private)
        kernel     % quadrature kernel built by mfc_siso.window_kernel in setupImpl
    end

    methods
        function obj = mfc_fhat_window_mimo_block(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = private)
        function ord = orderNum(obj)
            if strncmp(obj.model_order, 'First', 5), ord = 1; else, ord = 2; end
        end

        function sfx = ordinalSuffix(~, ord)
            if ord == 1, sfx = 'st'; else, sfx = 'nd'; end
        end
    end

    methods (Access = protected)

        function setupImpl(obj)
            obj.kernel = mfc_siso.window_kernel(orderNum(obj), obj.window_samples, obj.Ts);
        end

        function F_hat = stepImpl(obj, y, u_prev, t, varargin)
            if obj.use_live_alpha
                alpha_k = varargin{1};
            else
                alpha_k = obj.alpha;
            end

            state = struct('y_buf', obj.y_buf, 'u_buf', obj.u_buf);
            [F_hat, state] = mfc_fhat_sliding_window( ...
                y, u_prev, alpha_k, t, obj.kernel, state);
            obj.y_buf = state.y_buf;
            obj.u_buf = state.u_buf;
        end

        function resetImpl(obj)
            % One COLUMN per channel, with a codegen-constant length.
            obj.y_buf = zeros(obj.window_samples + 1, obj.n);
            obj.u_buf = zeros(obj.window_samples + 1, obj.n);
        end

        function validatePropertiesImpl(obj)
            if ~obj.use_live_alpha && ~isequal(size(obj.alpha), [obj.n obj.n])
                error('mfc:mimo:alphaSize', ...
                    'alpha must be %d-by-%d to match n.', obj.n, obj.n);
            end
        end

        function [sz, dt, cp] = getDiscreteStateSpecificationImpl(obj, ~)
            sz = [obj.window_samples + 1, obj.n];
            dt = 'double';  cp = false;
        end

        function sts = getSampleTimeImpl(obj)
            % Fixed discrete rate (do NOT inherit): the window buffers
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

        function num = getNumOutputsImpl(~), num = 1; end
        function varargout = getOutputNamesImpl(~),    varargout = {'F_hat'}; end
        function varargout = getOutputSizeImpl(obj),   varargout = {[obj.n 1]}; end
        function varargout = getOutputDataTypeImpl(~), varargout = {'double'}; end
        function varargout = isOutputComplexImpl(~),   varargout = {false}; end
        function varargout = isOutputFixedSizeImpl(~), varargout = {true}; end

        function icon = getIconImpl(obj)
            icon = sprintf('F-hat sliding window\n%d%s order, Tw = %g s\nmatrix alpha (%dx%d)', ...
                orderNum(obj), ordinalSuffix(obj, orderNum(obj)), ...
                obj.window_samples*obj.Ts, obj.n, obj.n);
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_fhat_window_mimo_block', ...
                'Title', 'MFC F-hat: sliding window (FIR), matrix alpha', ...
                'Text', sprintf(['Fixed-length FIR window estimator for the ultra-local ', ...
                    'model with n-by-1 signals and a square n-by-n alpha. The taps are ', ...
                    'scalar and shared across channels; alpha multiplies the input-window ', ...
                    'sum as a matrix.\n\n', ...
                    'Finite memory, no poles and no time origin to get wrong -- but no ', ...
                    'internal smoothing either: follow it with mfc_smoother_mimo_block if ', ...
                    'F_hat is noisy. Decoupled only (an FIR window has nowhere to fold ', ...
                    'poles): add explicit feedback into mfc_command_mimo_block''s fb ', ...
                    'input. F_hat is held at zero until the window fills.']));
        end
    end
end
