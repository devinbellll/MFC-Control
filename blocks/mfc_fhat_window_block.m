classdef mfc_fhat_window_block < matlab.System
    % mfc_fhat_window_block  Sliding-window (FIR) F estimator.
    %
    %   Estimates F in the ultra-local model
    %
    %       dot_y  = F + alpha*u      (model_order = 1, Eq. 11)
    %       ddot_y = F + alpha*u      (model_order = 2, Eq. 16)
    %
    %   by evaluating a FIXED-LENGTH weighted integral of the measurement
    %   and applied-input histories over the last Tw seconds:
    %
    %   The integral is precomputed into one fixed multiplier -- a TAP --
    %   per stored sample, so at run time this is a pure FIR filter:
    %
    %       F_hat = tap_y' * y_buf  +  alpha * ( tap_u_unit' * u_buf )
    %
    %   DECOUPLED BY CONSTRUCTION. This estimator sees only (y, u, alpha) --
    %   never the tracking error -- so F_hat is always the TRUE plant lumped
    %   dynamics and there is nothing to fold. There is no coupled variant
    %   of this block, and there cannot be one: stabilize with an explicit
    %   feedback law (a stock Discrete PID on the tracking error) into the
    %   fb input of mfc_command_block.
    %
    %   Unlike the algebraic estimators it has finite memory and no growing
    %   time weights, so it is INSENSITIVE to the choice of time origin --
    %   t is used only to hold the output at zero until the window fills.
    %   That makes it the right choice for a long or restarting run.
    %
    %   It also has no internal smoothing at all. If the estimate is noisy,
    %   follow it with an mfc_smoother_block on F_hat.
    %
    %   TAPS, not a quadrature rule. The weighting kernels are polynomials
    %   we wrote down ourselves, so each tap is the EXACT integral of the
    %   kernel against that sample's interpolation basis -- piecewise linear
    %   for y, piecewise constant for u (which is not an assumption at all:
    %   the command really does reach the plant through a zero-order hold,
    %   so the input term carries no quadrature error at all).
    %   window_samples means exactly what it says; nothing is rounded. See
    %   mfc_siso.window_kernel.
    %
    %   Ports
    %     In : y       plant measurement
    %          u_prev  command that actually reached the plant over the LAST
    %                  sample; in a loop assembled from separate blocks this
    %                  must come through an explicit unit delay
    %          t       clock (only used to hold the output until the window
    %                  fills -- unlike the algebraic blocks, the estimate
    %                  itself does not depend on the time origin)
    %          alpha   optional live gain (overrides the mask parameter)
    %     Out: F_hat
    %
    %   The math is mfc_fhat_sliding_window with a kernel precomputed once
    %   by mfc_siso.window_kernel in setupImpl; this class only maps
    %   parameters and Simulink state onto them.
    %
    %   See also mfc_fhat_sliding_window, mfc_siso,
    %   mfc_fhat_alg2_decoupled_block, mfc_smoother_block, mfc_siso_core.

    properties
        % alpha Ultra-local model input gain (ignored if the live alpha input is enabled)
        alpha = 1
    end

    properties (Nontunable)
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
        function obj = mfc_fhat_window_block(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = private)
        function n = orderNum(obj)
            if strncmp(obj.model_order, 'First', 5), n = 1; else, n = 2; end
        end

        function n_buf = bufferLength(obj)
            % Shared by resetImpl and getDiscreteStateSpecificationImpl so
            % the sizes cannot diverge. Built only from Nontunable
            % properties, so it is a compile-time constant under codegen.
            n_buf = obj.window_samples + 1;
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
            % Full-size zeros() with a codegen-constant length: code
            % generation types the discrete states from these assignments.
            n_buf = bufferLength(obj);
            obj.y_buf = zeros(n_buf, 1);
            obj.u_buf = zeros(n_buf, 1);
        end

        function [sz, dt, cp] = getDiscreteStateSpecificationImpl(obj, ~)
            sz = [bufferLength(obj), 1];  dt = 'double';  cp = false;
        end

        function sts = getSampleTimeImpl(obj)
            % Fixed discrete rate (do NOT inherit): the window buffers
            % assume they shift exactly once per Ts.
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
        function varargout = getOutputSizeImpl(~),     varargout = {[1 1]}; end
        function varargout = getOutputDataTypeImpl(~), varargout = {'double'}; end
        function varargout = isOutputComplexImpl(~),   varargout = {false}; end
        function varargout = isOutputFixedSizeImpl(~), varargout = {true}; end

        function icon = getIconImpl(obj)
            icon = sprintf('F-hat window\n%d order, FIR\nTw = %g s', ...
                           orderNum(obj), obj.window_samples*obj.Ts);
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_fhat_window_block', ...
                'Title', 'MFC F-hat: sliding window', ...
                'Text', sprintf(['Fixed-length FIR estimator for ', ...
                    'dot_y = F + alpha*u (1st order) or ddot_y = F + alpha*u (2nd).\n\n', ...
                    'Decoupled by construction: it never sees the tracking error, so ', ...
                    'F_hat is the true plant lumped dynamics and an explicit feedback ', ...
                    'law (a stock Discrete PID into the command block''s fb input) is ', ...
                    'required. Insensitive to the time origin (unlike the ', ...
                    'algebraic estimators) -- t only holds the output until the window ', ...
                    'fills.\n\nNo internal smoothing: add an mfc_smoother_block on F_hat ', ...
                    'if the estimate is noisy.']));
        end
    end
end
