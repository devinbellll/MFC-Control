classdef mfc_fhat_alg2_decoupled_block < matlab.System
    % mfc_fhat_alg2_decoupled_block  Algebraic F estimator, 2nd order, DECOUPLED.
    %
    %   Estimates F in the second-order ultra-local model
    %
    %       ddot_y = F + alpha*u
    %
    %   from the plant MEASUREMENT, by the algebraic / operational-calculus
    %   method: the model is written in the Laplace domain, differentiated
    %   twice w.r.t. s to annihilate both unknown initial conditions, and
    %   mapped back to the time domain where multiplication by t^n replaces
    %   (-d/ds)^n. The window grows from t = 0, so this estimator IS
    %   sensitive to the choice of time origin -- feed it a clock that
    %   starts with the run.
    %
    %   DECOUPLED means F_hat is the TRUE PLANT lumped dynamics
    %   ddot_y - alpha*u: no closed-loop dynamics are folded into it, and it
    %   is not stabilizing on its own. You must add an explicit feedback law
    %   -- a stock Discrete PID Controller on the tracking error -- into the
    %   fb input of mfc_command_block. See mfc_fhat_alg2_coupled_block for
    %   the alternative, where the feedback poles live inside the estimate
    %   instead.
    %
    %   Because the estimate is plant-only, it stays meaningful when the
    %   loop is opened, the setpoint jumps, or the feedback gains are
    %   retuned -- which is what makes this the variant worth logging.
    %
    %   Ports
    %     In : y       plant measurement
    %          u_prev  command that actually reached the plant over the LAST
    %                  sample; in a loop assembled from separate blocks this
    %                  must come through an explicit unit delay
    %          t       clock, starting with the run
    %          alpha   optional live gain (overrides the mask parameter)
    %     Out: F_hat
    %
    %   The math is mfc_fhat_algebraic_second_order with a_fold = b_fold = 0;
    %   this class only maps parameters and Simulink state onto it.
    %
    %   See also mfc_fhat_algebraic_second_order, mfc_fhat_alg2_coupled_block,
    %   mfc_fhat_alg1_decoupled_block, mfc_command_block, mfc_siso_core.

    properties
        % alpha Ultra-local model input gain (ignored if the live alpha input is enabled)
        alpha = 1
        % est_filter_window Internal num/den smoother memory [samples]
        est_filter_window = 10
        % est_hold_time F_hat held at zero until t exceeds this [s] (guards the near-zero denominator at startup)
        est_hold_time = 0.1
    end

    properties (Nontunable)
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
        function obj = mfc_fhat_alg2_decoupled_block(varargin)
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
            obj.z_km1        = 0;
            obj.z_km2        = 0;
            obj.num_filt_km1 = 0;
            obj.num_filt_km2 = 0;
            obj.den_filt_km1 = 0;
            obj.den_filt_km2 = 0;
        end

        function [sz, dt, cp] = getDiscreteStateSpecificationImpl(~, ~)
            sz = [1 1];  dt = 'double';  cp = false;
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
        function varargout = getOutputSizeImpl(~),     varargout = {[1 1]}; end
        function varargout = getOutputDataTypeImpl(~), varargout = {'double'}; end
        function varargout = isOutputComplexImpl(~),   varargout = {false}; end
        function varargout = isOutputFixedSizeImpl(~), varargout = {true}; end

        function icon = getIconImpl(~)
            icon = sprintf('F-hat algebraic\n2nd order\ndecoupled (y-driven)');
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_fhat_alg2_decoupled_block', ...
                'Title', 'MFC F-hat: algebraic 2nd order, decoupled', ...
                'Text', sprintf(['Growing-window operational-calculus estimator for ', ...
                    'ddot_y = F + alpha*u, driven by the plant measurement.\n\n', ...
                    'F_hat is the true-plant lumped dynamics -- nothing is folded in, so ', ...
                    'it does not stabilize by itself. Add an explicit feedback law (a ', ...
                    'stock Discrete PID on the tracking error) into the fb input of the ', ...
                    'command block.\n\nThe window grows from t = 0, so t must be a clock ', ...
                    'that starts with the run.']));
        end
    end
end
