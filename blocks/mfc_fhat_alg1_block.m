classdef mfc_fhat_alg1_block < matlab.System
    % mfc_fhat_alg1_block  Algebraic (growing-window) F estimator, 1st order.
    %
    %   Estimates F in the first-order ultra-local model
    %
    %       dot_z = F + alpha*u + b_fold*z
    %
    %   by the algebraic / operational-calculus method: the model is written
    %   in the Laplace domain, differentiated ONCE w.r.t. s to annihilate
    %   the single unknown initial condition, and mapped back to the time
    %   domain where multiplication by t replaces -d/ds. The window grows
    %   from t = 0, so this estimator IS sensitive to the choice of time
    %   origin -- feed it a clock that starts with the run.
    %
    %   DRIVE SIGNAL AND FOLDING:
    %     decoupled: z = measurement,    b_fold = 0
    %                F_hat is the pure-plant lumped dynamics dot_y - alpha*u;
    %                stabilize with an explicit mfc_feedback_block.
    %     coupled:   z = tracking error, b_fold = -Kp
    %                the closed-loop pole s + Kp is folded into the estimate.
    %
    %   NO DERIVATIVE ROOM AT FIRST ORDER. There is no a_fold here: the
    %   first-order model has no dot_z term to fold a derivative gain into.
    %   A coupled first-order loop therefore cannot supply damping at any
    %   gain, which is why it diverges on plants that need it (a double
    %   integrator, say). If you need Kd at first order, run DECOUPLED and
    %   let mfc_feedback_block apply it explicitly.
    %
    %   Ports
    %     In : z, u_prev, t  (+ optional alpha)
    %     Out: F_hat, valid  (+ optional num_raw, den_raw)
    %
    %   The math is mfc_fhat_algebraic_first_order; this class only maps
    %   parameters and Simulink state onto it.
    %
    %   See also mfc_fhat_algebraic_first_order, mfc_fhat_alg2_block,
    %   mfc_fhat_divide_block, mfc_smoother_block, mfc_siso_core.

    properties
        % alpha Ultra-local model input gain (ignored if the live alpha input is enabled)
        alpha = 1
        % b_fold Folded z coefficient: -Kp when coupled, 0 when decoupled
        b_fold = 0
    end

    properties (Nontunable)
        % Ts Sample time [s] (fixes the block's discrete rate)
        Ts = 0.01
        % est_filter_window Internal num/den smoother memory [samples]
        est_filter_window = 10
        % est_hold_time F_hat held at zero until t exceeds this [s] (guards the near-zero denominator at startup)
        est_hold_time = 0.1
    end

    properties (Nontunable, Logical)
        % internal_filter Smooth numerator and denominator inside the block (clear to emit them unsmoothed)
        internal_filter = true
        % expose_raw Add the num_raw and den_raw output ports
        expose_raw = false
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
        function obj = mfc_fhat_alg1_block(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)

        function varargout = stepImpl(obj, z, u_prev, t, varargin)
            if obj.use_live_alpha
                alpha_k = varargin{1};
            else
                alpha_k = obj.alpha;
            end

            % window = 0 is an exact pass-through in mfc_iir_smoother, so
            % "no internal filter" needs no separate code path.
            if obj.internal_filter
                win = obj.est_filter_window;
            else
                win = 0;
            end

            state = struct( ...
                'z_km1',        obj.z_km1, ...
                'z_km2',        obj.z_km2, ...
                'num_filt_km1', obj.num_filt_km1, ...
                'num_filt_km2', obj.num_filt_km2, ...
                'den_filt_km1', obj.den_filt_km1, ...
                'den_filt_km2', obj.den_filt_km2);

            [F_hat, state, dbg] = mfc_fhat_algebraic_first_order( ...
                z, u_prev, alpha_k, t, obj.Ts, win, obj.est_hold_time, ...
                obj.b_fold, state);

            obj.z_km1        = state.z_km1;
            obj.z_km2        = state.z_km2;
            obj.num_filt_km1 = state.num_filt_km1;
            obj.num_filt_km2 = state.num_filt_km2;
            obj.den_filt_km1 = state.den_filt_km1;
            obj.den_filt_km2 = state.den_filt_km2;

            varargout{1} = F_hat;
            varargout{2} = double(dbg.valid);
            if obj.expose_raw
                varargout{3} = dbg.num_raw;
                varargout{4} = dbg.den_raw;
            end
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
            names = {'z', 'u_prev', 't'};
            if obj.use_live_alpha, names{end+1} = 'alpha'; end
            varargout = names;
        end

        function num = getNumOutputsImpl(obj), num = 2 + 2*obj.expose_raw; end
        function varargout = getOutputNamesImpl(obj)
            if obj.expose_raw
                varargout = {'F_hat', 'valid', 'num_raw', 'den_raw'};
            else
                varargout = {'F_hat', 'valid'};
            end
        end
        function varargout = getOutputSizeImpl(obj)
            varargout = repmat({[1 1]}, 1, getNumOutputsImpl(obj));
        end
        function varargout = getOutputDataTypeImpl(obj)
            varargout = repmat({'double'}, 1, getNumOutputsImpl(obj));
        end
        function varargout = isOutputComplexImpl(obj)
            varargout = repmat({false}, 1, getNumOutputsImpl(obj));
        end
        function varargout = isOutputFixedSizeImpl(obj)
            varargout = repmat({true}, 1, getNumOutputsImpl(obj));
        end

        function icon = getIconImpl(obj)
            if obj.b_fold == 0
                mode = 'decoupled';
            else
                mode = 'coupled (folded)';
            end
            icon = sprintf('F-hat algebraic\n1st order\n%s', mode);
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_fhat_alg1_block', ...
                'Title', 'MFC F-hat: algebraic, 1st order', ...
                'Text', sprintf(['Growing-window operational-calculus estimator for ', ...
                    'dot_z = F + alpha*u + b_fold*z.\n\n', ...
                    'Decoupled: drive z with the measurement and leave b_fold = 0. ', ...
                    'Coupled: drive z with the tracking error and set b_fold = -Kp to ', ...
                    'fold the pole s + Kp into the estimate.\n\n', ...
                    'There is no a_fold at first order -- the model has no dot_z term to ', ...
                    'fold a derivative gain into, so a coupled first-order loop cannot ', ...
                    'supply damping. Use the decoupled structure if you need Kd.']));
        end
    end
end
