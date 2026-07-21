classdef mfc_fhat_alg2_block < matlab.System
    % mfc_fhat_alg2_block  Algebraic (growing-window) F estimator, 2nd order.
    %
    %   Estimates F in the second-order ultra-local model
    %
    %       ddot_z = F + alpha*u + a_fold*dot_z + b_fold*z
    %
    %   by the algebraic / operational-calculus method: the model is written
    %   in the Laplace domain, differentiated twice w.r.t. s to annihilate
    %   both unknown initial conditions, and mapped back to the time domain
    %   where multiplication by t^n replaces (-d/ds)^n. The window grows
    %   from t = 0, so this estimator IS sensitive to the choice of time
    %   origin -- feed it a clock that starts with the run.
    %
    %   DRIVE SIGNAL AND FOLDING (the coupled/decoupled choice, made here by
    %   what you wire into z and what you set a_fold/b_fold to):
    %     decoupled: z = measurement,   a_fold = b_fold = 0
    %                F_hat is the pure-plant lumped dynamics ddot_y - alpha*u;
    %                stabilize with an explicit mfc_feedback_block.
    %     coupled:   z = tracking error, a_fold = -Kd, b_fold = -Kp
    %                the closed-loop polynomial s^2 + Kd*s + Kp is folded
    %                into the estimate, so no explicit P/D feedback is needed.
    %
    %   DISSECTION. By default the block smooths numerator and denominator
    %   internally (with the SAME filter) before dividing -- that is
    %   mathematically essential, not a convenience. Two mask options open
    %   it up:
    %     expose_raw      adds num_raw / den_raw outputs so you can watch or
    %                     re-filter the ratio's ingredients
    %     internal_filter clear it to set the internal smoother memory to 0
    %                     (an exact pass-through), turning this block into a
    %                     pure num/den generator. Rebuild the chain as
    %                       num_raw -> mfc_smoother_block -.
    %                       den_raw -> mfc_smoother_block -+-> mfc_fhat_divide_block
    %                     with the SAME window on both smoothers.
    %
    %   Ports
    %     In : z, u_prev, t  (+ optional alpha)
    %     Out: F_hat, valid  (+ optional num_raw, den_raw)
    %          u_prev is the command that actually reached the plant over the
    %          LAST sample -- in a loop assembled from separate blocks this
    %          must come through an explicit unit delay.
    %
    %   The math is mfc_fhat_algebraic_second_order; this class only maps
    %   parameters and Simulink state onto it.
    %
    %   See also mfc_fhat_algebraic_second_order, mfc_fhat_alg1_block,
    %   mfc_fhat_divide_block, mfc_smoother_block, mfc_siso_core.

    properties
        % alpha Ultra-local model input gain (ignored if the live alpha input is enabled)
        alpha = 1
        % a_fold Folded dot_z coefficient: -Kd when coupled, 0 when decoupled
        a_fold = 0
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
        function obj = mfc_fhat_alg2_block(varargin)
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

            [F_hat, state, dbg] = mfc_fhat_algebraic_second_order( ...
                z, u_prev, alpha_k, t, obj.Ts, win, obj.est_hold_time, ...
                obj.a_fold, obj.b_fold, state);

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
            if obj.a_fold == 0 && obj.b_fold == 0
                mode = 'decoupled';
            else
                mode = 'coupled (folded)';
            end
            icon = sprintf('F-hat algebraic\n2nd order\n%s', mode);
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_fhat_alg2_block', ...
                'Title', 'MFC F-hat: algebraic, 2nd order', ...
                'Text', sprintf(['Growing-window operational-calculus estimator for ', ...
                    'ddot_z = F + alpha*u + a_fold*dot_z + b_fold*z.\n\n', ...
                    'Decoupled: drive z with the measurement and leave a_fold = b_fold = 0. ', ...
                    'Coupled: drive z with the tracking error and set a_fold = -Kd, ', ...
                    'b_fold = -Kp to fold s^2 + Kd*s + Kp into the estimate.\n\n', ...
                    'The window grows from t = 0, so t must be a clock that starts with ', ...
                    'the run. Enable expose_raw (and optionally clear internal_filter) to ', ...
                    'rebuild the num/den smoothing and division from separate blocks.']));
        end
    end
end
