classdef mfc_fhat_alg1_coupled_block < matlab.System
    % mfc_fhat_alg1_coupled_block  Algebraic F estimator, 1st order, COUPLED.
    %
    %   Estimates F in the first-order ultra-local model written about the
    %   TRACKING ERROR, with the desired closed-loop pole folded in:
    %
    %       dot_e = F + alpha*u - Kp*e,      e = y - y_sp
    %
    %   Same algebraic / operational-calculus machinery as the decoupled
    %   block -- Laplace domain, differentiated once w.r.t. s to annihilate
    %   the single unknown initial condition, mapped back with t for -d/ds
    %   -- but driven by the error and with the pole s + Kp absorbed into
    %   the estimate.
    %
    %   WHAT THAT BUYS YOU. The proportional action is already inside
    %   F_hat, so no explicit P feedback is needed: tie the command block's
    %   fb input to Ground, or to a stock discrete integrator if you want
    %   integral action (Ki is the ONLY gain that still has to be applied
    %   outside). Applying Kp again in an external PID doubles it.
    %
    %   NO DERIVATIVE ROOM AT FIRST ORDER. There is no Kd here, and that is
    %   structural, not an omission: the first-order model has no dot_e term
    %   to fold a derivative gain into, so this block cannot supply damping
    %   at any gain. It therefore diverges on plants that need damping (a
    %   double integrator, say). If you need D action at first order, use
    %   mfc_fhat_alg1_decoupled_block and put it in the external PID.
    %
    %   The price of coupling is that F_hat is no longer the plant's own
    %   dynamics: it is a closed-loop quantity that moves when you retune Kp
    %   and is meaningless with the loop open.
    %
    %   SIGN CONVENTION: e = y - y_sp (measurement minus filtered setpoint),
    %   matching mfc_siso.step and mfc_siso_core. Feeding y_sp - y inverts
    %   the folded pole and the loop goes unstable.
    %
    %   Ports
    %     In : err     tracking error e = y - y_sp
    %          u_prev  command that actually reached the plant over the LAST
    %                  sample; in a loop assembled from separate blocks this
    %                  must come through an explicit unit delay
    %          t       clock, starting with the run
    %          alpha   optional live gain (overrides the mask parameter)
    %     Out: F_hat
    %
    %   The math is mfc_fhat_algebraic_first_order with b_fold = -Kp; this
    %   class only maps parameters and Simulink state onto it.
    %
    %   See also mfc_fhat_algebraic_first_order, mfc_fhat_alg1_decoupled_block,
    %   mfc_fhat_alg2_coupled_block, mfc_command_block, mfc_siso_core.

    properties
        % alpha Ultra-local model input gain (ignored if the live alpha input is enabled)
        alpha = 1
        % Kp Proportional gain, folded in as b_fold = -Kp (closed-loop pole at -Kp). Do NOT also apply it externally.
        Kp = 25
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
        function obj = mfc_fhat_alg1_coupled_block(varargin)
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

            % Coupled: error-driven, with the pole s + Kp folded in.
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
            names = {'err', 'u_prev', 't'};
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
            icon = sprintf('F-hat algebraic\n1st order, coupled\nKp=%g folded (no Kd)', obj.Kp);
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_fhat_alg1_coupled_block', ...
                'Title', 'MFC F-hat: algebraic 1st order, coupled', ...
                'Text', sprintf(['Growing-window operational-calculus estimator for ', ...
                    'dot_e = F + alpha*u - Kp*e, driven by the tracking error ', ...
                    'e = y - y_sp.\n\n', ...
                    'The pole s + Kp is folded into F_hat, so NO external P feedback is ', ...
                    'required -- tie the command block''s fb input to Ground, or to a ', ...
                    'discrete integrator for Ki. Applying Kp again outside doubles it.\n\n', ...
                    'There is no Kd at first order: the model has no derivative term to ', ...
                    'fold one into, so this block cannot damp. Use the decoupled ', ...
                    'first-order block if the plant needs D action.']));
        end
    end
end
