classdef mfc_command_block < matlab.System
    % mfc_command_block  Ultra-local model inversion, the final MFC stage.
    %
    %   The one line that makes this model-free control:
    %
    %       u = ( -F_hat + ff - fb ) / alpha       then optionally clamped
    %
    %   Cancel the estimated lumped dynamics, inject the reference
    %   feedforward, subtract the feedback, scale by the model gain.
    %
    %   WHAT TO WIRE INTO fb depends on which estimator you chose, and that
    %   is the whole coupled/decoupled distinction:
    %
    %     decoupled estimator (y-driven, or the sliding window)
    %         fb <- a stock Discrete PID Controller driven by the tracking
    %               error. All of Kp, Kd, Ki are applied there.
    %
    %     coupled estimator (error-driven, poles folded)
    %         fb <- Ground, because Kp (and Kd at second order) are already
    %               inside F_hat. If you want integral action, wire a stock
    %               Discrete-Time Integrator on the error scaled by Ki --
    %               Ki is the only gain that survives coupling. Do NOT put a
    %               P or D term here; it would be applied twice.
    %
    %   SIGN CONVENTION: the feedback is SUBTRACTED and the error is
    %   e = y - y_sp (measurement minus filtered setpoint), matching
    %   mfc_siso.step. Drive the PID with that error, not y_sp - y.
    %
    %   ff is the feedforward derivative of the reference trajectory:
    %   ddot_sp for a second-order ultra-local model, dot_sp for first
    %   order. Take it from the matching output of an mfc_smoother_block
    %   with derivatives enabled. This is the ONLY place model order enters
    %   the command law.
    %
    %   ALPHA IS A DESIGN PARAMETER, not a plant identification. Too small
    %   over-drives the command, too large under-drives it; the estimator
    %   absorbs the mismatch into F_hat either way, so the loop still works
    %   over a wide range of alpha. It is the main knob to try when a loop
    %   is sluggish or twitchy. Enable the live input to sweep or schedule
    %   it -- but feed the SAME alpha to the estimator block, or F_hat and
    %   this inversion disagree about what the model is.
    %
    %   SATURATION HAS NO ANTI-WINDUP. The clamp is a plain limiter: an
    %   external integrator keeps accumulating while the command is pinned.
    %   With Ki = 0 (the usual MFC case) that is irrelevant. With integral
    %   action AND a real actuator limit, either use a PID block with its
    %   own back-calculation anti-windup, or use mfc_siso_core, which
    %   freezes its integrator in the same sample the clamp bites.
    %
    %   Stateless -- it holds no history and can run at any rate.
    %
    %   Ports
    %     In : F_hat, ff, fb  (+ optional alpha)
    %     Out: u
    %
    %   The math is mfc_siso.command.
    %
    %   See also mfc_siso, mfc_smoother_block, mfc_fhat_alg2_coupled_block,
    %   mfc_fhat_alg2_decoupled_block, mfc_siso_core.

    properties
        % alpha Ultra-local model input gain (ignored if the live alpha input is enabled)
        alpha = 1
        % u_min Lower command limit (only if saturation enabled)
        u_min = -600
        % u_max Upper command limit (only if saturation enabled)
        u_max = 600
    end

    properties (Nontunable, Logical)
        % use_control_sat Clamp the command to [u_min, u_max] (plain limiter, no anti-windup)
        use_control_sat = false
        % use_live_alpha Add the alpha input port (overrides the alpha parameter)
        use_live_alpha = false
    end

    methods
        function obj = mfc_command_block(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)

        function u = stepImpl(obj, F_hat, ff, fb, varargin)
            if obj.use_live_alpha
                alpha_k = varargin{1};
            else
                alpha_k = obj.alpha;
            end

            u = mfc_siso.command(F_hat, ff, fb, alpha_k);

            if obj.use_control_sat
                u = min(obj.u_max, max(obj.u_min, u));
            end
        end

        % ---- Ports -------------------------------------------------------
        function num = getNumInputsImpl(obj), num = 3 + obj.use_live_alpha; end
        function varargout = getInputNamesImpl(obj)
            names = {'F_hat', 'ff', 'fb'};
            if obj.use_live_alpha, names{end+1} = 'alpha'; end
            varargout = names;
        end

        function num = getNumOutputsImpl(~), num = 1; end
        function varargout = getOutputNamesImpl(~),     varargout = {'u'}; end
        function varargout = getOutputSizeImpl(~),      varargout = {[1 1]}; end
        function varargout = getOutputDataTypeImpl(~),  varargout = {'double'}; end
        function varargout = isOutputComplexImpl(~),    varargout = {false}; end
        function varargout = isOutputFixedSizeImpl(~),  varargout = {true}; end

        function icon = getIconImpl(obj)
            if obj.use_live_alpha
                head = sprintf('command\n(-F+ff-fb)/alpha\nalpha: live');
            else
                head = sprintf('command\n(-F+ff-fb)/alpha\nalpha = %g', obj.alpha);
            end
            if obj.use_control_sat
                icon = sprintf('%s\nsat [%g, %g]', head, obj.u_min, obj.u_max);
            else
                icon = head;
            end
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_command_block', ...
                'Title', 'MFC Command (model inversion)', ...
                'Text', sprintf(['u = ( -F_hat + ff - fb ) / alpha, optionally clamped ', ...
                    'to [u_min, u_max].\n\n', ...
                    'fb depends on the estimator: with a DECOUPLED estimator wire a stock ', ...
                    'Discrete PID on the tracking error (e = y - y_sp); with a COUPLED ', ...
                    'estimator wire Ground, since Kp (and Kd at 2nd order) are already ', ...
                    'folded into F_hat -- add only a Ki integrator if you want integral ', ...
                    'action.\n\n', ...
                    'ff is the reference feedforward derivative: ddot_sp for a ', ...
                    'second-order ultra-local model, dot_sp for first order.\n\n', ...
                    'alpha is a design parameter, not a plant identification; feed the ', ...
                    'same alpha to the estimator block. The clamp has no anti-windup.']));
        end
    end
end
