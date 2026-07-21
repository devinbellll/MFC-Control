classdef mfc_command_block < matlab.System
    % mfc_command_block  Ultra-local model inversion (MFC pipeline stage 4).
    %
    %   The one line that makes this model-free control:
    %
    %       u_raw = ( -F_hat + ff - fb ) / alpha
    %
    %   Cancel the estimated lumped dynamics, inject the reference
    %   feedforward, subtract the feedback, and scale by the model gain.
    %
    %   ff is the feedforward derivative of the reference trajectory:
    %   ddot_sp for a second-order ultra-local model, dot_sp for first
    %   order. Take it from the matching output of an mfc_smoother_block
    %   with derivatives enabled. This is the ONLY place model order enters
    %   the command law -- the coupled/decoupled split lives entirely in the
    %   estimator and feedback blocks.
    %
    %   ALPHA IS A DESIGN PARAMETER, not a plant identification. Too small
    %   over-drives the command, too large under-drives it; the estimator
    %   absorbs the mismatch into F_hat either way, so the loop still works
    %   over a wide range of alpha. It is the main knob to try when a loop
    %   is sluggish or twitchy. Enable the live input to sweep or schedule
    %   it -- but feed the SAME alpha to the estimator block, or F_hat and
    %   this inversion disagree about what the model is.
    %
    %   Stateless -- it holds no history and can run at any rate.
    %
    %   Ports
    %     In : F_hat, ff, fb  (+ optional alpha)
    %     Out: u_raw
    %
    %   The math is mfc_siso.command.
    %
    %   See also mfc_siso, mfc_smoother_block, mfc_feedback_block,
    %   mfc_command_filter_block, mfc_siso_core.

    properties
        % alpha Ultra-local model input gain (ignored if the live alpha input is enabled)
        alpha = 1
    end

    properties (Nontunable, Logical)
        % use_live_alpha Add the alpha input port (overrides the alpha parameter)
        use_live_alpha = false
    end

    methods
        function obj = mfc_command_block(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)

        function u_raw = stepImpl(obj, F_hat, ff, fb, varargin)
            if obj.use_live_alpha
                alpha_k = varargin{1};
            else
                alpha_k = obj.alpha;
            end
            u_raw = mfc_siso.command(F_hat, ff, fb, alpha_k);
        end

        % ---- Ports -------------------------------------------------------
        function num = getNumInputsImpl(obj), num = 3 + obj.use_live_alpha; end
        function varargout = getInputNamesImpl(obj)
            names = {'F_hat', 'ff', 'fb'};
            if obj.use_live_alpha, names{end+1} = 'alpha'; end
            varargout = names;
        end

        function num = getNumOutputsImpl(~), num = 1; end
        function varargout = getOutputNamesImpl(~), varargout = {'u_raw'}; end
        function varargout = getOutputSizeImpl(~),      varargout = {[1 1]}; end
        function varargout = getOutputDataTypeImpl(~),  varargout = {'double'}; end
        function varargout = isOutputComplexImpl(~),    varargout = {false}; end
        function varargout = isOutputFixedSizeImpl(~),  varargout = {true}; end

        function icon = getIconImpl(obj)
            if obj.use_live_alpha
                icon = sprintf('command\n(-F+ff-fb)/alpha\nalpha: live');
            else
                icon = sprintf('command\n(-F+ff-fb)/alpha\nalpha = %g', obj.alpha);
            end
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_command_block', ...
                'Title', 'MFC Command (model inversion)', ...
                'Text', sprintf(['u_raw = ( -F_hat + ff - fb ) / alpha.\n\n', ...
                    'ff is the reference feedforward derivative: ddot_sp for a ', ...
                    'second-order ultra-local model, dot_sp for first order -- the only ', ...
                    'place model order enters the command law.\n\n', ...
                    'alpha is a design parameter, not a plant identification; the ', ...
                    'estimator absorbs any mismatch into F_hat. Feed the same alpha to ', ...
                    'the estimator block.']));
        end
    end
end
