classdef mfc_command_mimo_block < matlab.System
    % mfc_command_mimo_block  Ultra-local model inversion with a MATRIX gain.
    %
    %   The vector-valued twin of mfc_command_block:
    %
    %       alpha * u = -F_hat + ff - fb      ->   u = alpha \ (-F_hat + ff - fb)
    %
    %   F_hat, ff, fb and u are n-by-1; alpha is a square, invertible n-by-n
    %   matrix, so the scalar division becomes a linear solve and the
    %   channels are cross-coupled through alpha alone.
    %
    %   Everything the SISO block says still applies: the feedback is
    %   SUBTRACTED and the error is e = y - y_sp; ff is ddot_sp for a
    %   second-order ultra-local model; alpha is a DESIGN parameter and must
    %   be the same matrix fed to the estimator; the clamp (element-wise,
    %   against scalar or n-by-1 limits) has NO anti-windup.
    %
    %   Stateless -- it holds no history and can run at any rate.
    %
    %   Ports
    %     In : F_hat, ff, fb  (n-by-1)  (+ optional alpha, n-by-n)
    %     Out: u              (n-by-1)
    %
    %   The math is mfc_siso.command_mimo.
    %
    %   See also mfc_command_block, mfc_siso, mfc_fhat_alg2_decoupled_mimo_block.

    properties
        % alpha Ultra-local model input gain, square n-by-n (ignored if the live alpha input is enabled)
        alpha = eye(2)
        % u_min Lower command limit, scalar or n-by-1 (only if saturation enabled)
        u_min = -600
        % u_max Upper command limit, scalar or n-by-1 (only if saturation enabled)
        u_max = 600
    end

    properties (Nontunable)
        % n Number of channels (alpha is n-by-n, signals are n-by-1)
        n = 2
    end

    properties (Nontunable, Logical)
        % use_control_sat Clamp the command to [u_min, u_max] (plain limiter, no anti-windup)
        use_control_sat = false
        % use_live_alpha Add the alpha input port (overrides the alpha parameter)
        use_live_alpha = false
    end

    methods
        function obj = mfc_command_mimo_block(varargin)
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

            u = mfc_siso.command_mimo(F_hat, ff, fb, alpha_k);

            if obj.use_control_sat
                u = min(obj.u_max, max(obj.u_min, u));
            end
        end

        function validatePropertiesImpl(obj)
            if ~obj.use_live_alpha && ~isequal(size(obj.alpha), [obj.n obj.n])
                error('mfc:mimo:alphaSize', ...
                    'alpha must be %d-by-%d to match n.', obj.n, obj.n);
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
        function varargout = getOutputSizeImpl(obj),    varargout = {[obj.n 1]}; end
        function varargout = getOutputDataTypeImpl(~),  varargout = {'double'}; end
        function varargout = isOutputComplexImpl(~),    varargout = {false}; end
        function varargout = isOutputFixedSizeImpl(~),  varargout = {true}; end

        function icon = getIconImpl(obj)
            if obj.use_live_alpha
                head = sprintf('command (MIMO)\nalpha\\(-F+ff-fb)\nalpha: live (%dx%d)', obj.n, obj.n);
            else
                head = sprintf('command (MIMO)\nalpha\\(-F+ff-fb)\nalpha: %dx%d', obj.n, obj.n);
            end
            if obj.use_control_sat
                icon = sprintf('%s\nsat', head);
            else
                icon = head;
            end
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_command_mimo_block', ...
                'Title', 'MFC Command (model inversion, matrix alpha)', ...
                'Text', sprintf(['u = alpha \\ ( -F_hat + ff - fb ), optionally clamped ', ...
                    'element-wise to [u_min, u_max].\n\n', ...
                    'F_hat, ff, fb and u are n-by-1; alpha is a square, invertible ', ...
                    'n-by-n matrix and must be the SAME matrix fed to the estimator.\n\n', ...
                    'fb depends on the estimator: with a DECOUPLED estimator wire an ', ...
                    'explicit feedback law on the tracking error (e = y - y_sp).\n\n', ...
                    'ff is the reference feedforward derivative: ddot_sp for a ', ...
                    'second-order ultra-local model. The clamp has no anti-windup.']));
        end
    end
end
