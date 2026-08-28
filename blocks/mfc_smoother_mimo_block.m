classdef mfc_smoother_mimo_block < matlab.System
    % mfc_smoother_mimo_block  Second-order IIR smoother, n channels.
    %
    %   The vector-valued twin of mfc_smoother_block: the same unity-DC-gain,
    %   critically damped, second-order IIR low-pass (repeated real pole at
    %   z = W/(W+1), memory ~W samples), applied CHANNEL-WISE to an n-by-1
    %   signal. One shared W, n independent recursions -- the smoother is a
    %   scalar filter and nothing about it mixes channels.
    %
    %   TWO ROLES, exactly as in the SISO case:
    %     1. reference trajectory filter -- enable the derivative outputs to
    %        get dot_sp/ddot_sp for mfc_command_mimo_block's ff input
    %     2. F_hat post-filter -- mainly after mfc_fhat_window_mimo_block,
    %        which has no internal smoothing of its own
    %
    %   The algebraic estimators smooth their own numerator and denominator
    %   internally (one shared window, which is what keeps the ratio
    %   unbiased); that is not something you wire up out here.
    %
    %   W = 0 is an exact pass-through, so there is no separate enable flag.
    %
    %   Ports
    %     In : x                          (n-by-1)
    %     Out: x_filt, and (optional) dot_x, ddot_x -- backward finite
    %          differences of the FILTERED history, all n-by-1
    %
    %   The math is mfc_siso.ref_traj / mfc_iir_smoother, both elementwise.
    %
    %   See also mfc_smoother_block, mfc_fhat_window_mimo_block,
    %   mfc_command_mimo_block, mfc_mimo_core.

    properties
        % window Smoother memory W [samples]; 0 = exact pass-through
        window = 10
    end

    properties (Nontunable)
        % n Number of channels (signals are n-by-1)
        n = 2
        % Ts Sample time [s] (fixes the block's discrete rate)
        Ts = 0.01
    end

    properties (Nontunable, Logical)
        % output_derivatives Add the dot_x and ddot_x output ports
        output_derivatives = false
    end

    properties (DiscreteState)
        x_km1
        x_km2
    end

    methods
        function obj = mfc_smoother_mimo_block(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)

        function varargout = stepImpl(obj, x)
            [x_filt, dot_x, ddot_x] = mfc_siso.ref_traj( ...
                x, obj.x_km1, obj.x_km2, obj.Ts, obj.window, true);

            obj.x_km2 = obj.x_km1;
            obj.x_km1 = x_filt;

            varargout{1} = x_filt;
            if obj.output_derivatives
                varargout{2} = dot_x;
                varargout{3} = ddot_x;
            end
        end

        function resetImpl(obj)
            obj.x_km1 = zeros(obj.n, 1);
            obj.x_km2 = zeros(obj.n, 1);
        end

        function [sz, dt, cp] = getDiscreteStateSpecificationImpl(obj, ~)
            sz = [obj.n 1];  dt = 'double';  cp = false;
        end

        function sts = getSampleTimeImpl(obj)
            % Fixed discrete rate (do NOT inherit): the recursion assumes it
            % advances exactly once per Ts.
            sts = createSampleTime(obj, 'Type', 'Discrete', 'SampleTime', obj.Ts);
        end

        % ---- Ports -------------------------------------------------------
        function num = getNumInputsImpl(~),  num = 1; end
        function varargout = getInputNamesImpl(~), varargout = {'x'}; end

        function num = getNumOutputsImpl(obj)
            num = 1 + 2*obj.output_derivatives;
        end
        function varargout = getOutputNamesImpl(obj)
            if obj.output_derivatives
                varargout = {'x_filt', 'dot_x', 'ddot_x'};
            else
                varargout = {'x_filt'};
            end
        end
        function varargout = getOutputSizeImpl(obj)
            varargout = repmat({[obj.n 1]}, 1, getNumOutputsImpl(obj));
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
            if obj.window == 0
                icon = sprintf('IIR smoother (%d ch)\nW = 0 (pass-through)', obj.n);
            else
                icon = sprintf('IIR smoother (%d ch)\nW = %g', obj.n, obj.window);
            end
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_smoother_mimo_block', ...
                'Title', 'MFC IIR Smoother, n channels', ...
                'Text', sprintf(['Unity-DC-gain critically damped 2nd-order IIR low-pass ', ...
                    '(repeated pole at W/(W+1), memory ~W samples), applied channel-wise ', ...
                    'to an n-by-1 signal.\n\n', ...
                    'Used as the reference-trajectory filter and as an F_hat post-filter ', ...
                    '-- the latter matters most after mfc_fhat_window_mimo_block, which ', ...
                    'does no smoothing of its own. Enable the derivative outputs when ', ...
                    'driving mfc_command_mimo_block''s ff input. W = 0 bypasses exactly.']));
        end
    end
end
