classdef mfc_smoother_block < matlab.System
    % mfc_smoother_block  Second-order IIR smoother (MFC pipeline stage 1).
    %
    %   Unity-DC-gain, critically damped, second-order IIR low-pass with a
    %   repeated real pole at z = W/(W+1), i.e. a memory of roughly W
    %   samples. W = 0 is an exact pass-through.
    %
    %   ONE BLOCK, FOUR ROLES. This is deliberately the same class wherever
    %   the MFC pipeline smooths a signal, so changing the smoother math is
    %   a single edit that propagates everywhere:
    %     1. reference trajectory filter (the "input smoother") -- enable
    %        the derivative outputs to get the feedforward dot_sp/ddot_sp
    %     2. algebraic estimator NUMERATOR filter   (dissected mode)
    %     3. algebraic estimator DENOMINATOR filter (dissected mode)
    %     4. F_hat post-filter -- mainly useful after mfc_fhat_window_block,
    %        which has no internal smoothing of its own
    %
    %   Roles 2 and 3 must use the SAME window: the algebraic estimators
    %   divide numerator by denominator, and only identical smoothing of
    %   both leaves the ratio unbiased.
    %
    %   Ports
    %     In : x
    %     Out: x_filt, and (optional) dot_x, ddot_x -- backward finite
    %          differences of the FILTERED signal history, so they are the
    %          derivatives of what actually leaves the block.
    %
    %   The math is mfc_siso.ref_traj / mfc_iir_smoother; this class only
    %   maps parameters and Simulink state onto them.
    %
    %   See also mfc_siso, mfc_iir_smoother, mfc_siso_core.

    properties
        % window Smoother memory W [samples]; 0 = exact pass-through
        window = 10
    end

    properties (Nontunable)
        % Ts Sample time [s] (fixes the block's discrete rate)
        Ts = 0.01
    end

    properties (Nontunable, Logical)
        % use_filter Apply the smoother (false = pass the input straight through)
        use_filter = true
        % output_derivatives Add the dot_x and ddot_x output ports
        output_derivatives = false
    end

    properties (DiscreteState)
        x_km1
        x_km2
    end

    methods
        function obj = mfc_smoother_block(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)

        function varargout = stepImpl(obj, x)
            [x_filt, dot_x, ddot_x] = mfc_siso.ref_traj( ...
                x, obj.x_km1, obj.x_km2, obj.Ts, obj.window, obj.use_filter);

            obj.x_km2 = obj.x_km1;
            obj.x_km1 = x_filt;

            varargout{1} = x_filt;
            if obj.output_derivatives
                varargout{2} = dot_x;
                varargout{3} = ddot_x;
            end
        end

        function resetImpl(obj)
            obj.x_km1 = 0;
            obj.x_km2 = 0;
        end

        function [sz, dt, cp] = getDiscreteStateSpecificationImpl(~, ~)
            sz = [1 1];  dt = 'double';  cp = false;
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
            if obj.use_filter
                icon = sprintf('IIR smoother\nW = %g', obj.window);
            else
                icon = sprintf('IIR smoother\n(bypassed)');
            end
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_smoother_block', ...
                'Title', 'MFC IIR Smoother', ...
                'Text', sprintf(['Unity-DC-gain critically damped 2nd-order IIR low-pass ', ...
                    '(repeated pole at W/(W+1), memory ~W samples; W = 0 passes through).\n\n', ...
                    'Used as the reference-trajectory filter, as the numerator and ', ...
                    'denominator filters of a dissected algebraic estimator, and as an ', ...
                    'F_hat post-filter. Enable the derivative outputs when driving the ', ...
                    'feedforward path of mfc_command_block.']));
        end
    end
end
