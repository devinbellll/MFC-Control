classdef mfc_fhat_divide_block < matlab.System
    % mfc_fhat_divide_block  Ratio and startup hold of a dissected algebraic estimator.
    %
    %   The tail end of an algebraic F estimator, split out so the smoothing
    %   in front of it can be replaced:
    %
    %       valid = (den ~= 0) && (t > hold_time)
    %       F_hat = valid ? num/den : 0
    %
    %   Use it to close the dissected chain:
    %
    %       [alg estimator, expose_raw] -- num_raw --> [smoother] --,
    %                                   -- den_raw --> [smoother] --+--> [divide]
    %
    %   with the SAME window on both smoothers. Smoothing numerator and
    %   denominator identically is what keeps the ratio unbiased; it is also
    %   why this stage cannot simply be an F_hat post-filter.
    %
    %   WHY THE HOLD. The denominator of the growing-window estimators is t
    %   (1st order) or t^2 (2nd order), which is zero at t = 0 and tiny just
    %   after. Dividing by it during the startup transient produces an
    %   enormous F_hat that a controller will happily convert into an
    %   enormous command. hold_time simply refuses to divide until the
    %   denominator has grown; 'valid' tells downstream logic when the
    %   estimate became real.
    %
    %   Stateless -- it holds no history and can run at any rate.
    %
    %   Ports
    %     In : num, den, t
    %     Out: F_hat, valid
    %
    %   See also mfc_fhat_alg1_block, mfc_fhat_alg2_block,
    %   mfc_smoother_block, mfc_siso.

    properties
        % hold_time F_hat forced to 0 while t <= this [s]
        hold_time = 0.1
    end

    methods
        function obj = mfc_fhat_divide_block(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)

        function [F_hat, valid] = stepImpl(obj, num, den, t)
            is_valid = (den ~= 0) && (t > obj.hold_time);
            if is_valid
                F_hat = num / den;
            else
                F_hat = 0;
            end
            valid = double(is_valid);
        end

        % ---- Ports -------------------------------------------------------
        function num = getNumInputsImpl(~),  num = 3; end
        function varargout = getInputNamesImpl(~), varargout = {'num', 'den', 't'}; end
        function num = getNumOutputsImpl(~), num = 2; end
        function varargout = getOutputNamesImpl(~), varargout = {'F_hat', 'valid'}; end
        function varargout = getOutputSizeImpl(~)
            varargout = {[1 1], [1 1]};
        end
        function varargout = getOutputDataTypeImpl(~)
            varargout = {'double', 'double'};
        end
        function varargout = isOutputComplexImpl(~)
            varargout = {false, false};
        end
        function varargout = isOutputFixedSizeImpl(~)
            varargout = {true, true};
        end

        function icon = getIconImpl(obj)
            icon = sprintf('F-hat divide\nnum/den\nhold %g s', obj.hold_time);
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_fhat_divide_block', ...
                'Title', 'MFC F-hat: divide and hold', ...
                'Text', sprintf(['Closes a dissected algebraic estimator: divides the ', ...
                    '(separately smoothed) numerator by the denominator and holds the ', ...
                    'result at zero until t exceeds hold_time.\n\n', ...
                    'Both inputs must be smoothed with the SAME window for the ratio to ', ...
                    'stay unbiased. The hold guards the near-zero denominator (t or t^2) ', ...
                    'during the startup transient.']));
        end
    end
end
