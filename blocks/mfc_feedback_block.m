classdef mfc_feedback_block < matlab.System
    % mfc_feedback_block  The explicit feedback law (MFC pipeline stage 3).
    %
    %   Everything the F estimator did NOT already absorb:
    %
    %     decoupled : fb = Kd*dot_err + Kp*err + Ki*int_err   (full iPD(I))
    %     coupled   : fb =                       Ki*int_err
    %
    %   In the coupled structure Kp is folded into F_hat at either model
    %   order and Kd is folded at second order, so only the integral term
    %   remains explicit here. Set 'coupled' to match the estimator block
    %   you wired up -- applying Kp twice (folded AND explicit) doubles the
    %   proportional gain and is the most common wiring mistake.
    %
    %   The integral is trapezoidal. dot_err is a raw backward difference
    %   and is the noise-sensitive term of the whole controller; it is
    %   output even in the coupled case (where fb ignores it) so it can be
    %   logged, or filtered with an mfc_smoother_block and fed back in as an
    %   external derivative.
    %
    %   ANTI-WINDUP. The optional 'freeze' input discards this sample's
    %   integration and holds the previous integral -- wire it from the
    %   'sat' output of mfc_command_filter_block. That signal must pass
    %   through a UNIT DELAY to avoid an algebraic loop, so a loop assembled
    %   from separate blocks freezes ONE SAMPLE LATER than the all-in-one
    %   mfc_siso_core does. This is expected; see
    %   Knowledge/block-library-signal-flow.md.
    %
    %   Ports
    %     In : err  (+ optional freeze)
    %     Out: fb, int_err, dot_err
    %
    %   The math is mfc_siso.feedback; this class only maps parameters and
    %   Simulink state onto it.
    %
    %   See also mfc_siso, mfc_command_block, mfc_command_filter_block,
    %   mfc_siso_core.

    properties
        % Kp Proportional gain (ignored when coupled -- it is folded into F_hat there)
        Kp = 25
        % Kd Derivative gain (ignored when coupled -- folded into F_hat at 2nd order, unavailable at 1st)
        Kd = 10
        % Ki Integral gain, always explicit (0 disables integral action)
        Ki = 0
    end

    properties (Nontunable)
        % Ts Sample time [s] (fixes the block's discrete rate)
        Ts = 0.01
    end

    properties (Nontunable, Logical)
        % coupled Estimator already absorbed Kp (and Kd at 2nd order): apply only Ki here
        coupled = false
        % use_freeze_input Add the freeze input port (anti-windup handshake from the command filter)
        use_freeze_input = false
    end

    properties (DiscreteState)
        err_km1
        int_err
    end

    methods
        function obj = mfc_feedback_block(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)

        function [fb, int_err, dot_err] = stepImpl(obj, err, varargin)
            [fb, int_err, dot_err] = mfc_siso.feedback( ...
                err, obj.err_km1, obj.int_err, obj.Ts, ...
                obj.Kp, obj.Kd, obj.Ki, obj.coupled);

            if obj.use_freeze_input && varargin{1} ~= 0
                int_err = obj.int_err;      % anti-windup: discard this sample
                if obj.coupled
                    fb = obj.Ki*int_err;
                else
                    fb = obj.Kd*dot_err + obj.Kp*err + obj.Ki*int_err;
                end
            end

            obj.err_km1 = err;
            obj.int_err = int_err;
        end

        function resetImpl(obj)
            obj.err_km1 = 0;
            obj.int_err = 0;
        end

        function [sz, dt, cp] = getDiscreteStateSpecificationImpl(~, ~)
            sz = [1 1];  dt = 'double';  cp = false;
        end

        function sts = getSampleTimeImpl(obj)
            % Fixed discrete rate (do NOT inherit): the trapezoidal integral
            % and the backward difference assume one advance per Ts.
            sts = createSampleTime(obj, 'Type', 'Discrete', 'SampleTime', obj.Ts);
        end

        % ---- Ports -------------------------------------------------------
        function num = getNumInputsImpl(obj), num = 1 + obj.use_freeze_input; end
        function varargout = getInputNamesImpl(obj)
            if obj.use_freeze_input
                varargout = {'err', 'freeze'};
            else
                varargout = {'err'};
            end
        end

        function num = getNumOutputsImpl(~), num = 3; end
        function varargout = getOutputNamesImpl(~)
            varargout = {'fb', 'int_err', 'dot_err'};
        end
        function varargout = getOutputSizeImpl(~)
            varargout = {[1 1], [1 1], [1 1]};
        end
        function varargout = getOutputDataTypeImpl(~)
            varargout = {'double', 'double', 'double'};
        end
        function varargout = isOutputComplexImpl(~)
            varargout = {false, false, false};
        end
        function varargout = isOutputFixedSizeImpl(~)
            varargout = {true, true, true};
        end

        function icon = getIconImpl(obj)
            if obj.coupled
                icon = sprintf('feedback\ncoupled: Ki only\nKi = %g', obj.Ki);
            else
                icon = sprintf('feedback\niPD(I) explicit\nKp=%g Kd=%g Ki=%g', ...
                               obj.Kp, obj.Kd, obj.Ki);
            end
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_feedback_block', ...
                'Title', 'MFC Feedback Law', ...
                'Text', sprintf(['Applies whatever the F estimator did not absorb.\n\n', ...
                    'Decoupled: fb = Kd*dot_err + Kp*err + Ki*int_err. ', ...
                    'Coupled: fb = Ki*int_err only, because Kp (and Kd at 2nd order) ', ...
                    'are folded into F_hat by the estimator -- set the coupled flag to ', ...
                    'match your estimator or the proportional gain is applied twice.\n\n', ...
                    'The optional freeze input implements anti-windup from the command ', ...
                    'filter''s sat output; it needs a unit delay to break the algebraic ', ...
                    'loop, so it acts one sample later than mfc_siso_core does.']));
        end
    end
end
