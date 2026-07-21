classdef mfc_command_filter_block < matlab.System
    % mfc_command_filter_block  Output EMA filter and saturation (MFC stage 5).
    %
    %   Exponential moving average on the raw command, then an optional
    %   clamp:
    %
    %       u = ( u_raw + (c - 1)*u_prev ) / c,      c = command_filter
    %       u = min(u_max, max(u_min, u))            if saturation enabled
    %
    %   c = 1 is an exact pass-through. Larger c trades command bandwidth
    %   for smoothness -- useful because the MFC command law feeds ddot_sp
    %   and the estimator's noise straight through 1/alpha.
    %
    %   ANTI-WINDUP HANDSHAKE. The 'sat' output is true on any sample where
    %   the clamp actually bit. Wire it to the 'freeze' input of
    %   mfc_feedback_block so the integrator stops accumulating while the
    %   actuator is pinned. That path needs a UNIT DELAY to avoid an
    %   algebraic loop, so a composed loop freezes ONE SAMPLE LATER than the
    %   all-in-one mfc_siso_core, which does both in the same sample. The
    %   difference is one sample of extra windup at the moment of
    %   saturation; mfc_siso_core remains the reference implementation.
    %
    %   PREVIOUS COMMAND. By default the block feeds back its own last
    %   output. Enable the u_prev input to supply the command that actually
    %   reached the plant instead -- from a measured actuator, or after an
    %   external limiter. Feed the same signal to the estimator's u_prev
    %   port, or the estimator and the filter disagree about what the plant
    %   was actually driven with.
    %
    %   Ports
    %     In : u_raw  (+ optional u_prev)
    %     Out: u, sat
    %
    %   The math is mfc_siso.limit.
    %
    %   See also mfc_siso, mfc_command_block, mfc_feedback_block,
    %   mfc_siso_core.

    properties
        % command_filter Output EMA constant c; 1 = disabled (pass-through)
        command_filter = 1
        % u_min Lower command limit (only if saturation enabled)
        u_min = -600
        % u_max Upper command limit (only if saturation enabled)
        u_max = 600
    end

    properties (Nontunable)
        % Ts Sample time [s] (fixes the block's discrete rate)
        Ts = 0.01
    end

    properties (Nontunable, Logical)
        % use_control_sat Clamp the command to [u_min, u_max] and emit the sat flag
        use_control_sat = false
        % use_external_command Add the u_prev input port (command that actually reached the plant)
        use_external_command = false
    end

    properties (DiscreteState)
        u_km1
    end

    methods
        function obj = mfc_command_filter_block(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)

        function [u, sat] = stepImpl(obj, u_raw, varargin)
            if obj.use_external_command
                u_prev = varargin{1};
            else
                u_prev = obj.u_km1;
            end

            [u, frozen] = mfc_siso.limit(u_raw, u_prev, obj.command_filter, ...
                                         obj.use_control_sat, obj.u_min, obj.u_max);
            obj.u_km1 = u;
            sat = double(frozen);
        end

        function resetImpl(obj)
            obj.u_km1 = 0;
        end

        function [sz, dt, cp] = getDiscreteStateSpecificationImpl(~, ~)
            sz = [1 1];  dt = 'double';  cp = false;
        end

        function sts = getSampleTimeImpl(obj)
            % Fixed discrete rate (do NOT inherit): the EMA assumes it
            % advances exactly once per Ts.
            sts = createSampleTime(obj, 'Type', 'Discrete', 'SampleTime', obj.Ts);
        end

        % ---- Ports -------------------------------------------------------
        function num = getNumInputsImpl(obj), num = 1 + obj.use_external_command; end
        function varargout = getInputNamesImpl(obj)
            if obj.use_external_command
                varargout = {'u_raw', 'u_prev'};
            else
                varargout = {'u_raw'};
            end
        end

        function num = getNumOutputsImpl(~), num = 2; end
        function varargout = getOutputNamesImpl(~), varargout = {'u', 'sat'}; end
        function varargout = getOutputSizeImpl(~),     varargout = {[1 1], [1 1]}; end
        function varargout = getOutputDataTypeImpl(~), varargout = {'double', 'double'}; end
        function varargout = isOutputComplexImpl(~),   varargout = {false, false}; end
        function varargout = isOutputFixedSizeImpl(~), varargout = {true, true}; end

        function icon = getIconImpl(obj)
            if obj.use_control_sat
                icon = sprintf('command filter\nEMA c = %g\nsat [%g, %g]', ...
                               obj.command_filter, obj.u_min, obj.u_max);
            else
                icon = sprintf('command filter\nEMA c = %g\n(no saturation)', ...
                               obj.command_filter);
            end
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_command_filter_block', ...
                'Title', 'MFC Command Filter and Saturation', ...
                'Text', sprintf(['EMA u = (u_raw + (c-1)*u_prev)/c (c = 1 passes ', ...
                    'through), then an optional clamp to [u_min, u_max].\n\n', ...
                    'The sat output is the anti-windup handshake: wire it through a unit ', ...
                    'delay to the freeze input of mfc_feedback_block. Because of that ', ...
                    'delay a composed loop freezes one sample later than mfc_siso_core.\n\n', ...
                    'Enable u_prev to feed back the command that actually reached the ', ...
                    'plant; give the estimator block the same signal.']));
        end
    end
end
