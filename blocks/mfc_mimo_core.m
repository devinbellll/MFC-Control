classdef mfc_mimo_core < matlab.System
    % mfc_mimo_core  Model-Free Control controller, n channels, all settings
    % on the mask. The DEVELOPMENT BENCH for the non-Riachy grid.
    %
    %   Everything mfc_siso_core is, with the channel count as a mask
    %   parameter: n-by-1 signals, a square n-by-n alpha inverted by a
    %   linear solve, and scalar-or-matrix gains. The variant is selected on
    %   four axes instead of three:
    %
    %     Model order          First order  (dot_y  = F + alpha*u)
    %                          Second order (ddot_y = F + alpha*u)
    %     Controller structure Coupled (error-driven estimator, the matrix
    %                          polynomial folded in) or Decoupled
    %                          (measurement-driven estimator, explicit
    %                          iP/iPD(I) feedback)
    %     Estimator type       Algebraic (growing window) or Sliding window
    %                          (FIR taps; decoupled only)
    %     Channel count n      n = 1 is the scalar controller and takes the
    %                          scalar code path exactly -- same division,
    %                          same scalar anti-windup test -- so it
    %                          reproduces mfc_siso_core bit for bit.
    %
    %   WHY THIS BLOCK EXISTS, given the library's rule that a structural
    %   choice is a different block. Two reasons, both about development
    %   rather than deployment:
    %
    %     1. It is the reference implementation on the vector side, the way
    %        mfc_siso_core is on the scalar side: one block that runs the
    %        whole mfc_siso.step pipeline, against which a hand-composed
    %        loop can be checked.
    %     2. It makes the decoupled, non-Riachy grid -- 1st/2nd order,
    %        algebraic/sliding window, any n -- reachable by changing a
    %        dropdown, so a variant sweep is a parameter sweep instead of a
    %        rewiring job.
    %
    %   SHIP THE SPECIFIC BLOCKS, NOT THIS ONE. Once a variant is chosen,
    %   the composed loop built from the specific estimator block (whose
    %   ports state the structure) plus stock Simulink parts is the thing to
    %   hand downstream. This block is where you find out which variant you
    %   want. Riachy's trick is deliberately NOT one of the settings: it
    %   changes what F_hat means (Fk = F + Kd*dot_y) and what the
    %   feedforward must carry, so it stays its own block --
    %   mfc_fhat_riachy2_mimo_block.
    %
    %   Command law, per sample:
    %     2nd order: alpha*u = -F_hat + ddot_sp - fb
    %     1st order: alpha*u = -F_hat + dot_sp  - fb
    %     fb = Kd*dot_err + Kp*err + Ki*int_err   (decoupled, either order)
    %     fb =                       Ki*int_err   (coupled: Kp always folded,
    %                                              Kd folded only at 2nd order)
    %   followed by an optional output EMA filter and saturation with
    %   PER-CHANNEL integrator-freeze anti-windup.
    %
    %   Ports
    %     In : y_sp (setpoint, n-by-1), y_m (measurement, n-by-1), t (clock,
    %          scalar)
    %          u_applied -- optional (n-by-1): the command that actually
    %          reached the plant during the last sample; replaces the
    %          internally stored previous command everywhere it is used.
    %          alpha -- optional (n-by-n): overrides the mask parameter.
    %     Out: u, F_hat, sp_filt, err, u_raw  (all n-by-1)
    %
    %   The numerics live in mfc_siso (config/init/step and the stages) and
    %   the plain estimator functions; this class only maps mask parameters
    %   and Simulink ports/states onto them.
    %
    %   FOR THE ESTIMATOR ALONE, use mfc_fhat_decoupled_dev_block: the same
    %   idea one stage down -- order, estimator and width on the mask, F_hat
    %   out, no smoother/feedback/command wrapped around it.
    %
    %   See also mfc_fhat_decoupled_dev_block, mfc_siso_core, mfc_siso,
    %   mfc_fhat_riachy2_mimo_block, mfc_command_mimo_block,
    %   mfc_fhat_window_mimo_block.

    % ---- Algorithm variant (mask dropdowns) -----------------------------
    properties (Nontunable)
        % n Number of channels (signals n-by-1, alpha n-by-n; n = 1 is the scalar controller)
        n = 2
        % model_order Order of the ultra-local model
        model_order = 'Second order (ddot_y = F + alpha*u)'
        % controller_structure Where the stabilizing dynamics live
        controller_structure = 'Decoupled (measurement-driven estimator, explicit iP/iPD)'
        % estimator_type F-hat estimation method
        estimator_type = 'Algebraic (growing window)'
    end

    properties (Hidden, Constant)
        model_orderSet = matlab.system.StringSet({ ...
            'First order (dot_y = F + alpha*u)', ...
            'Second order (ddot_y = F + alpha*u)'});
        controller_structureSet = matlab.system.StringSet({ ...
            'Coupled (error-driven estimator, poles folded)', ...
            'Decoupled (measurement-driven estimator, explicit iP/iPD)'});
        estimator_typeSet = matlab.system.StringSet({ ...
            'Algebraic (growing window)', ...
            'Sliding window (FIR taps)'});
    end

    % ---- Tuning ----------------------------------------------------------
    properties
        % alpha Ultra-local model input gain, square n-by-n (ignored if the live alpha input is enabled)
        alpha = eye(2)
        % Kp Proportional gain, scalar or n-by-n. Folded into F_hat whenever coupled; explicit whenever decoupled.
        Kp = 25
        % Kd Derivative gain, scalar or n-by-n. Folded when coupled at 2nd order; unused when coupled at 1st order; explicit whenever decoupled.
        Kd = 10
        % Ki Integral gain, scalar or n-by-n (0 disables integral action)
        Ki = 0
        % ref_filter_window Reference trajectory filter memory [samples]
        ref_filter_window = 10
        % command_filter Output EMA filter constant; 1 = disabled (pass-through)
        command_filter = 1
        % u_min Lower command limit, scalar or n-by-1 (only if saturation enabled)
        u_min = -600
        % u_max Upper command limit, scalar or n-by-1 (only if saturation enabled)
        u_max = 600
    end

    properties (Nontunable)
        % Ts Sample time [s] (fixes the block's discrete rate)
        Ts = 0.01
        % est_filter_window Estimator memory [samples] (algebraic: num/den smoother; sliding window: window length). Nontunable: sizes the window buffers and quadrature kernel.
        est_filter_window = 10
        % est_hold_time Algebraic estimator held at zero until t exceeds this [s]
        est_hold_time = 0.1
    end

    % ---- Options (mask checkboxes) ---------------------------------------
    properties (Nontunable, Logical)
        % use_ref_filter Filter the setpoint into a smooth reference trajectory
        use_ref_filter = true
        % use_control_sat Clamp the command to [u_min, u_max] (per-channel anti-windup)
        use_control_sat = false
        % use_external_command Add the u_applied input port (feed back the command that actually reached the plant)
        use_external_command = false
        % use_live_alpha Add the alpha input port (overrides the alpha parameter)
        use_live_alpha = false
    end

    properties (DiscreteState)
        sp_filt_km1
        sp_filt_km2
        z_km1
        z_km2
        num_filt_km1
        num_filt_km2
        den_filt_km1
        den_filt_km2
        y_buf
        u_buf
        err_km1
        int_err
        u_km1
    end

    properties (Access = private)
        cfg     % configuration struct built by mfc_siso.config in setupImpl
    end

    methods
        function obj = mfc_mimo_core(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = private)
        function cfg = buildConfig(obj)
            if strncmp(obj.model_order, 'First', 5), order = 1; else, order = 2; end
            if strncmp(obj.controller_structure, 'Coupled', 7)
                structure = 'coupled';
            else
                structure = 'decoupled';
            end
            if strncmp(obj.estimator_type, 'Algebraic', 9)
                estimator = 'algebraic';
            else
                estimator = 'sliding_window';
            end
            cfg = mfc_siso.config( ...
                'n',                 obj.n, ...
                'model_order',       order, ...
                'structure',         structure, ...
                'estimator',         estimator, ...
                'Ts',                obj.Ts, ...
                'alpha',             obj.alpha, ...
                'Kp',                obj.Kp, ...
                'Kd',                obj.Kd, ...
                'Ki',                obj.Ki, ...
                'ref_filter_window', obj.ref_filter_window, ...
                'est_filter_window', obj.est_filter_window, ...
                'est_hold_time',     obj.est_hold_time, ...
                'command_filter',    obj.command_filter, ...
                'use_ref_filter',    obj.use_ref_filter, ...
                'use_control_sat',   obj.use_control_sat, ...
                'u_min',             obj.u_min, ...
                'u_max',             obj.u_max);
        end

        function n_buf = bufferLength(obj)
            % Window buffer length: shared by resetImpl and
            % getDiscreteStateSpecificationImpl so the sizes cannot diverge.
            % Built only from Nontunable properties, so it is a compile-time
            % constant under code generation.
            if strncmp(obj.estimator_type, 'Algebraic', 9)
                n_buf = 1;                            % unused placeholder
            else
                n_buf = obj.est_filter_window + 1;
            end
        end

        function state = packState(obj)
            state = struct( ...
                'sp_filt_km1',  obj.sp_filt_km1, ...
                'sp_filt_km2',  obj.sp_filt_km2, ...
                'z_km1',        obj.z_km1, ...
                'z_km2',        obj.z_km2, ...
                'num_filt_km1', obj.num_filt_km1, ...
                'num_filt_km2', obj.num_filt_km2, ...
                'den_filt_km1', obj.den_filt_km1, ...
                'den_filt_km2', obj.den_filt_km2, ...
                'y_buf',        obj.y_buf, ...
                'u_buf',        obj.u_buf, ...
                'err_km1',      obj.err_km1, ...
                'int_err',      obj.int_err, ...
                'u_km1',        obj.u_km1);
        end

        function unpackState(obj, state)
            obj.sp_filt_km1  = state.sp_filt_km1;
            obj.sp_filt_km2  = state.sp_filt_km2;
            obj.z_km1        = state.z_km1;
            obj.z_km2        = state.z_km2;
            obj.num_filt_km1 = state.num_filt_km1;
            obj.num_filt_km2 = state.num_filt_km2;
            obj.den_filt_km1 = state.den_filt_km1;
            obj.den_filt_km2 = state.den_filt_km2;
            obj.y_buf        = state.y_buf;
            obj.u_buf        = state.u_buf;
            obj.err_km1      = state.err_km1;
            obj.int_err      = state.int_err;
            obj.u_km1        = state.u_km1;
        end
    end

    methods (Access = protected)

        function validatePropertiesImpl(obj)
            % mfc_siso.config validates the gains, the alpha size against n,
            % and rejects the undefined coupled + sliding-window
            % combination, so mask errors surface at once.
            buildConfig(obj);
        end

        function setupImpl(obj)
            obj.cfg = buildConfig(obj);
        end

        function [u, F_hat, sp_filt, err, u_raw] = stepImpl(obj, setpoint, measure, t, varargin)
            % Live tunable-parameter changes must reach the config
            c = obj.cfg;
            c.alpha             = obj.alpha;
            c.Kp                = obj.Kp;
            c.Kd                = obj.Kd;
            c.Ki                = obj.Ki;
            c.ref_filter_window = obj.ref_filter_window;
            c.command_filter    = obj.command_filter;
            c.u_min             = obj.u_min;
            c.u_max             = obj.u_max;

            % Optional ports, in declaration order: u_applied, then alpha
            idx = 1;
            if obj.use_external_command
                u_prev = varargin{idx};  idx = idx + 1;
            else
                u_prev = obj.u_km1;
            end
            if obj.use_live_alpha
                alpha_k = varargin{idx};
            else
                alpha_k = obj.alpha;
            end

            state = packState(obj);
            [out, state] = mfc_siso.step(setpoint, measure, t, u_prev, alpha_k, c, state);
            unpackState(obj, state);

            u       = out.u;
            F_hat   = out.F_hat;
            sp_filt = out.sp_filt;
            err     = out.err;
            u_raw   = out.u_raw;
        end

        function resetImpl(obj)
            % Zero the state directly instead of via mfc_siso.init: code
            % generation types the discrete states from these assignments,
            % so the buffers need a full-size zeros() with codegen-constant
            % dimensions (mfc_siso.init sizes them from run-time cfg values).
            n_buf = bufferLength(obj);
            obj.sp_filt_km1  = zeros(obj.n, 1);
            obj.sp_filt_km2  = zeros(obj.n, 1);
            obj.z_km1        = zeros(obj.n, 1);
            obj.z_km2        = zeros(obj.n, 1);
            obj.num_filt_km1 = zeros(obj.n, 1);
            obj.num_filt_km2 = zeros(obj.n, 1);
            obj.den_filt_km1 = 0;   % denominator is t or t^2: scalar, shared
            obj.den_filt_km2 = 0;
            obj.y_buf        = zeros(n_buf, obj.n);
            obj.u_buf        = zeros(n_buf, obj.n);
            obj.err_km1      = zeros(obj.n, 1);
            obj.int_err      = zeros(obj.n, 1);
            obj.u_km1        = zeros(obj.n, 1);
        end

        % ---- save/load of the locked object ----------------------------
        % cfg is PRIVATE and built in setupImpl, so the base class does not
        % carry it through a save/load of a LOCKED object -- see the same
        % pair in mfc_siso_core for what that breaks (fast restart and
        % array sim()).
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            if isLocked(obj)
                s.cfg = obj.cfg;
            end
        end

        function loadObjectImpl(obj, s, wasLocked)
            if wasLocked
                obj.cfg = s.cfg;
            end
            loadObjectImpl@matlab.System(obj, s, wasLocked);
        end

        function [sz, dt, cp] = getDiscreteStateSpecificationImpl(obj, name)
            switch name
                case {'y_buf', 'u_buf'}
                    sz = [bufferLength(obj), obj.n];
                case {'den_filt_km1', 'den_filt_km2'}
                    sz = [1 1];
                otherwise
                    sz = [obj.n 1];
            end
            dt = 'double';
            cp = false;
        end

        function sts = getSampleTimeImpl(obj)
            % Fixed discrete rate (do NOT inherit): the per-sample recursions
            % assume they advance exactly once per Ts.
            sts = createSampleTime(obj, 'Type', 'Discrete', 'SampleTime', obj.Ts);
        end

        % ---- Ports -------------------------------------------------------
        function num = getNumInputsImpl(obj)
            num = 3 + obj.use_external_command + obj.use_live_alpha;
        end

        function varargout = getInputNamesImpl(obj)
            names = {'y_sp', 'y_m', 't'};
            if obj.use_external_command, names{end+1} = 'u_applied'; end
            if obj.use_live_alpha,       names{end+1} = 'alpha';     end
            varargout = names;
        end

        function num = getNumOutputsImpl(~)
            num = 5;
        end

        function varargout = getOutputNamesImpl(~)
            varargout = {'u', 'F_hat', 'sp_filt', 'err', 'u_raw'};
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
            if strncmp(obj.model_order, 'First', 5), o = '1st order'; else, o = '2nd order'; end
            if strncmp(obj.controller_structure, 'Coupled', 7), s = 'coupled'; else, s = 'decoupled'; end
            if strncmp(obj.estimator_type, 'Algebraic', 9), e = 'algebraic'; else, e = 'sliding window'; end
            icon = sprintf('MFC core (n = %d)\n%s, %s\n%s', obj.n, o, s, e);
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_mimo_core', ...
                'Title', 'MFC Controller, n channels (development bench)', ...
                'Text', sprintf(['Model-Free Control controller with n-by-1 signals and ', ...
                    'a square n-by-n alpha, with the whole non-Riachy variant grid on ', ...
                    'the mask: model order (first/second), structure (coupled: ', ...
                    'error-driven estimator with the matrix polynomial folded in; ', ...
                    'decoupled: measurement-driven estimator with explicit iP/iPD ', ...
                    'feedback), estimator (algebraic growing window or sliding-window ', ...
                    'FIR -- decoupled only) and the channel count.\n\n', ...
                    'n = 1 takes the scalar code path and reproduces mfc_siso_core ', ...
                    'exactly. Use this block to FIND the variant you want; ship the ', ...
                    'specific estimator block plus stock parts, whose ports state the ', ...
                    'structure. Riachy''s trick is not a setting here -- it changes what ', ...
                    'F_hat means, so it stays mfc_fhat_riachy2_mimo_block.\n\n', ...
                    'Saturation freezes the integrator per channel.']));
        end

        function groups = getPropertyGroupsImpl
            variantSection = matlab.system.display.Section( ...
                'Title', 'Algorithm variant', ...
                'PropertyList', {'n', 'model_order', 'controller_structure', 'estimator_type'});
            gainSection = matlab.system.display.Section( ...
                'Title', 'Model and feedback gains', ...
                'PropertyList', {'alpha', 'Kp', 'Kd', 'Ki'});
            filterSection = matlab.system.display.Section( ...
                'Title', 'Filters and estimator', ...
                'PropertyList', {'ref_filter_window', 'est_filter_window', ...
                                 'command_filter', 'est_hold_time', 'use_ref_filter'});
            outputSection = matlab.system.display.Section( ...
                'Title', 'Command saturation', ...
                'PropertyList', {'use_control_sat', 'u_min', 'u_max'});
            interfaceSection = matlab.system.display.Section( ...
                'Title', 'Optional input ports', ...
                'PropertyList', {'use_external_command', 'use_live_alpha'});
            executionSection = matlab.system.display.Section( ...
                'Title', 'Execution', ...
                'PropertyList', {'Ts'});

            groups = [ ...
                matlab.system.display.SectionGroup('Title', 'Algorithm', ...
                    'Sections', [variantSection, gainSection]), ...
                matlab.system.display.SectionGroup('Title', 'Signal conditioning', ...
                    'Sections', [filterSection, outputSection]), ...
                matlab.system.display.SectionGroup('Title', 'Interface', ...
                    'Sections', [interfaceSection, executionSection])];
        end
    end
end
