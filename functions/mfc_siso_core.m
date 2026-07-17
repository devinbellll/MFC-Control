classdef mfc_siso_core < matlab.System
    % mfc_siso_core  Model-Free Control (MFC) SISO controller block.
    %
    %   Manager block for the MFC SISO controller family. The algorithm
    %   variant is selected on the mask along three axes:
    %
    %     Model order          First order  (dot_y  = F + alpha*u)
    %                          Second order (ddot_y = F + alpha*u)
    %     Controller structure Coupled   -- the F estimator is driven by the
    %                          tracking error; at second order the closed-loop
    %                          polynomial s^2 + Kd*s + Kp is folded into the
    %                          estimate, at first order Kp stays explicit.
    %                          Decoupled -- the F estimator is driven by the
    %                          pure measurement (true-plant F) and an explicit
    %                          iP/iPD(I) feedback law stabilizes the error.
    %     Estimator type       Algebraic (growing-window operational-calculus
    %                          recursion) or Sliding window (fixed-length
    %                          Simpson-quadrature integral; decoupled only).
    %
    %   The numerics live in plain, individually testable functions:
    %   mfc_siso_step (per-sample law and dispatch), mfc_fhat_algebraic_
    %   first/second_order, mfc_fhat_sliding_window, mfc_window_kernel and
    %   mfc_iir_smoother. This class only maps mask parameters and Simulink
    %   ports/states onto those functions.
    %
    %   Command law (Kd is always active at second order):
    %     2nd order: u = ( -F_hat + ddot_sp - fb ) / alpha
    %                fb = Kd*dot_err + Kp*err + Ki*int_err   (decoupled)
    %                fb =                       Ki*int_err   (coupled)
    %     1st order: u = ( -F_hat + dot_sp - Kp*err - Ki*int_err ) / alpha
    %   followed by an optional output EMA filter and saturation with
    %   integrator-freeze anti-windup.
    %
    %   Ports
    %     In : y_sp (setpoint), y_m (measurement), t (clock)
    %          u_applied -- optional (External applied-command input): the
    %          command that actually reached the plant during the last
    %          sample, e.g. from a measured actuator or after an external
    %          limiter. It replaces the internally stored previous command
    %          in the estimator, the window buffers and the EMA filter.
    %          alpha -- optional (Live alpha input): overrides the alpha
    %          mask parameter sample by sample.
    %     Out: u, F_hat, sp_filt, err, u_raw (pre-filter/saturation command),
    %          F_valid (1 once the estimator is past its startup hold).
    %
    %   See also mfc_siso_config, mfc_siso_step, mfc_siso_init.

    % ---- Algorithm variant (mask dropdowns) -----------------------------
    properties (Nontunable)
        % model_order Order of the ultra-local model
        model_order = 'Second order (ddot_y = F + alpha*u)'
        % controller_structure Where the stabilizing dynamics live
        controller_structure = 'Coupled (error-driven estimator, poles folded)'
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
            'Sliding window (Simpson quadrature)'});
    end

    % ---- Tuning ----------------------------------------------------------
    properties
        % alpha Ultra-local model input gain (ignored if live alpha input is enabled)
        alpha = 1
        % Kp Proportional gain (2nd order: s^2 + Kd*s + Kp; double pole at -p -> Kp = p^2. 1st order: pole at -Kp)
        Kp = 25
        % Kd Derivative gain (2nd order only; double pole at -p -> Kd = 2p)
        Kd = 10
        % Ki Integral gain (0 disables integral action)
        Ki = 0
        % ref_filter_window Reference trajectory filter memory [samples]
        ref_filter_window = 10
        % command_filter Output EMA filter constant; 1 = disabled (pass-through)
        command_filter = 1
        % u_min Lower command limit (only if saturation enabled)
        u_min = -600
        % u_max Upper command limit (only if saturation enabled)
        u_max = 600
    end

    properties (Nontunable)
        % Ts Sample time [s] (fixes the block's discrete rate)
        Ts = 0.01
        % est_filter_window Estimator memory [samples] (algebraic: num/den smoother; sliding window: window length, rounded to even). Nontunable: sizes the window buffers and quadrature kernel.
        est_filter_window = 10
        % est_hold_time Algebraic estimator held at zero until t exceeds this [s]
        est_hold_time = 0.1
    end

    % ---- Options (mask checkboxes) ---------------------------------------
    properties (Nontunable, Logical)
        % use_ref_filter Filter the setpoint into a smooth reference trajectory
        use_ref_filter = true
        % use_control_sat Clamp the command to [u_min, u_max] (with anti-windup)
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
        cfg     % configuration struct built by mfc_siso_config in setupImpl
    end

    methods
        function obj = mfc_siso_core(varargin)
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
            cfg = mfc_siso_config( ...
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
                n     = obj.est_filter_window;
                n     = n + mod(n, 2);                % even (Simpson)
                n_buf = n + 1;
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
            % mfc_siso_config validates gains and rejects the undefined
            % coupled + sliding-window combination, so mask errors surface
            % at once.
            buildConfig(obj);
        end

        function setupImpl(obj)
            obj.cfg = buildConfig(obj);
        end

        function [u, F_hat, sp_filt, err, u_raw, F_valid] = stepImpl(obj, setpoint, measure, t, varargin)
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
            [out, state] = mfc_siso_step(setpoint, measure, t, u_prev, alpha_k, c, state);
            unpackState(obj, state);

            u       = out.u;
            F_hat   = out.F_hat;
            sp_filt = out.sp_filt;
            err     = out.err;
            u_raw   = out.u_raw;
            F_valid = double(out.est_valid);
        end

        function resetImpl(obj)
            % Zero the state directly instead of via mfc_siso_init: code
            % generation types the discrete states from these assignments,
            % so the buffers need a full-size zeros() with a codegen-constant
            % length (mfc_siso_init sizes them from a run-time cfg value).
            n_buf = bufferLength(obj);
            obj.sp_filt_km1  = 0;
            obj.sp_filt_km2  = 0;
            obj.z_km1        = 0;
            obj.z_km2        = 0;
            obj.num_filt_km1 = 0;
            obj.num_filt_km2 = 0;
            obj.den_filt_km1 = 0;
            obj.den_filt_km2 = 0;
            obj.y_buf        = zeros(n_buf, 1);
            obj.u_buf        = zeros(n_buf, 1);
            obj.err_km1      = 0;
            obj.int_err      = 0;
            obj.u_km1        = 0;
        end

        function [sz, dt, cp] = getDiscreteStateSpecificationImpl(obj, name)
            switch name
                case {'y_buf', 'u_buf'}
                    sz = [bufferLength(obj), 1];
                otherwise
                    sz = [1 1];
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
            num = 6;
        end

        function varargout = getOutputNamesImpl(~)
            varargout = {'u', 'F_hat', 'sp_filt', 'err', 'u_raw', 'F_valid'};
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
            if strncmp(obj.model_order, 'First', 5), o = '1st order'; else, o = '2nd order'; end
            if strncmp(obj.controller_structure, 'Coupled', 7), s = 'coupled'; else, s = 'decoupled'; end
            if strncmp(obj.estimator_type, 'Algebraic', 9), e = 'algebraic'; else, e = 'sliding window'; end
            icon = sprintf('MFC SISO\n%s, %s\n%s', o, s, e);
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_siso_core', ...
                'Title', 'MFC SISO Controller', ...
                'Text', sprintf(['Model-Free Control SISO controller.\n\n', ...
                    'Estimates the ultra-local model term F online and generates the ', ...
                    'command u = (-F_hat + reference feedforward - feedback)/alpha. ', ...
                    'Select the model order (first/second), the controller structure ', ...
                    '(coupled: error-driven estimator with folded poles; decoupled: ', ...
                    'measurement-driven estimator with explicit iP/iPD feedback) and ', ...
                    'the estimator type (algebraic growing window, or sliding-window ', ...
                    'Simpson quadrature -- decoupled only).\n\n', ...
                    'Optional ports: u_applied feeds back the command that actually ', ...
                    'reached the plant (e.g. after an external limiter or from a ', ...
                    'measured actuator); alpha overrides the mask parameter live.']));
        end

        function groups = getPropertyGroupsImpl
            variantSection = matlab.system.display.Section( ...
                'Title', 'Algorithm variant', ...
                'PropertyList', {'model_order', 'controller_structure', 'estimator_type'});
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
