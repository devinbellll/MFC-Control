classdef mfc_siso_core < matlab.System
    % mfc_siso_core  Second-order Model-Free Control (MFC) SISO controller.
    %
    %   Direct port of mfc_core.c (Paparazzi firmware).  Folded-estimator
    %   variant: the control poles are folded into the F estimator numerator.
    %
    %   Ultra-local model:  ddot_y = F + alpha*u
    %   Closed-loop poles:
    %     use_kd=false: double pole at -kp  (a = -2*kp,    b = -kp^2)
    %     use_kd=true:  omega_n=kp, zeta=kd (a = -2*kd*kp, b = -kp^2)
    %
    %   Reference trajectory filter (time_trajec = T, dimensionless):
    %     sp_traj[k] = (sp + (2T^2+2T)*sp[k-1] + (-T^2)*sp[k-2]) / (T^2+2T+1)
    %   Same IIR structure for the F estimator (int_window = W).
    %
    %   Command EMA filter: U = (raw + (command_filter-1)*U_prev) / command_filter
    %   command_filter=1 disables it (pass-through).
    %
    %   Output saturation to [u_min, u_max] is applied only when
    %   use_control_sat is enabled; otherwise U passes through unclamped.
    %
    %   NOTE on use_ref_filter=false: the feedforward ddot_sp is still computed
    %   from the raw setpoint history (sp_traj = setpoint passed through).
    %   The firmware zeroes it instead — that is a firmware bug.
    %
    %   Inputs : setpoint, measure, t     Outputs: U, F_k, sp_traj, err

    properties
        alpha          = 10      % Model gain: ddot_y = F + alpha*u
        kp             = 1       % Proportional gain (natural frequency omega_n if use_kd)
        kd             = 0       % Damping ratio zeta; only used when use_kd = true
        time_trajec    = 50      % Reference trajectory filter time constant [samples]
        int_window     = 5       % F-estimator averaging window [samples]
        command_filter = 1       % Output EMA filter constant; 1 = disabled (pass-through)
        u_min          = -9600   % Minimum command output; applied only if use_control_sat
        u_max          =  9600   % Maximum command output; applied only if use_control_sat
    end

    properties (Nontunable)
        Ts = 0.002                 % Sample time [s]
    end

    properties (Nontunable, Logical)
        use_control_sat = false    % Clamp U to [u_min, u_max]
        use_kd          = false    % false: double pole at -kp | true: 2nd-order pole (omega_n=kp, zeta=kd)
        use_ref_filter  = true     % false: raw setpoint pass-through (e.g. guidance axes) | true: filtered trajectory
    end

    properties (DiscreteState)
        setpoint_trajec_km1
        setpoint_trajec_km2
        error_km1
        error_km2
        estimator_num_km1
        estimator_num_km2
        estimator_den_km1
        estimator_den_km2
        command_km1
    end

    methods
        function obj = mfc_siso_core(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)

        function [U, F_k, sp_traj, err] = stepImpl(obj, setpoint, measure, t)
            Ts = obj.Ts;
            T  = obj.time_trajec;
            W  = obj.int_window;
            al = obj.alpha;

            % 1) Reference trajectory filter (or raw pass-through)
            if obj.use_ref_filter
                sp_traj = (setpoint + (2*T^2 + 2*T)*obj.setpoint_trajec_km1 + ...
                           (-T^2)*obj.setpoint_trajec_km2) / (T^2 + 2*T + 1);
            else
                sp_traj = setpoint;
            end
            % ddot always from sp_traj history (use_ref_filter=false uses raw setpoint history)
            ddot_sp = (sp_traj - 2*obj.setpoint_trajec_km1 + obj.setpoint_trajec_km2) / Ts^2;

            % 2) Control poles
            if obj.use_kd
                a = -2 * obj.kd * obj.kp;
            else
                a = -2 * obj.kp;
            end
            b = -obj.kp^2;

            err = measure - sp_traj;

            sde   = -(t*err - (t-Ts)*obj.error_km1) / Ts;
            s2d2e =  (t^2*err - 2*(t-Ts)^2*obj.error_km1 + (t-2*Ts)^2*obj.error_km2) / Ts^2;
            sd2e  =  (t^2*err - (t-Ts)^2*obj.error_km1) / Ts;
            de    = -t*err;
            d2e   =  t^2*err;
            d2u   =  t^2*obj.command_km1;

            num = 2*err + 4*sde + s2d2e - a*(2*de + sd2e) - b*d2e - al*d2u;
            den = t^2;

            % 3) F estimator filter
            Fnum = (num + (2*W^2 + 2*W)*obj.estimator_num_km1 + (-W^2)*obj.estimator_num_km2) ...
                   / (W^2 + 2*W + 1);
            Fden = (den + (2*W^2 + 2*W)*obj.estimator_den_km1 + (-W^2)*obj.estimator_den_km2) ...
                   / (W^2 + 2*W + 1);

            F_k = 0;
            if (Fden ~= 0) && (t > 0.1)
                F_k = Fnum / Fden;
            end

            % 4) Command + EMA filter + saturation
            raw_cmd = -F_k/al + ddot_sp/al;
            U = (raw_cmd + (obj.command_filter - 1)*obj.command_km1) / obj.command_filter;
            if obj.use_control_sat
                if U > obj.u_max; U = obj.u_max; end
                if U < obj.u_min; U = obj.u_min; end
            end

            % Shift history
            obj.setpoint_trajec_km2 = obj.setpoint_trajec_km1;
            obj.setpoint_trajec_km1 = sp_traj;
            obj.error_km2           = obj.error_km1;
            obj.error_km1           = err;
            obj.estimator_num_km2   = obj.estimator_num_km1;
            obj.estimator_num_km1   = Fnum;
            obj.estimator_den_km2   = obj.estimator_den_km1;
            obj.estimator_den_km1   = Fden;
            obj.command_km1         = U;
        end

        function resetImpl(obj)
            obj.setpoint_trajec_km1 = 0;
            obj.setpoint_trajec_km2 = 0;
            obj.error_km1           = 0;
            obj.error_km2           = 0;
            obj.estimator_num_km1   = 0;
            obj.estimator_num_km2   = 0;
            obj.estimator_den_km1   = 0;
            obj.estimator_den_km2   = 0;
            obj.command_km1         = 0;
        end

        function [sz, dt, cp] = getDiscreteStateSpecificationImpl(~, ~)
            sz = [1 1];
            dt = 'double';
            cp = false;
        end

        function sts = getSampleTimeImpl(obj)
            sts = createSampleTime(obj, 'Type', 'Discrete', 'SampleTime', obj.Ts);
        end

        function [o1, o2, o3, o4] = getOutputSizeImpl(~)
            o1 = [1 1];  o2 = [1 1];  o3 = [1 1];  o4 = [1 1];
        end
        function [o1, o2, o3, o4] = getOutputDataTypeImpl(~)
            o1 = 'double';  o2 = 'double';  o3 = 'double';  o4 = 'double';
        end
        function [o1, o2, o3, o4] = isOutputComplexImpl(~)
            o1 = false;  o2 = false;  o3 = false;  o4 = false;
        end
        function [o1, o2, o3, o4] = isOutputFixedSizeImpl(~)
            o1 = true;  o2 = true;  o3 = true;  o4 = true;
        end

        function num = getNumInputsImpl(~)
            num = 3;
        end

        function icon = getIconImpl(~)
            icon = 'mfc\_siso\_core';
        end
    end
end
