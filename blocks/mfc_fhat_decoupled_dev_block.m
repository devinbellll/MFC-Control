classdef mfc_fhat_decoupled_dev_block < matlab.System
    % mfc_fhat_decoupled_dev_block  F ESTIMATOR development bench: every
    % decoupled, non-Riachy variant on one mask.
    %
    %   One estimator block whose ports never change, covering the whole
    %   decoupled non-Riachy grid:
    %
    %     Model order      First order  (dot_y  = F + alpha*u)
    %                      Second order (ddot_y = F + alpha*u)
    %     Estimator type   Algebraic (growing-window operational calculus)
    %                      Sliding window (fixed-length FIR taps)
    %     Channel count n  n = 1 is the scalar estimator and takes the
    %                      scalar code path; n > 1 makes y, u_prev and F_hat
    %                      n-by-1 with a square n-by-n alpha.
    %
    %   That is four estimator classes plus their four vector twins reduced
    %   to two dropdowns and a width -- so sweeping the grid is a parameter
    %   change, not a rewiring job. F_hat is the only output either way, so
    %   nothing downstream of this block has to change when the dropdown
    %   does.
    %
    %   DECOUPLED ONLY, and that is the point of the block rather than a
    %   limitation of it. All four variants share one wiring diagram
    %   (measurement in, true-plant F_hat out, explicit feedback added
    %   outside), which is exactly what makes them interchangeable behind
    %   fixed ports. The coupled estimators take err instead of y and carry
    %   the gains they fold, so they are genuinely different blocks and are
    %   deliberately NOT reachable from here. Neither is Riachy's trick: it
    %   estimates Fk = F + Kd*dot_y from an auxiliary output, so its F_hat
    %   means something else and its loop needs D = 0 and a different ff --
    %   see mfc_fhat_riachy2_block.
    %
    %   THIS IS A BENCH, NOT A SHIPPING BLOCK. Once the variant is chosen,
    %   swap in the specific block: its class name and its mask then state
    %   which variant the model runs, instead of hiding it in a dropdown.
    %   The mapping is exactly one to one:
    %
    %       1st + algebraic       mfc_fhat_alg1_decoupled_block
    %       2nd + algebraic       mfc_fhat_alg2_decoupled_block
    %       either + window       mfc_fhat_window_block
    %       ... and the _mimo_block twin of each, when n > 1
    %
    %   For the whole loop rather than the estimator alone -- smoother,
    %   feedback, command, saturation -- the equivalent bench is
    %   mfc_mimo_core, which also reaches the coupled variants.
    %
    %   Ports (fixed, whatever the mask says)
    %     In : y       plant measurement                       (n-by-1)
    %          u_prev  command that actually reached the plant over the LAST
    %                  sample, through an explicit unit delay  (n-by-1)
    %          t       clock, starting with the run            (scalar)
    %          alpha   optional live gain                      (n-by-n)
    %     Out: F_hat   true-plant lumped dynamics              (n-by-1)
    %
    %   The math is mfc_fhat_algebraic_first_order / _second_order (both
    %   with the folds at zero) or mfc_fhat_sliding_window; this class only
    %   maps parameters and Simulink state onto them.
    %
    %   See also mfc_fhat_alg1_decoupled_block, mfc_fhat_alg2_decoupled_block,
    %   mfc_fhat_window_block, mfc_fhat_alg1_decoupled_mimo_block,
    %   mfc_fhat_alg2_decoupled_mimo_block, mfc_fhat_window_mimo_block,
    %   mfc_mimo_core, mfc_command_mimo_block.

    properties
        % alpha Ultra-local model input gain, square n-by-n (ignored if the live alpha input is enabled)
        alpha = 1
    end

    properties (Nontunable)
        % n Number of channels (alpha is n-by-n, signals are n-by-1; n = 1 is the scalar estimator)
        n = 1
        % model_order Order of the ultra-local model
        model_order = 'Second order (ddot_y = F + alpha*u)'
        % estimator Estimation method
        estimator = 'Algebraic (growing window)'
        % Ts Sample time [s] (fixes the block's discrete rate)
        Ts = 0.01
        % est_filter_window Algebraic only: internal num/den smoother memory [samples]
        est_filter_window = 10
        % est_hold_time Algebraic only: F_hat held at zero until t exceeds this [s]
        est_hold_time = 0.1
        % window_samples Sliding window only: window length [intervals]; Tw = window_samples*Ts
        window_samples = 10
    end

    properties (Hidden, Constant)
        model_orderSet = matlab.system.StringSet({ ...
            'First order (dot_y = F + alpha*u)', ...
            'Second order (ddot_y = F + alpha*u)'});
        estimatorSet = matlab.system.StringSet({ ...
            'Algebraic (growing window)', ...
            'Sliding window (FIR taps)'});
    end

    properties (Nontunable, Logical)
        % use_live_alpha Add the alpha input port (overrides the alpha parameter)
        use_live_alpha = false
    end

    properties (DiscreteState)
        % algebraic estimators
        z_km1
        z_km2
        num_filt_km1
        num_filt_km2
        den_filt_km1
        den_filt_km2
        % sliding-window estimator (1-row placeholders when algebraic)
        y_buf
        u_buf
    end

    properties (Access = private)
        kernel     % quadrature kernel built by mfc_siso.window_kernel in setupImpl
    end

    methods
        function obj = mfc_fhat_decoupled_dev_block(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = private)
        function tf = isAlgebraic(obj)
            tf = strncmp(obj.estimator, 'Algebraic', 9);
        end

        function ord = orderNum(obj)
            if strncmp(obj.model_order, 'First', 5), ord = 1; else, ord = 2; end
        end

        function n_buf = bufferLength(obj)
            % 1-row placeholder when algebraic: the dead sliding-window
            % branch is still compiled under code generation, so the buffers
            % must exist and be typed either way. Built only from Nontunable
            % properties, so it is a compile-time constant.
            if isAlgebraic(obj)
                n_buf = 1;
            else
                n_buf = obj.window_samples + 1;
            end
        end
    end

    methods (Access = protected)

        function setupImpl(obj)
            % Always built: harmless when algebraic, and keeping it
            % unconditional keeps the kernel a fixed type under codegen.
            obj.kernel = mfc_siso.window_kernel(orderNum(obj), obj.window_samples, obj.Ts);
        end

        function F_hat = stepImpl(obj, y, u_prev, t, varargin)
            if obj.use_live_alpha
                alpha_k = varargin{1};
            else
                alpha_k = obj.alpha;
            end

            state = struct( ...
                'z_km1',        obj.z_km1, ...
                'z_km2',        obj.z_km2, ...
                'num_filt_km1', obj.num_filt_km1, ...
                'num_filt_km2', obj.num_filt_km2, ...
                'den_filt_km1', obj.den_filt_km1, ...
                'den_filt_km2', obj.den_filt_km2, ...
                'y_buf',        obj.y_buf, ...
                'u_buf',        obj.u_buf);

            % Decoupled: measurement-driven, nothing folded (the folds are
            % what the coupled blocks exist for).
            if isAlgebraic(obj)
                if orderNum(obj) == 2
                    [F_hat, state] = mfc_fhat_algebraic_second_order( ...
                        y, u_prev, alpha_k, t, obj.Ts, obj.est_filter_window, ...
                        obj.est_hold_time, 0, 0, state);
                else
                    [F_hat, state] = mfc_fhat_algebraic_first_order( ...
                        y, u_prev, alpha_k, t, obj.Ts, obj.est_filter_window, ...
                        obj.est_hold_time, 0, state);
                end
            else
                [F_hat, state] = mfc_fhat_sliding_window( ...
                    y, u_prev, alpha_k, t, obj.kernel, state);
            end

            obj.z_km1        = state.z_km1;
            obj.z_km2        = state.z_km2;
            obj.num_filt_km1 = state.num_filt_km1;
            obj.num_filt_km2 = state.num_filt_km2;
            obj.den_filt_km1 = state.den_filt_km1;
            obj.den_filt_km2 = state.den_filt_km2;
            obj.y_buf        = state.y_buf;
            obj.u_buf        = state.u_buf;
        end

        function resetImpl(obj)
            obj.z_km1        = zeros(obj.n, 1);
            obj.z_km2        = zeros(obj.n, 1);
            obj.num_filt_km1 = zeros(obj.n, 1);
            obj.num_filt_km2 = zeros(obj.n, 1);
            obj.den_filt_km1 = 0;   % denominator is t or t^2: scalar, shared
            obj.den_filt_km2 = 0;
            % One COLUMN per channel, with a codegen-constant length: code
            % generation types the discrete states from these assignments.
            n_buf     = bufferLength(obj);
            obj.y_buf = zeros(n_buf, obj.n);
            obj.u_buf = zeros(n_buf, obj.n);
        end

        function validatePropertiesImpl(obj)
            if ~obj.use_live_alpha && ~isequal(size(obj.alpha), [obj.n obj.n])
                error('mfc:mimo:alphaSize', ...
                    'alpha must be %d-by-%d to match n.', obj.n, obj.n);
            end
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
            dt = 'double';  cp = false;
        end

        function sts = getSampleTimeImpl(obj)
            % Fixed discrete rate (do NOT inherit): the backward-difference
            % recursion and the window buffers both assume they advance
            % exactly once per Ts.
            sts = createSampleTime(obj, 'Type', 'Discrete', 'SampleTime', obj.Ts);
        end

        % ---- Ports -------------------------------------------------------
        function num = getNumInputsImpl(obj), num = 3 + obj.use_live_alpha; end
        function varargout = getInputNamesImpl(obj)
            names = {'y', 'u_prev', 't'};
            if obj.use_live_alpha, names{end+1} = 'alpha'; end
            varargout = names;
        end

        function num = getNumOutputsImpl(~), num = 1; end
        function varargout = getOutputNamesImpl(~),    varargout = {'F_hat'}; end
        function varargout = getOutputSizeImpl(obj),   varargout = {[obj.n 1]}; end
        function varargout = getOutputDataTypeImpl(~), varargout = {'double'}; end
        function varargout = isOutputComplexImpl(~),   varargout = {false}; end
        function varargout = isOutputFixedSizeImpl(~), varargout = {true}; end

        function icon = getIconImpl(obj)
            if isAlgebraic(obj)
                kind = 'algebraic';
            else
                kind = sprintf('window, Tw = %g s', obj.window_samples*obj.Ts);
            end
            icon = sprintf('F-hat dev bench\n%d%s order, decoupled\n%s\nn = %d', ...
                orderNum(obj), suffix(obj), kind, obj.n);
        end

        function s = suffix(obj)
            if orderNum(obj) == 1, s = 'st'; else, s = 'nd'; end
        end
    end

    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('mfc_fhat_decoupled_dev_block', ...
                'Title', 'MFC F-hat: decoupled development bench', ...
                'Text', sprintf(['Every DECOUPLED, non-Riachy F estimator behind one ', ...
                    'set of ports: model order (first/second), estimator (algebraic ', ...
                    'growing window or sliding-window FIR) and channel count on the ', ...
                    'mask, F_hat out.\n\n', ...
                    'Use it to sweep the grid without rewiring. Once the variant is ', ...
                    'chosen, swap in the specific block -- its name then states which ', ...
                    'variant the model runs. The coupled estimators are not reachable ', ...
                    'here (they take err and carry the gains they fold), and neither is ', ...
                    'Riachy''s trick (its F_hat is Fk = F + Kd*dot_y, and its loop needs ', ...
                    'D = 0 and a different ff).\n\n', ...
                    'F_hat is the true-plant lumped dynamics: add explicit feedback into ', ...
                    'the command block''s fb input.']));
        end

        function groups = getPropertyGroupsImpl
            variantSection = matlab.system.display.Section( ...
                'Title', 'Variant', ...
                'PropertyList', {'n', 'model_order', 'estimator'});
            gainSection = matlab.system.display.Section( ...
                'Title', 'Model gain', ...
                'PropertyList', {'alpha', 'use_live_alpha'});
            algSection = matlab.system.display.Section( ...
                'Title', 'Algebraic settings', ...
                'PropertyList', {'est_filter_window', 'est_hold_time'});
            winSection = matlab.system.display.Section( ...
                'Title', 'Sliding-window settings', ...
                'PropertyList', {'window_samples'});
            execSection = matlab.system.display.Section( ...
                'Title', 'Execution', ...
                'PropertyList', {'Ts'});

            groups = [ ...
                matlab.system.display.SectionGroup('Title', 'Estimator', ...
                    'Sections', [variantSection, gainSection]), ...
                matlab.system.display.SectionGroup('Title', 'Numerics', ...
                    'Sections', [algSection, winSection, execSection])];
        end
    end
end
