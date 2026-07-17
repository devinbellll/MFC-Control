function cfg = mfc_siso_config(varargin)
%MFC_SISO_CONFIG Build and validate a configuration for MFC_SISO_STEP.
%
%   cfg = MFC_SISO_CONFIG('Name', Value, ...)
%
%   Assembles the configuration struct consumed by MFC_SISO_STEP /
%   MFC_SISO_INIT, applies defaults, validates the variant selection, and
%   precomputes the sliding-window quadrature kernel when needed. This is
%   the single place where an MFC SISO controller variant is fully
%   specified; the Simulink block MFC_SISO_CORE builds its configuration
%   through this function too.
%
%   Variant selection (2 x 2 x 2 grid, minus the undefined corner)
%     'model_order'  : 1 | 2 (default 2)
%         1: ultra-local model  dot_y = F + alpha*u
%         2: ultra-local model ddot_y = F + alpha*u
%     'structure'    : 'coupled' | 'decoupled' (default 'coupled')
%         coupled  : the F estimator is driven by the tracking ERROR. At
%                    2nd order the closed-loop polynomial s^2 + Kd*s + Kp is
%                    folded into the estimate; at 1st order the closed-loop
%                    pole s + Kp is folded (Kd has no derivative room to
%                    fold at 1st order, so it is unused there).
%         decoupled: the F estimator is driven by the pure MEASUREMENT
%                    (true-plant F); feedback is an explicit iP/iPD(I) law,
%                    at either model order.
%     'estimator'    : 'algebraic' | 'sliding_window' (default 'algebraic')
%         algebraic     : growing-window operational-calculus recursion.
%         sliding_window: fixed-length Simpson-quadrature integral
%                         (decoupled only -- coupled+sliding_window errors).
%
%   Tuning
%     'Ts'                : sample time [s] (default 0.01)
%     'alpha'             : ultra-local model input gain (default 1)
%     'Kp', 'Kd', 'Ki'    : feedback gains (defaults 25, 10, 0).
%         2nd order: characteristic polynomial s^2 + Kd*s + Kp
%                    (double pole at -p  <=>  Kd = 2p, Kp = p^2).
%         1st order: single pole at -Kp.
%         Kp is folded into F_hat whenever coupled (either order); explicit
%         whenever decoupled (either order). Kd is folded into F_hat when
%         coupled at 2nd order, unused when coupled at 1st order (no
%         derivative room to fold), and explicit whenever decoupled.
%         Ki is always applied explicitly (0 disables integral action).
%     'ref_filter_window' : reference trajectory filter memory [samples]
%                           (default 10; see MFC_IIR_SMOOTHER)
%     'est_filter_window' : estimator memory [samples] (default 10)
%                           algebraic: num/den IIR smoother memory;
%                           sliding_window: window length (rounded to even)
%     'est_hold_time'     : algebraic estimator held at 0 for t <= this [s]
%                           (default 0.1; ignored by sliding_window, which
%                           holds until its window fills)
%     'command_filter'    : output EMA constant, u = (raw + (c-1)*u_prev)/c;
%                           1 disables (default 1)
%     'use_ref_filter'    : true/false (default true). When false the raw
%                           setpoint passes through; feedforward derivatives
%                           are then finite differences of the raw setpoint.
%     'use_control_sat'   : true/false (default false); clamp to [u_min, u_max]
%                           with integrator freeze anti-windup
%     'u_min', 'u_max'    : saturation limits (defaults -600, 600)
%
%   Output: cfg struct with the fields above (structure/estimator resolved
%   to logicals cfg.coupled / cfg.algebraic) plus cfg.kernel (sliding-window
%   quadrature kernel; also present but unused for algebraic variants, so
%   cfg has a single concrete type under code generation).
%
%   See also MFC_SISO_INIT, MFC_SISO_STEP, MFC_SISO_CORE.

p = struct( ...
    'model_order',       2, ...
    'structure',         'coupled', ...
    'estimator',         'algebraic', ...
    'Ts',                0.01, ...
    'alpha',             1, ...
    'Kp',                25, ...
    'Kd',                10, ...
    'Ki',                0, ...
    'ref_filter_window', 10, ...
    'est_filter_window', 10, ...
    'est_hold_time',     0.1, ...
    'command_filter',    1, ...
    'use_ref_filter',    true, ...
    'use_control_sat',   false, ...
    'u_min',             -600, ...
    'u_max',             600);

assert(mod(numel(varargin), 2) == 0, 'mfc_siso_config: name-value pairs expected.');
for i = 1:2:numel(varargin)
    name = varargin{i};
    assert(isfield(p, name), 'mfc_siso_config: unknown option ''%s''.', name);
    p.(name) = varargin{i+1};
end

% --- validate the variant selection -------------------------------------
assert(p.model_order == 1 || p.model_order == 2, ...
    'mfc_siso_config: model_order must be 1 or 2.');
assert(any(strcmp(p.structure, {'coupled', 'decoupled'})), ...
    'mfc_siso_config: structure must be ''coupled'' or ''decoupled''.');
assert(any(strcmp(p.estimator, {'algebraic', 'sliding_window'})), ...
    'mfc_siso_config: estimator must be ''algebraic'' or ''sliding_window''.');

coupled   = strcmp(p.structure, 'coupled');
algebraic = strcmp(p.estimator, 'algebraic');

assert(~(coupled && ~algebraic), ...
    ['mfc_siso_config: the coupled (error-driven, pole-folded) structure is ', ...
     'only defined for the algebraic estimator. Use structure=''decoupled'' ', ...
     'with estimator=''sliding_window''.']);

% --- validate tuning -----------------------------------------------------
validateattributes(p.Ts, {'numeric'}, {'scalar', 'positive'}, '', 'Ts');
validateattributes(p.command_filter, {'numeric'}, {'scalar', '>=', 1}, '', 'command_filter');
validateattributes(p.ref_filter_window, {'numeric'}, {'scalar', 'nonnegative'}, '', 'ref_filter_window');
validateattributes(p.est_filter_window, {'numeric'}, {'scalar', 'positive'}, '', 'est_filter_window');
assert(p.u_max > p.u_min, 'mfc_siso_config: u_max must exceed u_min.');

% --- precompute the sliding-window kernel --------------------------------
% Computed for EVERY variant: under code generation both estimator branches
% of MFC_SISO_STEP are compiled (cfg.algebraic is a run-time struct field
% there), so cfg.kernel must always be a struct with one concrete field set.
% Algebraic variants never evaluate it; their window value is sanitized so
% any positive smoother memory remains accepted.
if algebraic
    win = max(2, round(p.est_filter_window));
else
    win = p.est_filter_window;
end
kernel = mfc_window_kernel(p.model_order, win, p.Ts);

% Assemble cfg in one pass: MATLAB Coder forbids adding struct fields after
% the struct has been read, so every field must exist before first use.
cfg           = p;
cfg.coupled   = coupled;
cfg.algebraic = algebraic;
cfg.kernel    = kernel;
end
