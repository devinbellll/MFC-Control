function build_mfc_lib(varargin)
%BUILD_MFC_LIB Generate the MFC Simulink block library.
%
%   >> setup;  build_mfc_lib
%   >> build_mfc_lib('save', false)     % build and leave it open, do not save
%
%   Creates library/mfc_lib.mdl containing a MATLABSystem block for every
%   class in blocks/, laid out in pipeline order.
%
%   The library is deliberately small. Every choice that changes the
%   CONTROL STRUCTURE is a different block, not a parameter you have to
%   remember to set:
%
%     coupled vs decoupled   different estimator blocks, with different
%                            input port names (err vs y) and only the gains
%                            that variant actually folds on the mask
%     model order            different estimator blocks (1st vs 2nd)
%
%   The one exception is mfc_fhat_riachy2_block, whose estimator kind
%   (algebraic or sliding window) IS a mask parameter: both choices give
%   the same ports, the same wiring and the same Fk -- only the numerics
%   differ, so it is a tuning knob and not a structure.
%
%   Anything Simulink already does well is NOT wrapped: the explicit
%   feedback law is a stock Discrete PID Controller into the command
%   block's fb input, and output filtering is a stock Discrete Filter.
%
%   Run this ONCE in MATLAB and commit the generated .mdl. It exists as a
%   script rather than a hand-written model file because a MATLABSystem
%   block's ports and mask are derived from the class at build time -- so
%   the library can always be regenerated from the classes, and the classes
%   stay the single source of truth.
%
%   Requires MATLAB + Simulink. The repo convention is text-format .mdl
%   (not .slx), which is what this saves.
%
%   See also SLBLOCKS, MFC_SISO_CORE.

    p = inputParser;
    addParameter(p, 'save', true, @islogical);
    parse(p, varargin{:});
    do_save = p.Results.save;

    here    = fileparts(mfilename('fullpath'));
    libname = 'mfc_lib';
    target  = fullfile(here, [libname '.mdl']);

    % Start clean -- a stale library open in memory silently wins otherwise.
    close_system(libname, 0);   %#ok<*NASGU>  (no-op if not loaded)

    new_system(libname, 'Library');
    open_system(libname);
    set_param(libname, 'Lock', 'off');

    % ---- blocks, laid out in pipeline order -----------------------------
    % Row 1: the all-in-one controller and the smoother.
    % Row 2: decoupled estimators (F_hat is the true plant dynamics).
    % Row 3: coupled estimators (feedback poles folded into F_hat).
    % Row 4: the model inversion that ends every composed loop.
    % Row 5: the matrix-alpha (MIMO) pair -- same math on n-by-1 signals,
    %        cross-coupled through a square alpha only.
    %
    % {class, display name, [x y], annotation}
    B = { ...
      'mfc_siso_core',                  'MFC SISO Controller',   [ 40  40], 'All-in-one (start here)'; ...
      'mfc_smoother_block',             'IIR Smoother',          [340  40], 'Reference trajectory + ff derivatives, or an F_hat post-filter'; ...
      ...
      'mfc_fhat_alg1_decoupled_block',  'F-hat Alg 1st (decoupled)', [ 40 200], 'in: y  -- add an external PID'; ...
      'mfc_fhat_alg2_decoupled_block',  'F-hat Alg 2nd (decoupled)', [340 200], 'in: y  -- add an external PID'; ...
      'mfc_fhat_window_block',          'F-hat Sliding Window',      [640 200], 'in: y  -- decoupled only, no coupled variant exists'; ...
      'mfc_fhat_riachy2_block',         'F-hat Riachy 2nd',          [940 200], 'in: y  -- Y = y + Kd*int y; PID needs D = 0, ff = ddot_sp + Kd*dot_sp'; ...
      ...
      'mfc_fhat_alg1_coupled_block',    'F-hat Alg 1st (coupled)',   [ 40 380], 'in: err -- Kp folded, no Kd available at 1st order'; ...
      'mfc_fhat_alg2_coupled_block',    'F-hat Alg 2nd (coupled)',   [340 380], 'in: err -- Kp and Kd folded, fb -> Ground'; ...
      ...
      'mfc_command_block',              'Command (inversion)',       [ 40 560], 'u = (-F_hat + ff - fb)/alpha, optional clamp'; ...
      ...
      'mfc_fhat_alg2_decoupled_mimo_block', 'F-hat Alg 2nd (decoupled, matrix alpha)', [ 40 740], 'in: y (n-by-1) -- vector signals, square n-by-n alpha'; ...
      'mfc_fhat_riachy2_mimo_block',        'F-hat Riachy 2nd (matrix gains)',          [340 740], 'in: y (n-by-1) -- Y = y + Kd*int y, square n-by-n Kd and alpha'; ...
      'mfc_command_mimo_block',             'Command (inversion, matrix alpha)',       [640 740], 'u = alpha\(-F_hat + ff - fb), vector signals'};

    for i = 1:size(B, 1)
        cls  = B{i, 1};
        name = B{i, 2};
        pos  = B{i, 3};
        path = [libname '/' name];

        add_block('simulink/User-Defined Functions/MATLAB System', path, ...
                  'System', cls, ...
                  'Position', [pos(1), pos(2), pos(1)+180, pos(2)+90]);

        add_block('built-in/Note', [libname '/note' num2str(i)], ...
                  'Position', [pos(1)+90, pos(2)+105], ...
                  'Text', B{i, 4}, 'FontSize', '9');
    end

    % ---- how to close the loop, stated once, in the library --------------
    add_block('built-in/Note', [libname '/wiring'], 'Position', [700 560], ...
        'FontSize', 10, 'HorizontalAlignment', 'left', 'Text', sprintf([ ...
        'Closing a composed loop\n', ...
        '\n', ...
        'err = y - sp_filt   (measurement minus FILTERED setpoint, this sign)\n', ...
        'u_prev needs a real Unit Delay -- the estimator wants the command\n', ...
        '  that was actually applied over the last sample.\n', ...
        '\n', ...
        'decoupled estimator:  err -> [Discrete PID: Kp, Kd, Ki] -> command/fb\n', ...
        'Riachy estimator:     err -> [Discrete PID: Kp, Ki, D = 0] -> command/fb\n', ...
        '                      ff = ddot_sp + Kd*dot_sp (Kd is already in Y)\n', ...
        'coupled estimator:    command/fb -> Ground\n', ...
        '                      (or a Ki integrator on err, if you want one)']));

    set_param(libname, 'Lock', 'on');

    if do_save
        save_system(libname, target);
        fprintf('build_mfc_lib: wrote %s\n', target);
        fprintf('  Add library/ to the path (setup.m does) and refresh the\n');
        fprintf('  Library Browser with:  >> sl_refresh_customizations\n');
    else
        fprintf('build_mfc_lib: built %s in memory (not saved)\n', libname);
    end
end
