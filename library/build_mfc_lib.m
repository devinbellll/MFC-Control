function build_mfc_lib(varargin)
%BUILD_MFC_LIB Generate the MFC Simulink block library.
%
%   >> setup;  build_mfc_lib
%   >> build_mfc_lib('save', false)     % build and leave it open, do not save
%
%   Creates library/mfc_lib.mdl containing a MATLABSystem block for every
%   class in blocks/, laid out in pipeline order, plus a "Dissected chain"
%   subsystem showing an algebraic estimator taken apart into raw
%   numerator/denominator, two smoothers and a divide block.
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
    % {class, display name, [x y], annotation}
    B = { ...
      'mfc_siso_core',            'MFC SISO Controller',     [40   40], 'All-in-one (start here)'; ...
      'mfc_smoother_block',       'IIR Smoother',            [340  40], 'Stage 1 / num / den / F post-filter'; ...
      'mfc_fhat_alg1_block',      'F-hat Algebraic 1st',     [40  200], 'Stage 2: growing window'; ...
      'mfc_fhat_alg2_block',      'F-hat Algebraic 2nd',     [340 200], 'Stage 2: growing window'; ...
      'mfc_fhat_window_block',    'F-hat Sliding Window',    [640 200], 'Stage 2: Simpson, decoupled only'; ...
      'mfc_fhat_divide_block',    'F-hat Divide + Hold',     [940 200], 'Closes a dissected estimator'; ...
      'mfc_feedback_block',       'Feedback Law',            [40  360], 'Stage 3: iPD(I) / Ki-only'; ...
      'mfc_command_block',        'Command (inversion)',     [340 360], 'Stage 4: (-F+ff-fb)/alpha'; ...
      'mfc_command_filter_block', 'Command Filter + Sat',    [640 360], 'Stage 5: EMA + clamp'};

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

    % ---- a worked "dissected chain" the user can copy out ---------------
    build_dissected_demo(libname);

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


function build_dissected_demo(libname)
%BUILD_DISSECTED_DEMO A ready-made subsystem showing the estimator taken apart.
%
%   alg2 (expose_raw, internal_filter off) -> num_raw -> smoother -,
%                                          -> den_raw -> smoother -+-> divide -> F_hat
%
%   Both smoothers MUST share the same window; that is the whole point of
%   the arrangement and the easiest thing to get wrong by hand.

    sub = [libname '/Dissected Estimator (example)'];
    add_block('built-in/Subsystem', sub, 'Position', [940 40 1180 140]);

    W = 10;  Ts = 0.01;

    add_block('built-in/Inport',  [sub '/z'],      'Position', [ 30  50  60  64]);
    add_block('built-in/Inport',  [sub '/u_prev'], 'Position', [ 30 110  60 124]);
    add_block('built-in/Inport',  [sub '/t'],      'Position', [ 30 170  60 184]);

    add_block('simulink/User-Defined Functions/MATLAB System', [sub '/estimator'], ...
              'System', 'mfc_fhat_alg2_block', ...
              'Position', [130 40 280 190]);
    set_param([sub '/estimator'], 'expose_raw', 'on', ...
                                  'internal_filter', 'off', ...
                                  'Ts', num2str(Ts));

    add_block('simulink/User-Defined Functions/MATLAB System', [sub '/num filter'], ...
              'System', 'mfc_smoother_block', 'Position', [350 40 470 90]);
    set_param([sub '/num filter'], 'window', num2str(W), 'Ts', num2str(Ts));

    add_block('simulink/User-Defined Functions/MATLAB System', [sub '/den filter'], ...
              'System', 'mfc_smoother_block', 'Position', [350 120 470 170]);
    set_param([sub '/den filter'], 'window', num2str(W), 'Ts', num2str(Ts));

    add_block('simulink/User-Defined Functions/MATLAB System', [sub '/divide'], ...
              'System', 'mfc_fhat_divide_block', 'Position', [540 60 660 160]);

    add_block('built-in/Outport', [sub '/F_hat'], 'Position', [720  80 750  94]);
    add_block('built-in/Outport', [sub '/valid'], 'Position', [720 130 750 144]);

    % estimator ports: in (z, u_prev, t) / out (F_hat, valid, num_raw, den_raw)
    add_line(sub, 'z/1',      'estimator/1');
    add_line(sub, 'u_prev/1', 'estimator/2');
    add_line(sub, 't/1',      'estimator/3');
    add_line(sub, 'estimator/3', 'num filter/1');     % num_raw
    add_line(sub, 'estimator/4', 'den filter/1');     % den_raw
    add_line(sub, 'num filter/1', 'divide/1');
    add_line(sub, 'den filter/1', 'divide/2');
    add_line(sub, 't/1',          'divide/3', 'autorouting', 'on');
    add_line(sub, 'divide/1', 'F_hat/1');
    add_line(sub, 'divide/2', 'valid/1');

    add_block('built-in/Note', [sub '/note'], 'Position', [400 220], 'FontSize', 9, ...
        'Text', sprintf(['Both smoothers must use the SAME window (%d here).\n', ...
                         'Identical smoothing of numerator and denominator is\n', ...
                         'what keeps the ratio unbiased.'], W));
end
