function test_golden_mimo()
%TEST_GOLDEN_MIMO Assert the 2x2 MIMO variants reproduce the reference traces,
%and cross-checks the diagonal-alpha case against the existing SISO trace
%it must reduce to.
%
%   octave --no-gui -q tests/test_golden_mimo.m     (also runs in MATLAB)
%
%   Five variants: the two hand-wired ones described below, plus
%   'mimo_1st_decoupled_alg', 'mimo_2nd_decoupled_win' and
%   'mimo_2nd_coupled_alg', which cover the rest of the n-channel estimator
%   grid and are driven through mfc_siso.config/step -- the same pipeline
%   mfc_mimo_core runs -- on a plant with a genuine matrix input gain. See
%   MFC_GOLDEN_TRACE_MIMO for what each one pins, including why the coupled
%   variant needs identical per-channel dynamics to stay stable.
%
%   Two properties matter most for a firmware port and are what this test
%   pins:
%
%     * 'mimo_diag' uses alpha = eye(2) and two identical channels, so
%       mfc_command_mimo's matrix solve and the estimator's alpha*d2u term
%       must reduce to the same result as the scalar SISO path -- checked
%       directly against tests/golden/2nd_decoupled_alg.csv, not only
%       against its own capture.
%     * 'mimo_cross' uses an off-diagonal alpha and two DIFFERENT plants, so
%       F_hat is genuinely different per channel (num_raw is per-element)
%       while the estimator's den_raw = t^2 stays scalar and shared -- see
%       MFC_GOLDEN_TRACE_MIMO for how the single 'valid' column pins that.
%
%   See also MFC_GOLDEN_TRACE_MIMO, GOLDEN_CAPTURE_MIMO, TEST_GOLDEN.

here   = fileparts(mfilename('fullpath'));
golden = fullfile(here, 'golden');

assert(exist(golden, 'dir') == 7, ...
    'test_golden_mimo: %s missing -- run golden_capture_mimo.m on a known-good tree first.', golden);

[traces, names] = mfc_golden_trace_mimo();

n_fail = 0;

% --- 1) self-consistency: bit-identical against the captured csv ---------
for c = 1:numel(traces)
    f = fullfile(golden, [names{c} '.csv']);
    if exist(f, 'file') ~= 2
        fprintf('  FAIL  %-12s  no golden file (%s)\n', names{c}, f);
        n_fail = n_fail + 1;
        continue;
    end

    ref = dlmread(f, ',', 1, 0);
    got = traces{c};

    if ~isequal(size(ref), size(got))
        fprintf('  FAIL  %-12s  size %dx%d, expected %dx%d\n', names{c}, ...
                size(got, 1), size(got, 2), size(ref, 1), size(ref, 2));
        n_fail = n_fail + 1;
        continue;
    end

    if isequaln(ref, got)
        fprintf('  PASS  %-12s  %d samples, bit-identical\n', names{c}, size(got, 1));
    else
        n_fail = n_fail + 1;
        [~, worst] = max(max(abs(ref - got), [], 1));
        fprintf('  FAIL  %-12s  worst column %d (max|d|=%.3g)\n', ...
                names{c}, worst, max(abs(ref(:, worst) - got(:, worst))));
    end
end

% --- 2) cross-check: diagonal alpha reduces to the SISO golden trace -----
% alpha = eye(2) and two identical channels means the matrix solve in
% mfc_command_mimo, and the alpha*d2u matrix-vector product inside the
% estimator, must degenerate to the plain scalar path -- so channel 1 and
% channel 2 of 'mimo_diag' should each reproduce tests/golden/2nd_decoupled_alg.csv's
% u column. A tolerance is used rather than bit-exact equality: the matrix
% solve (LAPACK) and the vector estimator arithmetic take a different --
% if mathematically equivalent -- path than the scalar SISO code.
siso_f = fullfile(golden, '2nd_decoupled_alg.csv');
assert(exist(siso_f, 'file') == 2, ...
    'test_golden_mimo: %s missing -- required for the diagonal-alpha cross-check.', siso_f);
siso   = dlmread(siso_f, ',', 1, 0);
u_siso = siso(:, 1);

idx = find(strcmp(names, 'mimo_diag'), 1);
T   = traces{idx};
tol = 1e-9;
for ch = 1:2
    d  = max(abs(T(:, ch) - u_siso));
    ok = d < tol;
    if ~ok, n_fail = n_fail + 1; end
    fprintf('  %s  mimo_diag channel %d matches 2nd_decoupled_alg.csv (max|d|=%.3g, tol=%.0e)\n', ...
            pass_str(ok), ch, d, tol);
end

fprintf('\ntest_golden_mimo: %d check(s) failed\n', n_fail);
if n_fail > 0
    error('test_golden_mimo: %d check(s) failed', n_fail);
end
end

function s = pass_str(ok)
    if ok, s = 'PASS'; else, s = 'FAIL'; end
end
