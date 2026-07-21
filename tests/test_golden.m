function test_golden()
%TEST_GOLDEN Assert the MFC controller core still reproduces the reference traces.
%
%   octave --no-gui -q tests/test_golden.m     (also runs in MATLAB)
%
%   Re-runs MFC_GOLDEN_TRACE and demands a BIT-IDENTICAL match against the
%   csv files captured by GOLDEN_CAPTURE from the pre-refactor tree. This is
%   the no-regression gate for the decomposition: any change that alters a
%   single output sample of any supported variant fails here.
%
%   A tolerance is deliberately NOT used. The refactor only moves code
%   between files; it performs no algebraic rearrangement, so the floating
%   point result must be unchanged. If this ever needs a tolerance, the
%   change is bigger than advertised and deserves a look.
%
%   See also MFC_GOLDEN_TRACE, GOLDEN_CAPTURE.

here   = fileparts(mfilename('fullpath'));
golden = fullfile(here, 'golden');

assert(exist(golden, 'dir') == 7, ...
    'test_golden: %s missing -- run golden_capture.m on a known-good tree first.', golden);

[traces, names] = mfc_golden_trace();
cols = {'u', 'F_hat', 'sp_filt', 'err', 'u_raw', 'valid'};

n_fail = 0;
for c = 1:numel(traces)
    f = fullfile(golden, [names{c} '.csv']);
    if exist(f, 'file') ~= 2
        fprintf('  FAIL  %-20s  no golden file (%s)\n', names{c}, f);
        n_fail = n_fail + 1;
        continue;
    end

    ref = dlmread(f, ',', 1, 0);       % skip the header row
    got = traces{c};

    if ~isequal(size(ref), size(got))
        fprintf('  FAIL  %-20s  size %dx%d, expected %dx%d\n', names{c}, ...
                size(got, 1), size(got, 2), size(ref, 1), size(ref, 2));
        n_fail = n_fail + 1;
        continue;
    end

    if isequaln(ref, got)
        fprintf('  PASS  %-20s  %d samples, bit-identical\n', names{c}, size(got, 1));
    else
        n_fail = n_fail + 1;
        bad = find(any(ref ~= got & ~(isnan(ref) & isnan(got)), 2), 1);
        [~, worst] = max(max(abs(ref - got), [], 1));
        fprintf('  FAIL  %-20s  first mismatch at sample %d, worst column %s (max|d|=%.3g)\n', ...
                names{c}, bad, cols{worst}, max(abs(ref(:, worst) - got(:, worst))));
    end
end

fprintf('\ntest_golden: %d/%d variants match\n', numel(traces) - n_fail, numel(traces));
if n_fail > 0
    error('test_golden: %d variant(s) regressed', n_fail);
end
end
