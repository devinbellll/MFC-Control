function test_golden_riachy()
%TEST_GOLDEN_RIACHY Assert Riachy's trick reproduces its reference traces.
%
%   octave --no-gui -q tests/test_golden_riachy.m     (also runs in MATLAB)
%
%   Bit-identical comparison against tests/golden/riachy_*.csv, plus one
%   structural cross-check that is worth more than the equality: the
%   auxiliary output Y must be exactly y + Kd*int y at every sample, where y
%   is recovered from the trace as err + sp_filt. If a future change ever
%   makes the block compute its own transform, or slips a forward Euler in
%   where the trapezoid belongs, the Y columns move and this fails before
%   anyone notices the loop drifting.
%
%   The four variants and why each is there: see MFC_GOLDEN_TRACE_RIACHY.
%
%   See also MFC_GOLDEN_TRACE_RIACHY, GOLDEN_CAPTURE_RIACHY, TEST_RIACHY,
%   TEST_GOLDEN_MIMO.

here   = fileparts(mfilename('fullpath'));
golden = fullfile(here, 'golden');

assert(exist(golden, 'dir') == 7, ...
    'test_golden_riachy: %s missing -- run golden_capture_riachy.m on a known-good tree first.', golden);

[traces, names] = mfc_golden_trace_riachy();

n_fail = 0;

% --- 1) bit-identical against the captured csv ---------------------------
for c = 1:numel(traces)
    f = fullfile(golden, [names{c} '.csv']);
    if exist(f, 'file') ~= 2
        fprintf('  FAIL  %-16s  no golden file (%s)\n', names{c}, f);
        n_fail = n_fail + 1;
        continue;
    end

    ref = dlmread(f, ',', 1, 0);
    got = traces{c};

    if ~isequal(size(ref), size(got))
        fprintf('  FAIL  %-16s  size %dx%d, expected %dx%d\n', names{c}, ...
                size(got, 1), size(got, 2), size(ref, 1), size(ref, 2));
        n_fail = n_fail + 1;
        continue;
    end

    if isequaln(ref, got)
        fprintf('  PASS  %-16s  %d samples, bit-identical\n', names{c}, size(got, 1));
    else
        n_fail = n_fail + 1;
        [~, worst] = max(max(abs(ref - got), [], 1));
        fprintf('  FAIL  %-16s  worst column %d (max|d|=%.3g)\n', ...
                names{c}, worst, max(abs(ref(:, worst) - got(:, worst))));
    end
end

% --- 2) Y really is y + Kd*int y, re-derived from the trace --------------
% The gains used by MFC_GOLDEN_TRACE_RIACHY, restated here on purpose: if
% they change there, this check is supposed to fail rather than follow.
Ts   = 0.01;
Kd_s = 8;                                  % 2p, p = 4
Kd_m = [Kd_s 1.5; -0.8 Kd_s];              % non-diagonal, as in the trace

for c = 1:numel(traces)
    T = traces{c};
    n = (size(T, 2) - 1) / 5;
    if n == 1, Kd = Kd_s; else, Kd = Kd_m; end

    Y       = T(:, 2*n + (1:n));
    sp_filt = T(:, 3*n + (1:n));
    err     = T(:, 4*n + (1:n));
    y       = err + sp_filt;                       % err = y - sp_filt

    % Re-run the transform on the recovered measurement
    rs = struct('int_km1', zeros(n,1), 'y_km1', zeros(n,1));
    Yr = zeros(size(Y));
    for k = 1:size(T, 1)
        [Yk, rs] = mfc_riachy_transform(y(k, :).', Kd, Ts, rs);
        Yr(k, :) = Yk.';
    end

    d  = max(max(abs(Y - Yr)));
    ok = d < 1e-9;
    if ~ok, n_fail = n_fail + 1; end
    fprintf('  %s  %-16s  Y = y + Kd*int y re-derived (max|d|=%.3g)\n', ...
            pass_str(ok), names{c}, d);
end

fprintf('\ntest_golden_riachy: %d check(s) failed\n', n_fail);
if n_fail > 0
    error('test_golden_riachy: %d check(s) failed', n_fail);
end
end

function s = pass_str(ok)
    if ok, s = 'PASS'; else, s = 'FAIL'; end
end
