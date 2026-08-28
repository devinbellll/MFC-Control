function golden_capture_mimo(varargin)
%GOLDEN_CAPTURE_MIMO Write the reference 2x2 MIMO traces to tests/golden/*.csv.
%
%   Run ONCE, from a known-good tree, before refactoring:
%       octave --no-gui -q tests/golden_capture_mimo.m
%
%   golden_capture_mimo('name1', 'name2', ...) writes ONLY those variants.
%   Use it when ADDING a variant, so the existing files are not silently
%   rewritten -- rewriting a golden is how a regression gets blessed by
%   accident.
%
%   Regenerating these files is an explicit decision to change the reference
%   behaviour -- do not re-run it to make TEST_GOLDEN_MIMO pass.
%
%   Values are written with %.17g, which round-trips a double exactly, so
%   TEST_GOLDEN_MIMO can demand bit-identical equality against its own
%   capture (it also cross-checks 'mimo_diag' against the existing
%   tests/golden/2nd_decoupled_alg.csv, to a tolerance).
%
%   See also MFC_GOLDEN_TRACE_MIMO, TEST_GOLDEN_MIMO, GOLDEN_CAPTURE.

here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, 'golden');
if ~exist(outdir, 'dir'), mkdir(outdir); end

[traces, names] = mfc_golden_trace_mimo();

if isempty(varargin)
    wanted = names;
else
    wanted = varargin;
    for i = 1:numel(wanted)
        assert(any(strcmp(names, wanted{i})), ...
            'golden_capture_mimo: unknown variant ''%s''.', wanted{i});
    end
end

n_written = 0;
for c = 1:numel(traces)
    if ~any(strcmp(wanted, names{c})), continue; end
    n_written = n_written + 1;
    T = traces{c};
    n = (size(T, 2) - 1) / 5;
    cols = mimo_cols(n);

    f = fullfile(outdir, [names{c} '.csv']);
    fid = fopen(f, 'w');
    fprintf(fid, '%s\n', strjoin(cols, ','));
    fmt = [repmat('%.17g,', 1, numel(cols) - 1), '%.17g\n'];
    for k = 1:size(T, 1)
        fprintf(fid, fmt, T(k, :));
    end
    fclose(fid);
    fprintf('  wrote %-16s  %d samples  max|u|=%.6g\n', ...
            [names{c} '.csv'], size(T, 1), max(max(abs(T(:, 1:n)))));
end
fprintf('golden_capture_mimo: %d variants written to %s\n', n_written, outdir);
end


function cols = mimo_cols(n)
%MIMO_COLS Column header names for an n-channel MIMO trace csv.
    base = {'u', 'F_hat', 'sp_filt', 'err', 'u_raw'};
    cols = {};
    for b = 1:numel(base)
        for i = 1:n
            cols{end+1} = sprintf('%s%d', base{b}, i); %#ok<AGROW>
        end
    end
    cols{end+1} = 'valid';
end
