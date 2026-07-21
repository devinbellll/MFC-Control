function golden_capture()
%GOLDEN_CAPTURE Write the reference MFC traces to tests/golden/*.csv.
%
%   Run ONCE, from a known-good tree, before refactoring:
%       octave --no-gui -q tests/golden_capture.m
%
%   Regenerating these files is an explicit decision to change the reference
%   behaviour -- do not re-run it to make TEST_GOLDEN pass.
%
%   Values are written with %.17g, which round-trips a double exactly, so
%   TEST_GOLDEN can demand bit-identical equality rather than a tolerance.
%
%   See also MFC_GOLDEN_TRACE, TEST_GOLDEN.

here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, 'golden');
if ~exist(outdir, 'dir'), mkdir(outdir); end

[traces, names] = mfc_golden_trace();

cols = 'u,F_hat,sp_filt,err,u_raw,valid';
for c = 1:numel(traces)
    f = fullfile(outdir, [names{c} '.csv']);
    fid = fopen(f, 'w');
    fprintf(fid, '%s\n', cols);
    T = traces{c};
    for k = 1:size(T, 1)
        fprintf(fid, '%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n', T(k, :));
    end
    fclose(fid);
    fprintf('  wrote %-22s  %d samples  max|u|=%.6g\n', ...
            [names{c} '.csv'], size(T, 1), max(abs(T(:, 1))));
end
fprintf('golden_capture: %d variants written to %s\n', numel(traces), outdir);
end
