%FIG_DATA_COUPLING Coupling matrices for the paper's structure figure (no randomness, no estimation).
%   (a) Reproduction of Iudice et al. Fig. 5: Q = 1, M = 64, binary sequence {c0, c1}, K = 2, ideal single-tap
%       channels (A_s = A_d = 1), OpenRIS cell. Our coupling_matrix drops the T_u/T_s factor, so |H| here equals
%       their plotted (T_s/T_u)|{H}|.
%   (b) This paper's configuration: K = 8, D = 0.75 (pattern [1 1 1 1 1 1 2 2]), same cell and channels.
%   (c) The same with the ideal +-1 element, for the dispersion/non-ideality difference.
%   Writes results/fig_coupling_data.json.

here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; M = 64; bins = -M/2:M/2-1; one = ones(1, M);
R = struct('M', M, 'fc', fc, 'df', df);
H = coupling_matrix(bins, fc, df, [1 2], one, one, P, M-1);
R.iudice_K2 = struct('re', real(H), 'im', imag(H));
H = coupling_matrix(bins, fc, df, [ones(1, 6), 2*ones(1, 2)], one, one, P, M-1);
R.ours_K8 = struct('re', real(H), 'im', imag(H));
H = coupling_matrix(bins, fc, df, [ones(1, 6), 2*ones(1, 2)], one, one, 'ideal', M-1);
R.ideal_K8 = struct('re', real(H), 'im', imag(H));
fid = fopen(fullfile(out, 'fig_coupling_data.json'), 'w'); fprintf(fid, '%s', jsonencode(R)); fclose(fid);
fprintf('done\n');
