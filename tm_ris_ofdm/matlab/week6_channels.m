%WEEK6_CHANNELS (#4) 3GPP TDL channels, delay spread, larger surfaces, near field.
%   (a) TDL-A / TDL-C / TDL-D (38.901), delay spread 30 ns (as Iudice et al.),
%       100 ns and 300 ns. Two estimator settings: our fixed 0.3 us delay window,
%       and an assumption-free window equal to the NR cyclic prefix (2.34 us).
%   (b) larger surfaces: Q = 64, 256, 1024 with TDL-C, 30 ns.
%   (c) near field: Q = 1024 (1.3 m aperture, Fraunhofer ~40 m), TDL-D with the LoS
%       path computed from exact element distances, user at 3, 10, 30 m.
%   S1, D = 0.75, K = 8, M = 256, Np = 2, timing-only est_general. Gap to perfect H.

clear; rng(3100);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
K = 8; Np = 2;
base = [ones(1, 6), 2*ones(1, 2)];
[F03, l03] = dpss_basis(bins, df, 0.3e-6, 16);
Lcp = ceil(M*df*2.34e-6) + 12; [Fcp, lcp] = dpss_basis(bins, df, 2.34e-6, Lcp);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
snrs = [20 30]; ntrial = 10; R = struct();

%% (a) TDL models and delay spread
models = {'A', 'C', 'D'}; dss = [30e-9 100e-9 300e-9];
for im = 1:3
    for id = 1:3
        gap = zeros(ntrial, 2, 2);                       % [win 0.3us, win CP] x snr
        rp_all = zeros(ntrial, 2);
        for t = 1:ntrial
            [As, Ad] = tdl38901(16, bins, df, fc, models{im}, dss(id));
            H = coupling_matrix(bins, fc, df, repmat(base, 16, 1), As, Ad, P, M-1);
            [g, rp] = gaps(H, T, {F03, Fcp}, {l03, lcp}, Np, evalb, snrs, 16);
            gap(t, :, :) = g; rp_all(t, :) = rp;
        end
        key = sprintf('a_TDL%s_ds%dns', models{im}, round(dss(id)*1e9));
        R.(key) = struct('gap_bits', squeeze(mean(gap, 1)), 'rate_perfect', mean(rp_all, 1));
        fprintf('(a) TDL-%s DS %3d ns: gap @20 win0.3/winCP %.2f/%.2f | @30 %.2f/%.2f (perfect %.2f/%.2f)\n', ...
            models{im}, round(dss(id)*1e9), R.(key).gap_bits(1,1), R.(key).gap_bits(2,1), ...
            R.(key).gap_bits(1,2), R.(key).gap_bits(2,2), R.(key).rate_perfect);
    end
end

%% (b) larger surfaces
for Q = [64 256 1024]
    nt = 5; gap = zeros(nt, 2, 2); rp_all = zeros(nt, 2);
    for t = 1:nt
        [As, Ad] = tdl38901(Q, bins, df, fc, 'C', 30e-9);
        H = coupling_matrix(bins, fc, df, repmat(base, Q, 1), As, Ad, P, M-1);
        [g, rp] = gaps(H, T, {F03, Fcp}, {l03, lcp}, Np, evalb, snrs, Q);
        gap(t, :, :) = g; rp_all(t, :) = rp;
    end
    key = sprintf('b_Q%d', Q);
    R.(key) = struct('gap_bits', squeeze(mean(gap, 1)), 'rate_perfect', mean(rp_all, 1));
    fprintf('(b) Q = %4d (TDL-C): gap @20 %.2f/%.2f | @30 %.2f/%.2f\n', Q, R.(key).gap_bits(:,1), R.(key).gap_bits(:,2));
end

%% (c) near field, Q = 1024, TDL-D LoS from exact distances
for dist = [3 10 30]
    nt = 5; gap = zeros(nt, 2, 2); rp_all = zeros(nt, 2);
    for t = 1:nt
        [As, Ad] = tdl38901(1024, bins, df, fc, 'D', 30e-9, dist);
        H = coupling_matrix(bins, fc, df, repmat(base, 1024, 1), As, Ad, P, M-1);
        [g, rp] = gaps(H, T, {F03, Fcp}, {l03, lcp}, Np, evalb, snrs, 1024);
        gap(t, :, :) = g; rp_all(t, :) = rp;
    end
    key = sprintf('c_nearfield_%dm', dist);
    R.(key) = struct('gap_bits', squeeze(mean(gap, 1)), 'rate_perfect', mean(rp_all, 1));
    fprintf('(c) near field, user at %2d m: gap @20 %.2f/%.2f | @30 %.2f/%.2f\n', dist, R.(key).gap_bits(:,1), R.(key).gap_bits(:,2));
end
fid = fopen(fullfile(out, 'week6_channels.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function [g, rp] = gaps(H, T, Fs, lams, Np, evalb, snrs, Q)
% Gap (bits) between perfect-H MMSE and MMSE with the estimate, for each delay window.
M = size(H, 1); g = zeros(numel(Fs), 2); rp = zeros(1, 2);
Pc = mean(sum(abs(H(:, evalb)).^2, 1));
for is = 1:2
    s2 = Pc / 10^(snrs(is)/10);
    S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
    Y = H*S + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
    rp(is) = mmse_rate(H, H, s2, evalb);
    for iw = 1:numel(Fs)
        Hh = est_general(Y, S, T, ones(Q,1), Fs{iw}, lams{iw}, s2, 8, 20);
        g(iw, is) = rp(is) - mmse_rate(Hh, H, s2, evalb);
    end
end
end

function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
