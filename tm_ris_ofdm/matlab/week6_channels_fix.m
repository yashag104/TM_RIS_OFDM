%WEEK6_CHANNELS_FIX Rich 3GPP channels: delay window, rank and pilot budget.
%   week6_channels showed large gaps for TDL-A/C/D at 100-300 ns delay spread: those
%   channels have 14-24 taps from different directions (our earlier model had 6 paths),
%   so rank 8 and a 0.3 us window are too small. Here:
%     window  tau_w = min(CP, 10 x DS)   (DS-adaptive; DS is a slow, measurable statistic)
%     rank    r in {8, 16, 32}
%     pilots  Np in {2, 4, 8}
%   TDL-A, TDL-C, TDL-D at DS 100 and 300 ns. S1, Q = 16, M = 256, 6 trials, 20/30 dB.
%   Bar: some (r, Np) brings the gap below 0.3 bit/s/Hz at 20 dB for every model at 100 ns.

clear; rng(3700);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3; Tcp = 2.34e-6;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
K = 8; Q = 16;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
snrs = [20 30]; ntrial = 6; ranks = [8 16 32]; Nps = [2 4 8];
models = {'A', 'C', 'D'}; dss = [100e-9 300e-9];
R = struct('ranks', ranks, 'Np', Nps, 'snr_db', snrs);
for im = 1:3
    for id = 1:2
        tw = min(Tcp, 10*dss(id)); L = ceil(M*df*tw) + 10;
        [F, lamF] = dpss_basis(bins, df, tw, L);
        gap = zeros(ntrial, numel(ranks), numel(Nps), 2);
        for t = 1:ntrial
            [As, Ad] = tdl38901(Q, bins, df, fc, models{im}, dss(id));
            H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
            Pc = mean(sum(abs(H(:, evalb)).^2, 1));
            for is = 1:2
                s2 = Pc / 10^(snrs(is)/10); rp = mmse_rate(H, H, s2, evalb);
                for ip = 1:numel(Nps)
                    S = (sign(randn(M,Nps(ip))) + 1j*sign(randn(M,Nps(ip)))) / sqrt(2);
                    Y = H*S + sqrt(s2/2)*(randn(M,Nps(ip)) + 1j*randn(M,Nps(ip)));
                    for ir = 1:numel(ranks)
                        Hh = est_general(Y, S, T, ones(Q,1), F, lamF, s2, ranks(ir), 20);
                        gap(t, ir, ip, is) = rp - mmse_rate(Hh, H, s2, evalb);
                    end
                end
            end
        end
        key = sprintf('TDL%s_ds%dns', models{im}, round(dss(id)*1e9));
        g = squeeze(mean(gap, 1));                          % rank x Np x snr
        R.(key) = struct('window_s', tw, 'L', L, 'gap_bits', g);
        fprintf('%s (window %.2f us, L %d): gap @20 [rows r=8/16/32, cols Np=2/4/8]\n', key, tw*1e6, L);
        disp(round(g(:,:,1), 2));
        fprintf('   @30\n'); disp(round(g(:,:,2), 2));
    end
end
fid = fopen(fullfile(out, 'week6_channels_fix.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
