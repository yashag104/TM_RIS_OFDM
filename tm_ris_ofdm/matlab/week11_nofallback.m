%WEEK11_NOFALLBACK Review G3b (pre-registered 2026-10-22): est_general with and without the pilot-fit fallback,
%   fresh channels, shared-path 30 ns / TDL-A 100 ns / TDL-C 100 ns, Np 2/4/8, 30 channels, 20/30 dB.
clear; rng(11900);
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; K = 8;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
chans = {'phys30', 'TDLA100', 'TDLC100'}; Nps = [2 4 8]; snrs = [20 30]; ntrial = 30;
R = struct('Np', Nps, 'snr_db', snrs, 'ntrial', ntrial);
for ic = 1:3
    if ic == 1, tw = 0.3e-6; L = 16; r = 8; else, tw = 1e-6; L = ceil(M*df*tw) + 10; r = 16; end
    [F, lamF] = dpss_basis(bins, df, tw, L);
    gap = nan(ntrial, 2, numel(Nps), 2); fbu = nan(ntrial, numel(Nps), 2);
    for t = 1:ntrial
        if ic == 1, [As, Ad] = phys_channel(Q, bins, df, fc, popt);
        else, [As, Ad] = tdl38901(Q, bins, df, fc, chans{ic}(4), 100e-9); end
        H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1); Pc = mean(sum(abs(H(:, evalb)).^2, 1));
        for is = 1:2
            s2 = Pc / 10^(snrs(is)/10); Sall = qpsk(M, 8); Yall = H*Sall + noise(M, 8, s2); rp = mmse_rate(H, H, s2, evalb);
            for ip = 1:numel(Nps)
                [Hf, inf0] = est_general_ab(Yall(:, 1:Nps(ip)), Sall(:, 1:Nps(ip)), T, ones(Q,1), F, lamF, s2, r, 20, struct('fallback', true));
                gap(t, 1, ip, is) = rp - mmse_rate(Hf, H, s2, evalb);
                gap(t, 2, ip, is) = rp - mmse_rate(inf0.Hals, H, s2, evalb); fbu(t, ip, is) = inf0.used_fallback;
            end
        end
    end
    C = struct('gap_with_without', squeeze(mean(gap, 1)), 'fallback_frac', squeeze(mean(fbu, 1)));
    for ip = 1:numel(Nps)
        for is = 1:2
            [m, h] = tci(gap(:, 2, ip, is) - gap(:, 1, ip, is));
            C.(sprintf('without_minus_with_Np%d_%ddB', Nps(ip), snrs(is))) = [m h];
            fprintf('%s Np=%d %d dB: with %.3f without %.3f | diff %+.3f +-%.3f | fallback used %.0f%%\n', chans{ic}, Nps(ip), snrs(is), ...
                mean(gap(:, 1, ip, is)), mean(gap(:, 2, ip, is)), m, h, 100*mean(fbu(:, ip, is)));
        end
    end
    C.gap_trials = gap; R.(chans{ic}) = C;
    fid = fopen(fullfile(out, 'week11_nofallback.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);
end
function [m, h] = tci(x)
n = numel(x); nu = n - 1; tq = sqrt(nu*(1/betaincinv(0.05, nu/2, 0.5) - 1));
m = mean(x); h = tq*std(x)/sqrt(n);
end
function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
function W = noise(M, n, s2), W = sqrt(s2/2) * (randn(M,n) + 1j*randn(M,n)); end
function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
