%WEEK11_CDL_RANK Review D2 + D3 (pre-registered 2026-10-22).
%   D2  Spatially consistent 3GPP CDL-A/C/D at 100 ns (cascade within the CP, as TDL) versus the TDL +
%       random-direction model used so far: gap of ours (window 1 us, L 18, r 16, Np 2/4/8) to perfect CSI,
%       30 channels, 20/30 dB, Q = 16.
%   D3  Rank of the true element-domain matrix [a_1 ... a_Q] (M x Q, user hop, CDL-C 100 ns) versus Q in
%       {16, 64, 256, 1024} and bandwidth: narrow (256 x 30 kHz = 7.7 MHz) and wide (256 tones spread over
%       98.3 MHz), far field without / with beam squint, near field at 3 m and 10 m (exact, frequency-
%       dependent distances on the LoS ray of CDL-D). Effective rank = number of singular values holding
%       99% / 99.9% of the energy; 10 draws.
clear; rng(11500);
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; K = 8;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
tw = 1e-6; L = ceil(M*df*tw) + 10; [F, lamF] = dpss_basis(bins, df, tw, L);
snrs = [20 30]; Nps = [2 4 8]; ntrial = 30; models = {'A', 'C', 'D'};
gap = nan(ntrial, 3, 2, numel(Nps), 2);                    % model x {TDL, CDL} x Np x snr
for t = 1:ntrial
    for im = 1:3
        for ik = 1:2
            if ik == 1, [As, Ad] = tdl38901(Q, bins, df, fc, models{im}, 100e-9);
            else, [As, Ad] = cdl38901(Q, bins, df, fc, models{im}, 100e-9); end
            H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1); Pc = mean(sum(abs(H(:, evalb)).^2, 1));
            for is = 1:2
                s2 = Pc / 10^(snrs(is)/10); Sall = qpsk(M, 8); Yall = H*Sall + noise(M, 8, s2); rp = mmse_rate(H, H, s2, evalb);
                for ip = 1:numel(Nps)
                    Hh = est_general(Yall(:, 1:Nps(ip)), Sall(:, 1:Nps(ip)), T, ones(Q,1), F, lamF, s2, 16, 20);
                    gap(t, im, ik, ip, is) = rp - mmse_rate(Hh, H, s2, evalb);
                end
            end
        end
    end
    fprintf('trial %2d: Np=8 20 dB gap TDL/CDL A %.2f/%.2f C %.2f/%.2f D %.2f/%.2f\n', t, squeeze(gap(t, :, :, 3, 1)).');
end
R = struct('ntrial', ntrial, 'snr_db', snrs, 'Np', Nps, 'models', {models}, 'gap_mean', squeeze(mean(gap, 1)), ...
           'dims', 'model x {TDL,CDL} x Np x snr');
for im = 1:3
    for ip = 1:numel(Nps)
        for is = 1:2
            [m, h] = tci(gap(:, im, 2, ip, is) - gap(:, im, 1, ip, is));
            R.(sprintf('cdl_minus_tdl_%s_Np%d_%ddB', models{im}, Nps(ip), snrs(is))) = [m h];
        end
    end
    fprintf('CDL-%s gap (Np 2/4/8) @20 %s | @30 %s | TDL @20 %s\n', models{im}, sprintf('%.3f ', squeeze(mean(gap(:, im, 2, :, 1), 1))), ...
        sprintf('%.3f ', squeeze(mean(gap(:, im, 2, :, 2), 1))), sprintf('%.3f ', squeeze(mean(gap(:, im, 1, :, 1), 1))));
end
R.gap_trials = gap;
% ---------------- D3 rank
Qs = [16 64 256 1024]; cases = {'narrow_far', 'wide_far', 'wide_far_squint', 'wide_near3m', 'wide_near10m'};
dfs = [30e3, 98.28e6/M, 98.28e6/M, 98.28e6/M, 98.28e6/M];
er99 = nan(numel(Qs), numel(cases)); er999 = er99; sv = cell(numel(Qs), numel(cases));
for iq = 1:numel(Qs)
    for ic = 1:numel(cases)
        e1 = zeros(1, 10); e2 = e1; s_acc = 0;
        for k = 1:10
            o = struct('squint', contains(cases{ic}, 'squint') || contains(cases{ic}, 'near'), 'nearfield_m', 0);
            if contains(cases{ic}, 'near3m'), o.nearfield_m = 3; elseif contains(cases{ic}, 'near10m'), o.nearfield_m = 10; end
            mdl = 'C'; if contains(cases{ic}, 'near'), mdl = 'D'; end
            [~, Ad] = cdl38901(Qs(iq), bins, dfs(ic), fc, mdl, 100e-9, o);
            s = svd(Ad.'); s2_ = s.^2 / sum(s.^2); cs = cumsum(s2_);
            e1(k) = find(cs >= 0.99, 1); e2(k) = find(cs >= 0.999, 1);
            s_acc = s_acc + s2_(1:min(64, numel(s2_))).';
        end
        er99(iq, ic) = mean(e1); er999(iq, ic) = mean(e2); sv{iq, ic} = s_acc/10;
        fprintf('Q = %4d %-16s effective rank 99%% %.1f | 99.9%% %.1f\n', Qs(iq), cases{ic}, er99(iq, ic), er999(iq, ic));
    end
end
R.rank = struct('Q', Qs, 'cases', {cases}, 'df_eff', dfs, 'eff_rank_99', er99, 'eff_rank_999', er999, 'sv_energy_first64', {sv});
fid = fopen(fullfile(out, 'week11_cdl_rank.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

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
