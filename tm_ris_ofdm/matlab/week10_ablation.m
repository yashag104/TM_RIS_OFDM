%WEEK10_ABLATION Review step 4: C2 ablation + G3 ALS diagnostics (pre-registered 2026-10-22).
%   (a) linear first stage only, (b) ALS rank r in {2,4,6,8,16} without fallback, (c) full pipeline
%   (ALS r0 with fallback; r0 = 8 shared-path, 16 TDL), (d) full + one DD pass (6 data symbols,
%   16-QAM at 20 dB, 64-QAM at 30 dB). G3: fallback fraction, residual per iteration, 5 random restarts at r = 8.
%   Shared-path 30 ns and TDL-C 100 ns, Np in {2,4,8}, 30 channels, 20/30 dB.

clear; rng(10200);
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; K = 8;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16; Nd = 6;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
chans = {'phys30', 'TDLC100'}; tw = [0.3e-6 1e-6]; Ls = [16, ceil(M*df*1e-6) + 10]; r0 = [8 16];
ranks = [2 4 6 8 16]; Nps = [2 4 8]; snrs = [20 30]; qam = [16 64]; ntrial = 30; nrs = 5;
vars = [{'linear'}, arrayfun(@(r) sprintf('als_r%d', r), ranks, 'uni', 0), {'full', 'full_dd'}];
R = struct('Np', Nps, 'snr_db', snrs, 'ntrial', ntrial, 'variants', {vars});
for ic = 1:2
    [F, lamF] = dpss_basis(bins, df, tw(ic), Ls(ic));
    gap = nan(ntrial, numel(vars), numel(Nps), 2); fb = nan(ntrial, numel(Nps), 2);
    rs_gap = nan(ntrial, nrs, 2); resit = nan(ntrial, 20, 2);
    for t = 1:ntrial
        if ic == 1, [As, Ad] = phys_channel(Q, bins, df, fc, popt);
        else, [As, Ad] = tdl38901(Q, bins, df, fc, 'C', 100e-9); end
        H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
        Pc = mean(sum(abs(H(:, evalb)).^2, 1));
        for is = 1:2
            s2 = Pc / 10^(snrs(is)/10); rp = mmse_rate(H, H, s2, evalb); Mq = qam(is);
            Sall = qpsk(M, max(Nps)); Yall = H*Sall + noise(M, max(Nps), s2);
            Xd = qammod_(randi(Mq, M, Nd) - 1, Mq); Yd = H*Xd + noise(M, Nd, s2);
            for ip = 1:numel(Nps)
                S = Sall(:, 1:Nps(ip)); Y = Yall(:, 1:Nps(ip)); k = 1;
                [Hf, inf0] = est_general_ab(Y, S, T, ones(Q,1), F, lamF, s2, r0(ic), 20, struct('fallback', true));
                gap(t, 1, ip, is) = rp - mmse_rate(inf0.H0, H, s2, evalb);
                for ir = 1:numel(ranks)
                    [~, inf1] = est_general_ab(Y, S, T, ones(Q,1), F, lamF, s2, ranks(ir), 20, struct('fallback', false));
                    gap(t, 1 + ir, ip, is) = rp - mmse_rate(inf1.Hals, H, s2, evalb);
                    if ranks(ir) == 8 && ip == 1, resit(t, :, is) = inf1.res_iter; end
                end
                gap(t, end-1, ip, is) = rp - mmse_rate(Hf, H, s2, evalb); fb(t, ip, is) = inf0.used_fallback;
                Xh = qamdec_(((Hf'*Hf + s2*eye(M)) \ Hf') * Yd, Mq);
                Hdd = est_general_ab([Y Yd], [S Xh], T, ones(Q,1), F, lamF, s2, r0(ic), 20, struct('fallback', true));
                gap(t, end, ip, is) = rp - mmse_rate(Hdd, H, s2, evalb);
            end
            for k = 1:nrs                                       % random restarts, r = 8, Np = 2
                [~, inf2] = est_general_ab(Sall(:, 1:2)*0 + Yall(:, 1:2), Sall(:, 1:2), T, ones(Q,1), F, lamF, s2, 8, 20, ...
                                          struct('fallback', false, 'init', 'random'));
                rs_gap(t, k, is) = rp - mmse_rate(inf2.Hals, H, s2, evalb);
            end
        end
        fprintf('%s trial %2d: Np=2 30 dB gap lin/r2/r4/r6/r8/r16/full/dd %s\n', chans{ic}, t, sprintf('%.2f ', squeeze(gap(t, :, 1, 2))));
    end
    C = struct('gap_mean', squeeze(mean(gap, 1)), 'dims', 'variant x Np x snr', 'fallback_frac', squeeze(mean(fb, 1)), ...
               'restart_gap_mean', squeeze(mean(rs_gap, 1)), 'restart_gap_spread', squeeze(mean(max(rs_gap, [], 2) - min(rs_gap, [], 2), 1)), ...
               'res_iter_mean', squeeze(mean(resit, 1)));
    for ip = 1:numel(Nps)
        for is = 1:2
            [m1, h1] = tci(gap(:, 1, ip, is) - gap(:, end-1, ip, is));   % linear - full
            [m2, h2] = tci(gap(:, end, ip, is) - gap(:, end-1, ip, is)); % dd - full
            C.(sprintf('linear_minus_full_Np%d_%ddB', Nps(ip), snrs(is))) = [m1 h1];
            C.(sprintf('dd_minus_full_Np%d_%ddB', Nps(ip), snrs(is))) = [m2 h2];
            fprintf('%s Np=%d %d dB: gaps %s | lin-full %+.3f +-%.3f | dd-full %+.3f +-%.3f | fallback %.0f%%\n', chans{ic}, Nps(ip), snrs(is), ...
                sprintf('%.3f ', squeeze(mean(gap(:, :, ip, is), 1))), m1, h1, m2, h2, 100*mean(fb(:, ip, is)));
        end
    end
    C.gap_trials = gap; C.restart_trials = rs_gap;
    R.(chans{ic}) = C;
    fid = fopen(fullfile(out, 'week10_ablation.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);
end

function [m, h] = tci(x)
n = numel(x); nu = n - 1; tq = sqrt(nu*(1/betaincinv(0.05, nu/2, 0.5) - 1));
m = mean(x); h = tq*std(x)/sqrt(n);
end
function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
function W = noise(M, n, s2), W = sqrt(s2/2) * (randn(M,n) + 1j*randn(M,n)); end
function x = qammod_(d, Mq)
m = sqrt(Mq); lv = -(m-1):2:(m-1);
x = reshape((lv(mod(d, m) + 1) + 1j*lv(floor(d / m) + 1)) / sqrt(2*(Mq-1)/3), size(d));
end
function x = qamdec_(z, Mq)
m = sqrt(Mq); sc = sqrt(2*(Mq-1)/3); lv = -(m-1):2:(m-1);
q = @(v) lv(min(max(round((v + m - 1)/2), 0), m-1) + 1);
x = reshape((q(real(z)*sc) + 1j*q(imag(z)*sc)) / sc, size(z));
end
function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
