%WEEK11_MOTIVATION Review I1-I3 (pre-registered 2026-10-22).
%   Physically consistent normalisation: fixed transmit power and noise. sigma^2 is set once per channel so
%   that the optimised static RIS has the reference SNR (10/20/30 dB); every configuration of that channel
%   uses the same sigma^2, so received power differences count.
%   I1  legitimate user: no RIS (direct path only, at -10 dB of the static-RIS path; or blocked = 0);
%       static 1-bit RIS with states co-phased at the user; static RIS with random states;
%       TM-RIS designed for the user (per-element polarity co-phasing h = 0 + random cyclic shifts, as
%       week3_designed_user) with one-tap, with ours (timing-only, Np = 2 and 4) and with perfect CSI.
%   I2  eavesdropper at another position (own RIS -> eve hop, same BS -> RIS hop), same pilots:
%       (a) without the pattern: one-tap, CE-BEM |h| <= 40 (robust equaliser); (b) with the pattern: ours.
%       Secrecy rate [R_legit - R_eve]^+ for each eavesdropper.
%   I3  designed-user gain across duty D in {0.625, 0.75, 0.875} (K = 8) and K in {4, 16} (D = 0.75).
%   Shared-path channels, 30 channels per point.

clear; rng(11300);
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; c0 = M/2 + 1; Q = 16; f = fc + bins(:)*df;
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
L2 = ceil(M*df*0.6e-6) + 12; [F2, lam2] = dpss_basis(bins, df, 0.6e-6, L2);
G0 = [element_response(f, 1, P), element_response(f, 2, P)];
cfgs = struct('D', {0.625, 0.75, 0.875, 0.75, 0.75}, 'K', {8, 8, 8, 4, 16});
snrs = [10 20 30]; ntrial = 30;
names = {'noRIS_direct', 'static_opt', 'static_random', 'static_opt_direct', 'TM_onetap', 'TM_ours_Np2', 'TM_ours_Np4', 'TM_perfect', 'TM_direct_ours_Np4', ...
         'eve_onetap', 'eve_bem40', 'eve_ours_withpattern'};
R = struct('snr_db', snrs, 'ntrial', ntrial, 'names', {names});
for ic = 1:numel(cfgs)
    D = cfgs(ic).D; K = cfgs(ic).K;
    base = [ones(1, round(D*K)), 2*ones(1, K - round(D*K))];
    rate = nan(ntrial, numel(names), numel(snrs));
    for t = 1:ntrial
        [As, Ad] = phys_channel(Q, bins, df, fc, popt); [~, Ae] = phys_channel(Q, bins, df, fc, popt);
        [~, adir] = phys_channel(1, bins, df, fc, popt);            % direct BS -> user path (1 x M)
        % static RIS: co-phased states at the band centre (greedy over a common reference)
        cq = Ad(:, c0) .* As(:, c0); best = -inf;
        for psi = linspace(0, 2*pi, 33)
            sq = 1 + (real(cq .* G0(c0, 1) * exp(-1j*psi)) < real(cq .* G0(c0, 2) * exp(-1j*psi)));
            v = abs(sum(cq .* G0(c0, sq).'));
            if v > best, best = v; sopt = sq; end
        end
        hst = sum(Ad.' .* G0(:, sopt) .* As.', 2);                % static diagonal channel
        hrd = sum(Ad.' .* G0(:, randi(2, Q, 1)) .* As.', 2);
        Pst = mean(abs(hst(evalb)).^2);
        adir = adir.' * sqrt(Pst * 0.1 / mean(abs(adir(evalb)).^2));   % direct path at -10 dB of static RIS
        % TM-RIS designed for the user (as week3_designed_user)
        shifts = randi(K, Q, 1) - 1; b0 = harmonic_coeffs(fc, base, 0, P); best = -inf;
        for psi = linspace(0, 2*pi, 33)
            p = real(cq * b0 * exp(-1j*psi)) >= 0; v = abs(sum(cq .* b0 .* (2*p - 1)));
            if v > best, best = v; pol = p; end
        end
        seqs = zeros(Q, K);
        for q = 1:Q, s = circshift(base, shifts(q)); if ~pol(q), s = 3 - s; end, seqs(q, :) = s; end
        H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
        He = coupling_matrix(bins, fc, df, seqs, As, Ae, P, M-1);
        [useq, ~, grp] = unique(seqs, 'rows'); Tg = cell(size(useq, 1), 2);
        for g = 1:size(useq, 1)
            Tg{g,1} = coupling_matrix(bins, fc, df, useq(g, :), ones(1,M), ones(1,M), 'state1', M-1);
            Tg{g,2} = coupling_matrix(bins, fc, df, useq(g, :), ones(1,M), ones(1,M), 'state2', M-1);
        end
        Hd = H + diag(adir);
        for is = 1:numel(snrs)
            s2 = Pst / 10^(snrs(is)/10);
            rate(t, 1, is) = mean(log2(1 + abs(adir(evalb)).^2 / s2));
            rate(t, 2, is) = mean(log2(1 + abs(hst(evalb)).^2 / s2));
            rate(t, 3, is) = mean(log2(1 + abs(hrd(evalb)).^2 / s2));
            rate(t, 4, is) = mean(log2(1 + abs(hst(evalb) + adir(evalb)).^2 / s2));
            dg = abs(diag(H)).^2; rate(t, 5, is) = mean(log2(1 + dg(evalb) ./ (sum(abs(H(evalb, :)).^2, 2) - dg(evalb) + s2)));
            S = qpsk(M, 4); W = noise(M, 4, s2);
            Y = H*S + W; Ye = He*S + noise(M, 4, s2); Yd = Hd*S + W;
            rate(t, 6, is) = mmse_rate(est_general(Y(:, 1:2), S(:, 1:2), Tg, grp, F, lamF, s2, 8, 20), H, s2, evalb);
            rate(t, 7, is) = mmse_rate(est_general(Y, S, Tg, grp, F, lamF, s2, 8, 20), H, s2, evalb);
            rate(t, 8, is) = mmse_rate(H, H, s2, evalb);
            rate(t, 9, is) = mmse_rate(est_general(Yd, S, Tg, grp, F, lamF, s2, 8, 20), Hd, s2, evalb);
            dge = abs(diag(He)).^2; rate(t, 10, is) = mean(log2(1 + dge(evalb) ./ (sum(abs(He(evalb, :)).^2, 2) - dge(evalb) + s2)));
            Hb = est_bem_cg(Ye, S, F2, lam2, s2, -40:40); se = max(norm(Ye - Hb*S, 'fro')^2/(M*4) - s2, 0);
            rate(t, 11, is) = mmse_rate_se(Hb, He, s2, evalb, se);
            rate(t, 12, is) = mmse_rate(est_general(Ye, S, Tg, grp, F, lamF, s2, 8, 20), He, s2, evalb);
        end
    end
    key = sprintf('D%03d_K%d', round(1000*D), K);
    m = squeeze(mean(rate, 1));
    C = struct('rate_mean', m, 'dims', 'config x snr', 'n_distinct_patterns', size(useq, 1));
    for is = 1:numel(snrs)
        [g1, h1] = tci(rate(:, 7, is) - rate(:, 5, is)); [g2, h2] = tci(rate(:, 7, is) - rate(:, 2, is));
        C.(sprintf('ours_minus_onetap_%ddB', snrs(is))) = [g1 h1]; C.(sprintf('ours_minus_staticopt_%ddB', snrs(is))) = [g2 h2];
        sec = max(rate(:, 7, is) - rate(:, 10:12, is), 0);
        C.(sprintf('secrecy_vs_eve_%ddB', snrs(is))) = mean(sec, 1);
    end
    C.rate_trials = rate; R.(key) = C;
    fprintf('%s (G = %d): 20 dB static_opt %.2f | TM one-tap %.2f ours4 %.2f perfect %.2f | eve one-tap %.2f bem %.2f ours %.2f\n', key, size(useq, 1), ...
        m(2, 2), m(5, 2), m(7, 2), m(8, 2), m(10, 2), m(11, 2), m(12, 2));
    fid = fopen(fullfile(out, 'week11_motivation.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);
end

function r = mmse_rate_se(Hh, H, s2, evalb, se2)
M = size(H, 1);
W = (Hh'*Hh + (s2 + se2)*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2)); r = mean(log2(1 + sinr(evalb)));
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
