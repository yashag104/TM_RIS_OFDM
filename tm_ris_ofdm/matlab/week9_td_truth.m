%WEEK9_TD_TRUTH Review step 2, workstream A1 (pre-registered 2026-10-22).
%   A1a  Received pilots from the sampled-gating time-domain simulator (td_simulate, Nfft = 4096): per
%        element filter at absolute frequency, gate per slot, sum, FFT. No harmonic algebra. The true
%        matrix for scoring is obtained by probing the simulator with every unit vector.
%   A1b  Element memory beyond the quasi-static model. The cell is a single-pole system,
%        w' = -(xi_s + j 2 pi (fc - f0_s)) w + x,  r = -x + 2 xi_r,s w,  (Gamma_s = -1 + 2 xi_r,s G_s).
%        At a switch at time t_k the internal state w is continuous, so the solution is the new state's
%        steady state plus a transient (w_old(t_k) - w_new(t_k)) exp(a_new (t - t_k)). Its contribution to
%        received bin mbar is exact in closed form:
%          H_tr(mbar,m) = sum_q a^d_q(mbar) sum_{edges k} (2 xi_r,new / T) G_new(mbar)
%                         exp(-j 2 pi (mbar - m) df t_k) [G_old(m) - G_new(m)] a^s_q(m).
%   Estimator unchanged; Np = 2 pilots only; 30 channels, 20/30 dB; matched = closed form (eq. 2).

clear; rng(9800);
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; T = 1/df; K = 8; Nfft = 4096;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16; Np = 2; f = fc + bins(:)*df;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
Tn = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
      coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
xi = P.xr + P.xi; Gs = @(s, fr) 1 ./ (1j*2*pi*(fr - P.f0(s)) + xi(s));
snrs = [20 30]; ntrial = 30;
gap = nan(ntrial, 3, 2); rel = nan(ntrial, 2);           % [eq2 matched, td, eq2 + transients]
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H2 = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
    Htd = zeros(M);
    for m = 1:M, e = zeros(M, 1); e(m) = 1; Htd(:, m) = td_simulate(e, bins, Nfft, fc, df, seqs, As, Ad, P); end
    Htr = zeros(M);
    for q = 1:Q
        for k = 1:K
            kp = mod(k - 2, K) + 1; so = seqs(q, kp); sn = seqs(q, k);
            if so == sn, continue; end
            tk = (k - 1)*T/K;
            outv = Ad(q, :).' .* (2*P.xr(sn)/T) .* Gs(sn, f);
            inv_ = (Gs(so, f) - Gs(sn, f)) .* As(q, :).';
            Htr = Htr + exp(-1j*2*pi*(bins(:) - bins)*df*tk) .* (outv * inv_.');
        end
    end
    Hmem = H2 + Htr;
    rel(t, :) = [norm(Htd - H2, 'fro'), norm(Htr, 'fro')] / norm(H2, 'fro');
    Pc = mean(sum(abs(H2(:, evalb)).^2, 1));
    for is = 1:2
        s2 = Pc / 10^(snrs(is)/10); S = qpsk(M, Np); W = noise(M, Np, s2);
        Ytd = zeros(M, Np);
        for p = 1:Np, Ytd(:, p) = td_simulate(S(:, p), bins, Nfft, fc, df, seqs, As, Ad, P); end
        Hs = {H2, Htd, Hmem}; Ys = {H2*S + W, Ytd + W, Hmem*S + W};
        for c = 1:3
            gap(t, c, is) = mmse_rate(Hs{c}, Hs{c}, s2, evalb) - mmse_rate(est_general(Ys{c}, S, Tn, ones(Q,1), F, lamF, s2, 8, 20), Hs{c}, s2, evalb);
        end
    end
    fprintf('trial %2d: |Htd-H2|/|H2| %.1e |Htr|/|H2| %.1e | 30 dB gap eq2/td/mem %.3f/%.3f/%.3f\n', t, rel(t, :), squeeze(gap(t, :, 2)));
end
R = struct('ntrial', ntrial, 'snr_db', snrs, 'cases', {{'eq2_matched', 'timedomain', 'eq2_plus_transients'}}, ...
    'rel_td_vs_eq2', mean(rel(:, 1)), 'rel_transient', mean(rel(:, 2)), 'gap_mean', squeeze(mean(gap, 1)));
for is = 1:2
    [m1, h1] = tci(gap(:, 2, is) - gap(:, 1, is)); [m2, h2] = tci(gap(:, 3, is) - gap(:, 1, is));
    R.(sprintf('td_minus_matched_%ddB', snrs(is))) = [m1 h1]; R.(sprintf('mem_minus_matched_%ddB', snrs(is))) = [m2 h2];
    fprintf('%d dB: td - matched %+.3f +-%.3f | transients - matched %+.3f +-%.3f\n', snrs(is), m1, h1, m2, h2);
end
fprintf('relative size: td vs eq2 %.1e, transient term %.1e\n', R.rel_td_vs_eq2, R.rel_transient);
R.gap_trials = gap; R.rel_trials = rel;
fid = fopen(fullfile(out, 'week9_td_truth.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

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
