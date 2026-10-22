%WEEK10_RISETIME Review step 4: C5 finite switching rise time (pre-registered 2026-10-22).
%   Linear ramp of duration rho (fraction of a slot) between states, in the truth and in our calibrated kernels.
%   Reports the harmonic energy outside |h| <= 12, banded MMSE (+-4, +-16, +-32) with the TRUE channel,
%   CE-BEM |h| <= 12 and 40 with Np = 8 (naive and oracle-robust equaliser), ours with Np = 2.
%   Shared-path channel, 30 channels, 20/30 dB.
clear; rng(10400);
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; K = 8; Nt = 32768;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16;
base = [ones(1, 6), 2*ones(1, 2)];
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
L2 = ceil(M*df*0.6e-6) + 12; [F2, lam2] = dpss_basis(bins, df, 0.6e-6, L2);
rhos = [0 0.01 0.05 0.1 0.3]; snrs = [20 30]; ntrial = 30;
names = {'banded4_trueH', 'banded16_trueH', 'banded32_trueH', 'bem12_Np8', 'bem40_Np8', 'bem12_Np8_robust', 'bem40_Np8_robust', 'ours_Np2', 'perfect'};
tail = zeros(1, numel(rhos)); rate = nan(ntrial, numel(names), numel(rhos), 2);
for ir = 1:numel(rhos)
    W1 = ramp_weights(base, 1, Nt, rhos(ir)); w = squeeze(W1(1, :, 1)) - squeeze(W1(1, :, 2));
    X = fft(w)/Nt; hh = [0:Nt/2-1, -Nt/2:-1]; E = abs(X).^2; tail(ir) = sum(E(abs(hh) > 12))/sum(E);
end
fprintf('energy outside |h|<=12 (ideal-element pattern): %s\n', sprintf('%.2e ', tail));
for ir = 1:numel(rhos)
    Tk = ramp_kernels(bins, fc, df, base, rhos(ir), M, Nt);
    Wq = ramp_weights(base, Q, Nt, rhos(ir));
    for t = 1:ntrial
        [As, Ad] = phys_channel(Q, bins, df, fc, popt);
        H = coupling_waveform(bins, fc, df, Wq, As, Ad, P, M-1);
        Pc = mean(sum(abs(H(:, evalb)).^2, 1));
        for is = 1:2
            s2 = Pc / 10^(snrs(is)/10);
            S = qpsk(M, 8); Y = H*S + noise(M, 8, s2);
            for b = 1:3, rate(t, b, ir, is) = eq_banded([4 16 32] * (b == 1:3).', H, H, s2, evalb); end
            Hb = {est_bem_cg(Y, S, F2, lam2, s2, -12:12), est_bem_cg(Y, S, F2, lam2, s2, -40:40)};
            for b = 1:2
                rate(t, 3 + b, ir, is) = mmse_rate(Hb{b}, H, s2, evalb);
                rate(t, 5 + b, ir, is) = mmse_rate_se(Hb{b}, H, s2, evalb, norm(Hb{b} - H, 'fro')^2/M);
            end
            rate(t, 8, ir, is) = mmse_rate(est_general(Y(:, 1:2), S(:, 1:2), Tk, ones(Q,1), F, lamF, s2, 8, 20), H, s2, evalb);
            rate(t, 9, ir, is) = mmse_rate(H, H, s2, evalb);
        end
    end
    fprintf('rho %.2f: 30 dB rates %s\n', rhos(ir), sprintf('%.2f ', squeeze(mean(rate(:, :, ir, 2), 1))));
end
R = struct('rho', rhos, 'snr_db', snrs, 'ntrial', ntrial, 'names', {names}, 'energy_outside_h12', tail, ...
    'rate_mean', squeeze(mean(rate, 1)), 'dims', 'method x rho x snr');
R.rate_trials = rate;
fid = fopen(fullfile(out, 'week10_risetime.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function r = eq_banded(Wd, He, H, s2, evalb)
M = size(H, 1); W = zeros(M);
for m = 1:M
    rows = max(1, m-Wd):min(M, m+Wd); cols = max(1, m-2*Wd):min(M, m+2*Wd);
    Hw = He(rows, cols); k = find(cols == m);
    w = (Hw*Hw' + s2*eye(numel(rows))) \ Hw(:, k); W(m, rows) = w';
end
G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2)); r = mean(log2(1 + sinr(evalb)));
end
function r = mmse_rate_se(Hh, H, s2, evalb, se2)
M = size(H, 1);
W = (Hh'*Hh + (s2 + se2)*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2)); r = mean(log2(1 + sinr(evalb)));
end
function W = ramp_weights(base, Q, Nt, rho)
K = numel(base); t = ((0:Nt-1) + 0.5) / Nt * K;
w1 = double(base(floor(t) + 1) == 1);
if rho > 0
    for k = 0:K-1
        a = base(mod(k-1, K) + 1); b = base(k + 1);
        if a == b, continue; end
        d = t - k; d(d > K/2) = d(d > K/2) - K; d(d < -K/2) = d(d < -K/2) + K;
        in = abs(d) < rho/2; ramp = (d(in) + rho/2) / rho;
        if b == 1, w1(in) = ramp; else, w1(in) = 1 - ramp; end
    end
end
W = zeros(Q, Nt, 2); W(:, :, 1) = repmat(w1, Q, 1); W(:, :, 2) = repmat(1 - w1, Q, 1);
end

function T = ramp_kernels(bins, fc, df, base, rho, M, Nt)
% Timing kernels from the (possibly ramped) switching waveform, state indicators.
W = ramp_weights(base, 1, Nt, rho);
T = {coupling_waveform(bins, fc, df, W, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_waveform(bins, fc, df, W, ones(1,M), ones(1,M), 'state2', M-1)};
end


function [m, h] = tci(x)
x = x(~isnan(x)); n = numel(x); nu = n - 1; tq = sqrt(nu*(1/betaincinv(0.05, nu/2, 0.5) - 1));
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
