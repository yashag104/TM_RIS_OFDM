%WEEK9_SYNC_WINDOW Step 1b, test F1 (pre-registered 2026-10-22, fresh seeds).
%   Mode (ii) timing (RIS locked to the BS clock) plus a receiver sync error delta_s ~ U(-150,150) ns.
%   Evidence search with (a) default window [0, 0.3 us], (b) margin window [-150 ns, 0.3 us],
%   (c) window start t0 in {-150,-75,0,+75} ns chosen jointly with d by the evidence.
%   Reference: perfect sync, known offset, default window. 30 channels, 20/30 dB, Np = 2, L = 16.
clear; rng(9500);
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; Tu = 1/df; K = 8; Tslot = Tu/K;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16;
base = [ones(1, 6), 2*ones(1, 2)]; st = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
E0 = (0:K)/K;
Tnom = {coupling_edges(bins, fc, df, E0, base, ones(1,M), ones(1,M), 'state1', M-1), ...
        coupling_edges(bins, fc, df, E0, base, ones(1,M), ones(1,M), 'state2', M-1)};
ramp = @(tau) exp(-1j*2*pi*bins*df*tau);

[Fd, Ld] = dpss_basis(bins, df, 0.3e-6, 16);
[Fb, Lb] = dpss_basis(bins, df, 0.45e-6, 16); Fb = exp(-1j*2*pi*bins(:)*df*(-150e-9)) .* Fb;
t0s = [-150 -75 0 75]*1e-9; Fc = cell(1, 4);
for k = 1:4, Fc{k} = exp(-1j*2*pi*bins(:)*df*t0s(k)) .* Fd; end
snrs = [20 30]; ntrial = 30; Np = 2; grid = -0.5:0.025:0.5;
gap = nan(ntrial, 4, 2); t0pick = nan(ntrial, 2); dsy = nan(ntrial, 1);
for t = 1:ntrial
    [As0, Ad0] = phys_channel(Q, bins, df, fc, popt);
    tb = 50e-9 + 250e-9*rand; tu = 30e-9 + 270e-9*rand; dsy(t) = (2*rand - 1)*150e-9;
    Hs = conj(ramp(tb + tu)).' .* coupling_edges(bins, fc, df, repmat(E0, Q, 1), st, As0 .* ramp(tb), Ad0 .* ramp(tu), P, M-1);
    He = conj(ramp(dsy(t))).' .* Hs;                     % receiver window late by delta_s
    dtr = -tb/Tslot; Pc = mean(sum(abs(Hs(:, evalb)).^2, 1));
    for is = 1:2
        s2 = Pc / 10^(snrs(is)/10); S = qpsk(M, Np); W = noise(M, Np, s2);
        Ys = Hs*S + W; Ye = He*S + W; rps = mmse_rate(Hs, Hs, s2, evalb); rpe = mmse_rate(He, He, s2, evalb);
        gap(t, 1, is) = rps - mmse_rate(est_general(Ys, S, shiftk(Tnom, bins, dtr, K), ones(Q,1), Fd, Ld, s2, 8, 20), Hs, s2, evalb);
        [d1, ~] = esearch(Ye, S, Tnom, bins, K, Fd, Ld, s2, grid);
        gap(t, 2, is) = rpe - mmse_rate(est_general(Ye, S, shiftk(Tnom, bins, d1, K), ones(Q,1), Fd, Ld, s2, 8, 20), He, s2, evalb);
        [d2, ~] = esearch(Ye, S, Tnom, bins, K, Fb, Lb, s2, grid);
        gap(t, 3, is) = rpe - mmse_rate(est_general(Ye, S, shiftk(Tnom, bins, d2, K), ones(Q,1), Fb, Lb, s2, 8, 20), He, s2, evalb);
        dk = zeros(1, 4); ek = zeros(1, 4);
        for k = 1:4, [dk(k), ek(k)] = esearch(Ye, S, Tnom, bins, K, Fc{k}, Ld, s2, grid); end
        [~, kb] = max(ek); t0pick(t, is) = t0s(kb);
        gap(t, 4, is) = rpe - mmse_rate(est_general(Ye, S, shiftk(Tnom, bins, dk(kb), K), ones(Q,1), Fc{kb}, Ld, s2, 8, 20), He, s2, evalb);
    end
    fprintf('trial %2d: delta_s %+4.0f ns | 30 dB gaps ref/a/b/c %.3f/%.3f/%.3f/%.3f (t0 %+4.0f ns)\n', t, dsy(t)*1e9, squeeze(gap(t, :, 2)), t0pick(t, 2)*1e9);
end
names = {'reference', 'a_default', 'b_margin', 'c_evidence_t0'};
R = struct('ntrial', ntrial, 'snr_db', snrs, 'names', {names}, 'gap_mean', squeeze(mean(gap, 1)));
for is = 1:2
    for v = 2:4
        [m, h] = tci(gap(:, v, is) - gap(:, 1, is));
        R.(sprintf('%s_%ddB', names{v}, snrs(is))) = [m h m+h];
        fprintf('%d dB %-14s loss vs reference %+.3f +-%.3f (upper %.3f)\n', snrs(is), names{v}, m, h, m+h);
    end
end
R.gap_trials = gap; R.t0pick = t0pick; R.sync_err = dsy;
fid = fopen(fullfile(out, 'week9_sync_window.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function [d, e] = esearch(Y, S, Tnom, bins, K, F, L, s2, grid)
fe = @(d) -evidence(d, Y, S, Tnom, bins, K, F, L, s2);
v = arrayfun(fe, grid); [~, i0] = min(v);
d = golden(fe, grid(max(i0-1,1)), grid(min(i0+1,numel(grid))), 14); e = -fe(d);
end
%% ---------------------------------------------------------------------
function [m, h] = tci(x)
% mean and half-width of the t-based 95% CI (base MATLAB: t quantile via betaincinv)
n = numel(x); nu = n - 1; tq = sqrt(nu*(1/betaincinv(0.05, nu/2, 0.5) - 1));
m = mean(x); h = tq*std(x)/sqrt(n);
end
function le = evidence(d, Y, S, Tnom, bins, K, F, lam, s2)
T = shiftk(Tnom, bins, d, K); [M, Np] = size(S); L = size(F, 2); y = Y(:);
Phi = zeros(M*Np, 2*L*L);
for s = 1:2
    for p = 1:Np
        Gs = T{s} * (S(:,p) .* F);
        for i = 1:L, Phi((p-1)*M+(1:M), (s-1)*L*L+(i-1)*L+(1:L)) = F(:, i) .* Gs; end
    end
end
pl = kron(lam(:), lam(:)); w = [pl; pl];
sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
Dw = sqrt(sc*w); B = (Dw .* (Phi'*Phi)) .* Dw.' / s2 + eye(numel(w));
Rc = chol((B + B')/2); z = Rc' \ (Dw .* (Phi'*y)) / s2;
le = -(real(y'*y)/s2 - real(z'*z)) - 2*sum(log(real(diag(Rc))));
end
function T = shiftk(Tnom, bins, d, K)
D = exp(-1j*2*pi*bins(:)*d/K); T = {D .* Tnom{1} .* D', D .* Tnom{2} .* D'};
end
function x = golden(f, a, b, n)
g = (sqrt(5) - 1)/2; c = b - g*(b - a); d = a + g*(b - a); fc = f(c); fd = f(d);
for i = 1:n
    if fc < fd, b = d; d = c; fd = fc; c = b - g*(b - a); fc = f(c);
    else, a = c; c = d; fc = fd; d = a + g*(b - a); fd = f(d); end
end
x = (a + b)/2;
end
function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
function W = noise(M, n, s2), W = sqrt(s2/2) * (randn(M,n) + 1j*randn(M,n)); end
function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
