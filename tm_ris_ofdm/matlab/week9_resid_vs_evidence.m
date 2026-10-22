%WEEK9_RESID_VS_EVIDENCE Step 1b, test F2 (pre-registered 2026-10-22, fresh seeds).
%   Mode (ii) physical timing (offset = BS-RIS delay), 30 channels, 20/30 dB. Paired gap difference
%   residual - evidence (t-CI) in regimes R1 (week-7 coarse residual procedure, Np = 2), R2 (fine-grid
%   residual, Np = 1), R3 (fine-grid residual, TDL-C 100 ns, window 1 us, r = 16, Np = 8), and R4 wall-clock
%   per search (fine residual vs evidence).
clear; rng(9600);
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; Tu = 1/df; K = 8; Tslot = Tu/K;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16;
base = [ones(1, 6), 2*ones(1, 2)]; st = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
E0 = (0:K)/K;
Tnom = {coupling_edges(bins, fc, df, E0, base, ones(1,M), ones(1,M), 'state1', M-1), ...
        coupling_edges(bins, fc, df, E0, base, ones(1,M), ones(1,M), 'state2', M-1)};
ramp = @(tau) exp(-1j*2*pi*bins*df*tau);

[F0, L0] = dpss_basis(bins, df, 0.3e-6, 16);
[Ft, Lt] = dpss_basis(bins, df, 1e-6, ceil(M*df*1e-6) + 10);
snrs = [20 30]; ntrial = 30; gc = -0.5:0.025:0.5; gf = -0.5:0.01:0.5;
reg = {'R1_coarse_Np2', 'R2_fine_Np1', 'R3_fine_TDLC100_Np8'};
gap = nan(ntrial, numel(reg), 3, 2); tm = nan(ntrial, numel(reg), 2);   % [resid, evidence, known]
for t = 1:ntrial
    tb = 50e-9 + 250e-9*rand; tu = 30e-9 + 270e-9*rand; dtr = -tb/Tslot;
    [As0, Ad0] = phys_channel(Q, bins, df, fc, popt);
    Hp = conj(ramp(tb + tu)).' .* coupling_edges(bins, fc, df, repmat(E0, Q, 1), st, As0 .* ramp(tb), Ad0 .* ramp(tu), P, M-1);
    [As1, Ad1] = tdl38901(Q, bins, df, fc, 'C', 100e-9);
    Ht = conj(ramp(tb + tu)).' .* coupling_edges(bins, fc, df, repmat(E0, Q, 1), st, As1 .* ramp(tb), Ad1 .* ramp(tu), P, M-1);
    for ir = 1:3
        if ir < 3, H = Hp; F = F0; L = L0; r = 8; else, H = Ht; F = Ft; L = Lt; r = 16; end
        Np = [2 1 8]; Np = Np(ir);
        Pc = mean(sum(abs(H(:, evalb)).^2, 1));
        for is = 1:2
            s2 = Pc / 10^(snrs(is)/10); S = qpsk(M, Np); Y = H*S + noise(M, Np, s2); rp = mmse_rate(H, H, s2, evalb);
            if ir == 1
                fr = @(d) resfun(d, Y, S, Tnom, bins, K, F, L, s2, Q, r, 8); g = gc;
            else
                fr = @(d) resfun(d, Y, S, Tnom, bins, K, F, L, s2, Q, r, 20); g = gf;
            end
            t0 = tic; v = arrayfun(fr, g); [~, i0] = min(v);
            dr = golden(fr, g(max(i0-1,1)), g(min(i0+1,numel(g))), 14); tr = toc(t0);
            fe = @(d) -evidence(d, Y, S, Tnom, bins, K, F, L, s2);
            t0 = tic; v = arrayfun(fe, gc); [~, i0] = min(v);
            de = golden(fe, gc(max(i0-1,1)), gc(min(i0+1,numel(gc))), 14); te = toc(t0);
            ds = [dr de dtr];
            for c = 1:3
                gap(t, ir, c, is) = rp - mmse_rate(est_general(Y, S, shiftk(Tnom, bins, ds(c), K), ones(Q,1), F, L, s2, r, 20), H, s2, evalb);
            end
            if is == 2, tm(t, ir, :) = [tr te]; end
        end
    end
    fprintf('trial %2d done | 30 dB resid-evid R1/R2/R3 %+.3f %+.3f %+.3f\n', t, squeeze(gap(t, :, 1, 2) - gap(t, :, 2, 2)));
end
R = struct('ntrial', ntrial, 'snr_db', snrs, 'regimes', {reg}, 'gap_mean_resid_evid_known', squeeze(mean(gap, 1)), ...
    'time_resid_evid_s', squeeze(mean(tm, 1)));
for ir = 1:3
    for is = 1:2
        [m, h] = tci(gap(:, ir, 1, is) - gap(:, ir, 2, is)); [m2, h2] = tci(gap(:, ir, 2, is) - gap(:, ir, 3, is));
        R.(sprintf('%s_%ddB', reg{ir}, snrs(is))) = struct('resid_minus_evid', [m h], 'evid_minus_known', [m2 h2]);
        fprintf('%-20s %d dB: resid-evid %+.3f +-%.3f | evid-known %+.3f +-%.3f | time resid %.1f s evid %.1f s\n', reg{ir}, snrs(is), m, h, m2, h2, mean(tm(:, ir, 1)), mean(tm(:, ir, 2)));
    end
end
R.gap_trials = gap; R.time_trials = tm;
fid = fopen(fullfile(out, 'week9_resid_vs_evidence.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function v = resfun(d, Y, S, Tnom, bins, K, F, L, s2, Q, r, it)
Hc = est_general(Y, S, shiftk(Tnom, bins, d, K), ones(Q,1), F, L, s2, r, it);
v = norm(Y - Hc*S, 'fro')^2;
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
