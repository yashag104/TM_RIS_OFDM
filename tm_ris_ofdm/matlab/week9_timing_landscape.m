%WEEK9_TIMING_LANDSCAPE Review step 1, test T3 (workstreams B2, B4). Pre-registered 2026-10-22.
%   Free off-grid offset d0 ~ U(-0.5, 0.5) slot (stress test, as week8_timing_evidence), shared-path
%   channel, 30 channels, 20/30 dB. Four pilot symbols are sent: 2 are used for estimation, 2 are held
%   out. For d on a grid over +-0.5 slot: in-sample pilot residual, held-out residual (estimate from the
%   2 in-sample pilots applied to the 2 held-out pilots), and evidence (2 in-sample pilots).
%   Then each criterion is refined by golden section. Hypothesis B4: wrong d fits the in-sample pilots
%   (flat residual) but predicts held-out pilots badly; the evidence tracks the held-out criterion.
%   PASS: held-out minimiser within 0.005 slot of the evidence maximiser in >= 90% of channels at 30 dB.
%   Also: projection error of shifted hop responses onto the DPSS basis (L = 16) vs delay 0..2 us.

clear; rng(9300);
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; Tu = 1/df; K = 8; Tslot = Tu/K;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16;
base = [ones(1, 6), 2*ones(1, 2)]; st = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
E0 = (0:K)/K;
Tnom = {coupling_edges(bins, fc, df, E0, base, ones(1,M), ones(1,M), 'state1', M-1), ...
        coupling_edges(bins, fc, df, E0, base, ones(1,M), ones(1,M), 'state2', M-1)};
snrs = [20 30]; ntrial = 30; dg = -0.5:0.01:0.5;
res_in = nan(ntrial, 2, numel(dg)); res_out = res_in; evid = res_in;
best = nan(ntrial, 2, 4); lossd = nan(ntrial, 2, 4);      % [resid_in, heldout, evidence, true]
for t = 1:ntrial
    d0 = rand - 0.5;
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H = coupling_edges(bins, fc, df, repmat(((0:K) + d0)/K, Q, 1), st, As, Ad, P, M-1);
    Pc = mean(sum(abs(H(:, evalb)).^2, 1));
    for is = 1:2
        s2 = Pc / 10^(snrs(is)/10);
        S = qpsk(M, 4); Y = H*S + noise(M, 4, s2);
        Si = S(:, 1:2); Yi = Y(:, 1:2); So = S(:, 3:4); Yo = Y(:, 3:4);
        fr = @(d) resid(d, Yi, Si, Yi, Si, Tnom, bins, K, F, lamF, s2, Q);
        fo = @(d) resid(d, Yi, Si, Yo, So, Tnom, bins, K, F, lamF, s2, Q);
        fe = @(d) -evidence(d, Yi, Si, Tnom, bins, K, F, lamF, s2);
        for k = 1:numel(dg), res_in(t, is, k) = fr(dg(k)); res_out(t, is, k) = fo(dg(k)); evid(t, is, k) = -fe(dg(k)); end
        crit = {squeeze(res_in(t, is, :)), squeeze(res_out(t, is, :)), -squeeze(evid(t, is, :))};
        fun = {fr, fo, fe};
        for c = 1:3
            [~, i0] = min(crit{c});
            best(t, is, c) = golden(fun{c}, dg(max(i0-1,1)), dg(min(i0+1,numel(dg))), 14);
        end
        best(t, is, 4) = d0; rp = mmse_rate(H, H, s2, evalb);
        for c = 1:4
            Hh = est_general(Yi, Si, shiftk(Tnom, bins, best(t, is, c), K), ones(Q,1), F, lamF, s2, 8, 20);
            lossd(t, is, c) = rp - mmse_rate(Hh, H, s2, evalb);
        end
    end
    fprintf('trial %2d d0 %+.3f | 30 dB: in %+.4f out %+.4f evid %+.4f (errors vs d0)\n', t, d0, squeeze(best(t, 2, 1:3)) - d0);
end
% projection error of shifted responses vs delay (output hop delayed, input hop advanced)
[As, Ad] = phys_channel(Q, bins, df, fc, popt);
dl = (0:0.05:2)*1e-6; pe = zeros(numel(dl), 2);
for k = 1:numel(dl)
    Dd = exp(-1j*2*pi*bins(:)*df*dl(k));
    Xo = Dd .* Ad.'; Xi = conj(Dd) .* As.';
    pe(k, :) = [norm(Xo - F*(F'*Xo), 'fro')/norm(Xo, 'fro'), norm(Xi - F*(F'*Xi), 'fro')/norm(Xi, 'fro')];
end
agree = abs(best(:, :, 2) - best(:, :, 3)) < 0.005;
R = struct('ntrial', ntrial, 'snr_db', snrs, 'd_grid', dg, 'criteria', {{'resid_in', 'heldout', 'evidence', 'true'}}, ...
    'rms_err_vs_true', squeeze(sqrt(mean((best(:, :, 1:3) - best(:, :, 4)).^2, 1))), ...
    'heldout_evidence_agree_frac', mean(agree, 1), 'gap_mean', squeeze(mean(lossd, 1)), ...
    'proj_delay_s', dl, 'proj_err_out_in', pe, ...
    'example_ch1', struct('res_in', squeeze(res_in(1, :, :)), 'res_out', squeeze(res_out(1, :, :)), 'evid', squeeze(evid(1, :, :)), 'd0', best(1, 1, 4)), ...
    'best', best, 'gap_trials', lossd);
for is = 1:2
    fprintf('%d dB: rms error vs true in/out/evid %.4f/%.4f/%.4f | heldout~evidence agree %.0f%% | gap in/out/evid/true %.3f/%.3f/%.3f/%.3f\n', ...
        snrs(is), R.rms_err_vs_true(is, :), 100*R.heldout_evidence_agree_frac(is), R.gap_mean(is, :));
end
fprintf('projection error (out, in) at 0 / 0.5 / 1 / 2 us: %s\n', sprintf('%.1e/%.1e  ', pe([1 11 21 41], :).'));
fid = fopen(fullfile(out, 'week9_timing_landscape.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function r = resid(d, Yi, Si, Yt, St, Tnom, bins, K, F, lamF, s2, Q)
Hc = est_general(Yi, Si, shiftk(Tnom, bins, d, K), ones(Q,1), F, lamF, s2, 8, 20);
r = norm(Yt - Hc*St, 'fro')^2;
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
