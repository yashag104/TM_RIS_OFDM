%WEEK9_PRIOR_SENSITIVITY Review step 1, test T4 (workstream B3). Pre-registered 2026-10-22.
%   Mode (ii) physical timing (RIS locked to the BS clock: offset d = -tau_b/T_slot, bulk delays on both
%   hops, receiver synced to the first cascaded arrival), 30 channels, 30 dB, Np = 2.
%   The evidence search and the estimator use a DPSS prior window [t0, t0 + tau_max]:
%     default t0 = 0, tau_max = 0.3 us; variants t0 = -100 ns, t0 = +100 ns, tau_max = 0.15 us, 0.45 us
%   (L = 16 frozen). Known-offset receiver uses the default prior.
%   PASS: rate loss vs known offset within 0.5 dB (0.165 bit/s/Hz at 30 dB) for all four variants.

clear; rng(9400);
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; Tu = 1/df; K = 8; Tslot = Tu/K;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16; Np = 2;
base = [ones(1, 6), 2*ones(1, 2)]; st = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
E0 = (0:K)/K;
Tnom = {coupling_edges(bins, fc, df, E0, base, ones(1,M), ones(1,M), 'state1', M-1), ...
        coupling_edges(bins, fc, df, E0, base, ones(1,M), ones(1,M), 'state2', M-1)};
ramp = @(tau) exp(-1j*2*pi*bins*df*tau);
vars = struct('name', {'default', 't0_m100ns', 't0_p100ns', 'tau_0p15', 'tau_0p45'}, ...
              't0', {0, -100e-9, 100e-9, 0, 0}, 'tau', {0.3e-6, 0.3e-6, 0.3e-6, 0.15e-6, 0.45e-6});
Fv = cell(1, numel(vars)); Lv = Fv;
for v = 1:numel(vars)
    [Fv{v}, Lv{v}] = dpss_basis(bins, df, vars(v).tau, 16);
    Fv{v} = exp(-1j*2*pi*bins(:)*df*vars(v).t0) .* Fv{v};     % delays [t0, t0 + tau]
end
snr = 30; ntrial = 30; grid = -0.5:0.025:0.5; slope = 0.33;
gap = nan(ntrial, numel(vars) + 2); dh = nan(ntrial, numel(vars)); dtr = nan(ntrial, 1);
for t = 1:ntrial
    [As0, Ad0] = phys_channel(Q, bins, df, fc, popt);
    tb = 50e-9 + 250e-9*rand; tu = 30e-9 + 270e-9*rand;
    sync = conj(ramp(tb + tu)).';
    H = sync .* coupling_edges(bins, fc, df, repmat(E0, Q, 1), st, As0 .* ramp(tb), Ad0 .* ramp(tu), P, M-1);
    dtr(t) = -tb/Tslot;
    Pc = mean(sum(abs(H(:, evalb)).^2, 1)); s2 = Pc / 10^(snr/10);
    S = qpsk(M, Np); Y = H*S + noise(M, Np, s2); rp = mmse_rate(H, H, s2, evalb);
    for v = 1:numel(vars)
        fe = @(d) -evidence(d, Y, S, Tnom, bins, K, Fv{v}, Lv{v}, s2);
        e = arrayfun(fe, grid); [~, i0] = min(e);
        dh(t, v) = golden(fe, grid(max(i0-1,1)), grid(min(i0+1,numel(grid))), 14);
        gap(t, v) = rp - mmse_rate(est_general(Y, S, shiftk(Tnom, bins, dh(t, v), K), ones(Q,1), Fv{v}, Lv{v}, s2, 8, 20), H, s2, evalb);
    end
    gap(t, end-1) = rp - mmse_rate(est_general(Y, S, shiftk(Tnom, bins, dtr(t), K), ones(Q,1), Fv{1}, Lv{1}, s2, 8, 20), H, s2, evalb);
    gap(t, end) = rp - mmse_rate(est_general(Y, S, Tnom, ones(Q,1), Fv{1}, Lv{1}, s2, 8, 20), H, s2, evalb);
    fprintf('trial %2d: d_true %+.4f | dhat %s\n', t, dtr(t), sprintf('%+.4f ', dh(t, :)));
end
R = struct('ntrial', ntrial, 'snr_db', snr, 'variants', {{vars.name}}, 'gap_mean', mean(gap, 1), ...
    'columns', {[{vars.name}, {'known', 'nominal_no_search'}]}, 'dhat_bias', mean(dh - dtr, 1), 'dhat_rms', sqrt(mean((dh - dtr).^2, 1)));
for v = 1:numel(vars)
    [m, h] = tci(gap(:, v) - gap(:, end-1));
    R.(['loss_vs_known_' vars(v).name]) = [m h];
    fprintf('%-10s loss vs known %+.3f +-%.3f bit/s/Hz (%+.2f dB) | dhat bias %+.4f rms %.4f slot\n', vars(v).name, m, h, m/slope, R.dhat_bias(v), R.dhat_rms(v));
end
fprintf('known %.3f | nominal (no search) %.3f\n', mean(gap(:, end-1)), mean(gap(:, end)));
R.gap_trials = gap; R.dhat = dh; R.dtrue = dtr;
fid = fopen(fullfile(out, 'week9_prior_sensitivity.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);
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
