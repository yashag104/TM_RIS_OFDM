%WEEK9_TIMING_PHYSICS Review step 1, tests T1 and T2 (workstreams A2, B1). Pre-registered in
%   results/PREREG_LOG.md (2026-10-22) before running.
%   Physical timing: each hop gets a bulk delay, tau_b ~ U(50,300) ns (BS -> RIS) and
%   tau_u ~ U(30,300) ns (RIS -> user), added to every path. The receiver's FFT window is
%   synced to the first cascaded arrival (tau_b + tau_u; perfect sync assumed). The ground
%   truth uses ABSOLUTE times: transmit symbol and switching edges referenced to t = 0, hop
%   responses with their bulk delays, then the sync ramp exp(+j 2 pi f_mbar (tau_b+tau_u)).
%   Two controller modes:
%     (i)  RIS locked to the incident signal: edges at (k/K + tau_b/T_u), optional per-edge jitter
%     (ii) RIS locked to the BS clock:        edges at k/K
%   T1: identities (i) synced H == no-bulk H; (ii) synced H == no-bulk H with offset d = -tau_b/T_slot.
%   T2: receivers nominal (no search), evidence (Alg. 2), known offset, perfect; 30 channels,
%       20/30 dB; paired t-intervals against the matched (no-bulk, nominal) result.
%   Frozen estimator: tau_max = 0.3 us, L = 16, r = 8, 20 iterations, Np = 2.

clear; rng(9100);
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; Tu = 1/df; K = 8; Tslot = Tu/K;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16; Np = 2;
base = [ones(1, 6), 2*ones(1, 2)]; st = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
E0 = (0:K)/K;
Tnom = {coupling_edges(bins, fc, df, E0, base, ones(1,M), ones(1,M), 'state1', M-1), ...
        coupling_edges(bins, fc, df, E0, base, ones(1,M), ones(1,M), 'state2', M-1)};
ramp = @(tau) exp(-1j*2*pi*bins*df*tau);                    % 1 x M delay ramp (baseband)
snrs = [20 30]; ntrial = 30; grid = -0.5:0.025:0.5;
modes = {'i', 'i_jit10ps', 'i_jit1ns', 'ii'}; jit = [0 10e-12 1e-9 0];
recv = {'nominal', 'evidence', 'known', 'perfect'};
t1err = zeros(ntrial, 2);
rate = nan(ntrial, numel(modes), numel(recv), 2); rmatch = nan(ntrial, 2); dhat = nan(ntrial, numel(modes), 2);
for t = 1:ntrial
    [As0, Ad0] = phys_channel(Q, bins, df, fc, popt);
    tb = 50e-9 + 250e-9*rand; tu = 30e-9 + 270e-9*rand;
    As = As0 .* ramp(tb); Ad = Ad0 .* ramp(tu);
    sync = conj(ramp(tb + tu)).';                           % M x 1, advance by tau_b + tau_u
    H0 = coupling_edges(bins, fc, df, repmat(E0, Q, 1), st, As0, Ad0, P, M-1);
    % ---- T1 identities
    Hi = sync .* coupling_edges(bins, fc, df, repmat(E0 + tb/Tu, Q, 1), st, As, Ad, P, M-1);
    Hii = sync .* coupling_edges(bins, fc, df, repmat(E0, Q, 1), st, As, Ad, P, M-1);
    d_ii = -tb/Tslot;
    H0off = coupling_edges(bins, fc, df, repmat(((0:K) + d_ii)/K, Q, 1), st, As0, Ad0, P, M-1);
    t1err(t, :) = [norm(Hi - H0, 'fro'), norm(Hii - H0off, 'fro')] / norm(H0, 'fro');
    % ---- T2 performance
    Hm = cell(1, numel(modes)); dtrue = zeros(1, numel(modes));
    for im = 1:numel(modes)
        switch modes{im}
            case 'ii', Hm{im} = Hii; dtrue(im) = d_ii;
            otherwise
                Ej = repmat(E0 + tb/Tu, Q, 1) + jit(im)/Tu*randn(Q, K+1);
                Ej(:, end) = Ej(:, 1) + 1;                  % keep the period exactly one
                Hm{im} = sync .* coupling_edges(bins, fc, df, Ej, st, As, Ad, P, M-1); dtrue(im) = 0;
        end
    end
    Pc = mean(sum(abs(H0(:, evalb)).^2, 1));
    for is = 1:2
        s2 = Pc / 10^(snrs(is)/10);
        S = qpsk(M, Np); W = noise(M, Np, s2);
        Y0 = H0*S + W; rp0 = mmse_rate(H0, H0, s2, evalb);
        rmatch(t, is) = mmse_rate(est_general(Y0, S, Tnom, ones(Q,1), F, lamF, s2, 8, 20), H0, s2, evalb) - rp0;
        for im = 1:numel(modes)
            H = Hm{im}; Y = H*S + W; rp = mmse_rate(H, H, s2, evalb);
            fe = @(d) -evidence(d, Y, S, Tnom, bins, K, F, lamF, s2);
            e = arrayfun(fe, grid); [~, i0] = min(e);
            de = golden(fe, grid(max(i0-1,1)), grid(min(i0+1,numel(grid))), 14);
            dhat(t, im, is) = de;
            Hh = {est_general(Y, S, Tnom, ones(Q,1), F, lamF, s2, 8, 20), ...
                  est_general(Y, S, shiftk(Tnom, bins, de, K), ones(Q,1), F, lamF, s2, 8, 20), ...
                  est_general(Y, S, shiftk(Tnom, bins, dtrue(im), K), ones(Q,1), F, lamF, s2, 8, 20), H};
            for ir = 1:4, rate(t, im, ir, is) = mmse_rate(Hh{ir}, H, s2, evalb) - rp; end   % = -gap
        end
    end
    fprintf('trial %2d: T1 err %.1e / %.1e | tau_b %3.0f ns | dhat(ii) %.4f vs %.4f @30\n', t, t1err(t,:), tb*1e9, dhat(t,4,2), d_ii);
end
gap = -rate; gm = -rmatch;
R = struct('ntrial', ntrial, 'snr_db', snrs, 'modes', {modes}, 'receivers', {recv}, ...
    't1_rel_err_max', max(t1err, [], 1), 'gap_matched', mean(gm, 1));
for is = 1:2
    for im = 1:numel(modes)
        key = sprintf('%s_%ddB', modes{im}, snrs(is));
        g = squeeze(gap(:, im, :, is));                           % ntrial x 4
        [dn, cn] = tci(g(:,1) - gm(:, is)); [de_, ce] = tci(g(:,2) - gm(:, is)); [dk, ck] = tci(g(:,2) - g(:,3));
        R.(key) = struct('gap_mean', mean(g, 1), 'nominal_minus_matched', [dn cn], 'evidence_minus_matched', [de_ ce], ...
                         'evidence_minus_known', [dk ck], 'dhat_mean', mean(dhat(:, im, is)), 'dhat_sd', std(dhat(:, im, is)));
        fprintf('%-10s %d dB: gap nominal/evid/known/perf %.3f/%.3f/%.3f/%.3f | matched %.3f | evid-matched %+.3f +-%.3f | evid-known %+.3f +-%.3f\n', ...
            modes{im}, snrs(is), mean(g, 1), mean(gm(:, is)), de_, ce, dk, ck);
    end
end
R.gap_trials = gap; R.gap_matched_trials = gm; R.dhat = dhat;
fid = fopen(fullfile(out, 'week9_timing_physics.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

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
