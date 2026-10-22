%WEEK10_FAIR_BASELINES Review step 3: C1, C4, E2, E3, E6 (pre-registered 2026-10-22).
%   Every estimator gets the same three equalisers: naive MMSE, robust MMSE with the error power estimated
%   from the pilot residual, and robust MMSE with the oracle error power. CE-BEM without a harmonic cap
%   (matrix-free LMMSE, est_bem_cg). Ideal-element model-aware LMMSE replaces the Toeplitz straw man.
%   Rates on the true H with two detector models: known G_mm (metric 13) and estimated gain (E6).
%   Channels: shared-path 30 ns and TDL-C 100 ns. 30 channels, 20/30 dB, Np = 2/4/8/16/32, pilots only.

clear; rng(10100);
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; K = 8;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
chans = {'phys30', 'TDLC100'};
cfg.phys30 = struct('tw', 0.3e-6, 'L', 16, 'r', 8, 'tw2', 0.6e-6);
cfg.TDLC100 = struct('tw', 1e-6, 'L', ceil(M*df*1e-6) + 10, 'r', 16, 'tw2', 2e-6);
Nps = [2 4 8 16 32]; snrs = [20 30]; ntrial = 30;
ests = {'onetap', 'ideal', 'bem12', 'bem40', 'bem128', 'bemall', 'omp', 'ours'};
eqs = {'naive', 'robust_resid', 'robust_oracle'};
% ---- check: matrix-free CE-BEM equals the dense estimator on a small case
[Fc, lc] = dpss_basis(bins, df, 0.6e-6, ceil(M*df*0.6e-6) + 12);
[As, Ad] = phys_channel(Q, bins, df, fc, popt); H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
S = qpsk(M, 4); Y = H*S + noise(M, 4, 1e-3);
chk = norm(est_bem_cg(Y, S, Fc, lc, 1e-3, -12:12, 1e-10, 2000) - est_bem_dense(Y, S, Fc, lc, 1e-3, -12:12), 'fro') ...
      / norm(est_bem_dense(Y, S, Fc, lc, 1e-3, -12:12), 'fro');
fprintf('matrix-free vs dense CE-BEM: relative difference %.1e\n', chk);
R = struct('Np', Nps, 'snr_db', snrs, 'ntrial', ntrial, 'estimators', {ests}, 'equalisers', {eqs}, 'bem_cg_check', chk);
for ic = 1:numel(chans)
    cf = cfg.(chans{ic});
    [F, lamF] = dpss_basis(bins, df, cf.tw, cf.L);
    L2 = ceil(M*df*cf.tw2) + 12; [F2, lam2] = dpss_basis(bins, df, cf.tw2, L2);
    rate = nan(ntrial, numel(ests), numel(eqs), numel(Nps), 2, 2);   % ... x snr x detector[known G, est gain]
    rperf = nan(ntrial, 2);
    for t = 1:ntrial
        if ic == 1, [As, Ad] = phys_channel(Q, bins, df, fc, popt);
        else, [As, Ad] = tdl38901(Q, bins, df, fc, 'C', 100e-9); end
        H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
        Pc = mean(sum(abs(H(:, evalb)).^2, 1));
        for is = 1:2
            s2 = Pc / 10^(snrs(is)/10); rperf(t, is) = rates(H, H, s2, evalb, 0) * [1; 0];
            Sall = qpsk(M, max(Nps)); Yall = H*Sall + noise(M, max(Nps), s2);
            for ip = 1:numel(Nps)
                Np = Nps(ip); S = Sall(:, 1:Np); Y = Yall(:, 1:Np);
                for ie = 1:numel(ests)
                    switch ests{ie}
                        case 'onetap', Hh = diag(mean(Y ./ S, 2));
                        case 'ideal',  Hh = est_ideal(Y, S, T, F, lamF, s2);
                        case 'bem12',  Hh = est_bem_cg(Y, S, F2, lam2, s2, -12:12);
                        case 'bem40',  Hh = est_bem_cg(Y, S, F2, lam2, s2, -40:40);
                        case 'bem128', Hh = est_bem_cg(Y, S, F2, lam2, s2, -128:128);
                        case 'bemall'
                            if Np < 8, continue; end
                            Hh = est_bem_cg(Y, S, F2, lam2, s2, -(M-1):(M-1));
                        case 'omp',    Hh = est_omp(Y, S, s2, -16:16, 0:1/(2*M*df):cf.tw2, bins, df);
                        case 'ours',   Hh = est_general(Y, S, T, ones(Q,1), F, lamF, s2, cf.r, 20);
                    end
                    se_r = max(norm(Y - Hh*S, 'fro')^2/(M*Np) - s2, 0);
                    se_o = norm(Hh - H, 'fro')^2 / M;
                    sev = [0 se_r se_o];
                    for iq = 1:3, rate(t, ie, iq, ip, is, :) = rates(Hh, H, s2, evalb, sev(iq)); end
                end
            end
        end
        fprintf('%s trial %2d: 20 dB Np=2 rate (robust_resid) %s | perfect %.2f\n', chans{ic}, t, ...
            sprintf('%.2f ', squeeze(rate(t, :, 2, 1, 1, 1))), rperf(t, 1));
    end
    G = permute(rperf, [1 3 4 5 2]) - rate(:, :, :, :, :, 1);   % gap, known-G detector
    G2 = permute(rperf, [1 3 4 5 2]) - rate(:, :, :, :, :, 2);  % gap, estimated-gain detector
    C = struct('rate_perfect', mean(rperf, 1), 'gap_knownG', squeeze(mean(G, 1, 'omitnan')), ...
               'gap_estgain', squeeze(mean(G2, 1, 'omitnan')), 'dims', 'estimator x equaliser x Np x snr');
    % E3: minimal Np with mean gap <= 0.3 at 20 dB, robust_resid, known-G detector
    npmin = cell(1, numel(ests));
    for ie = 1:numel(ests)
        g20 = squeeze(mean(G(:, ie, 2, :, 1), 1, 'omitnan'));
        k = find(g20 <= 0.3, 1); if isempty(k), npmin{ie} = '>32'; else, npmin{ie} = Nps(k); end
    end
    C.Np_min_0p3 = npmin;
    % C4: ours minus best baseline (best over CE-BEM sets, OMP, ideal and equalisers) per Np, snr
    for ip = 1:numel(Nps)
        for is = 1:2
            bl = squeeze(G(:, 2:7, :, ip, is));                % trial x 6 baselines x 3 eq
            best = min(reshape(bl, ntrial, []), [], 2, 'omitnan');
            [m, h] = tci(squeeze(G(:, 8, 2, ip, is)) - best);
            C.(sprintf('ours_minus_bestbaseline_Np%d_%ddB', Nps(ip), snrs(is))) = [m h];
            % C1 sanity: oracle-robust CE-BEM (best set) rate 20 -> 30 dB
        end
    end
    bemr = squeeze(mean(rate(:, 3:6, 3, :, :, 1), 1, 'omitnan'));   % sets x Np x snr
    C.bem_oracle_rate_20 = squeeze(bemr(:, :, 1)); C.bem_oracle_rate_30 = squeeze(bemr(:, :, 2));
    C.c1_sanity_bem_not_worse_at_30dB = all(bemr(:, :, 2) >= bemr(:, :, 1) - 0.05 | isnan(bemr(:, :, 1)), 'all');
    C.rate_trials = rate; C.rate_perfect_trials = rperf;
    R.(chans{ic}) = C;
    fprintf('%s: Np_min(0.3 @20 dB) %s\n', chans{ic}, strjoin(cellfun(@(e, n) sprintf('%s=%s', e, num2str(n)), ests, npmin, 'uni', 0), ', '));
    fprintf('%s: C1 sanity (oracle-robust CE-BEM not worse at 30 dB): %d\n', chans{ic}, C.c1_sanity_bem_not_worse_at_30dB);
    fid = fopen(fullfile(out, 'week10_fair_baselines.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);
end

%% ---------------------------------------------------------------------
function r = rates(Hh, H, s2, evalb, se2)
% [known-G, estimated-gain] rates of the (robust) MMSE built from Hh, on the true H
M = size(H, 1);
W = (Hh'*Hh + (s2 + se2)*eye(M)) \ Hh'; G = W*H; g = diag(G); gh = diag(W*Hh);
it = sum(abs(G).^2, 2) - abs(g).^2 + s2*sum(abs(W).^2, 2);
s1 = abs(g).^2 ./ it; s0 = abs(gh).^2 ./ (it + abs(g - gh).^2);
r = [mean(log2(1 + s1(evalb))), mean(log2(1 + s0(evalb)))];
end
function Hh = est_ideal(Y, S, T, F, lam, s2)
% model-aware LMMSE with the ideal +-1 element: H = (T1 - T2) o (F C F^T), one L x L delay matrix
[M, Np] = size(S); L = size(F, 2); y = Y(:); Td = T{1} - T{2};
Phi = zeros(M*Np, L*L);
for p = 1:Np
    Gs = Td * (S(:,p) .* F);
    for i = 1:L, Phi((p-1)*M+(1:M), (i-1)*L+(1:L)) = F(:, i) .* Gs; end
end
w = kron(lam(:), lam(:));
sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
Dw = sqrt(sc*w);
c = Dw .* (((Dw .* (Phi'*Phi)) .* Dw.' + s2*eye(numel(w))) \ (Dw .* (Phi'*y)));
Hh = Td .* (F*reshape(c, L, L).'*F.');
end
function Hh = est_bem_dense(Y, S, F2, lam, s2, hs)
[M, Np] = size(S); L = size(F2, 2); Nh = numel(hs); y = Y(:);
Phi = zeros(M*Np, Nh*L);
for k = 1:Nh, for p = 1:Np, Phi((p-1)*M+(1:M), (k-1)*L+(1:L)) = shr(S(:,p) .* F2, hs(k)); end, end
w = repmat(lam(:), Nh, 1);
sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
Dw = sqrt(sc*w);
c = Dw .* (((Dw .* (Phi'*Phi)) .* Dw.' + s2*eye(numel(w))) \ (Dw .* (Phi'*y)));
Hh = zeros(M);
for k = 1:Nh, Hh = Hh + diag(ones(M-abs(hs(k)),1), -hs(k)) * diag(F2 * c((k-1)*L+(1:L))); end
end
function Hh = est_omp(Y, S, s2, hs, taus, bins, df)
[M, Np] = size(S); y = Y(:);
E = exp(-1j*2*pi*bins(:)*df*taus); nh = numel(hs); nt = numel(taus);
D = zeros(M*Np, nh*nt);
for k = 1:nh, for p = 1:Np, D((p-1)*M+(1:M), (k-1)*nt+(1:nt)) = shr(S(:,p) .* E, hs(k)); end, end
Dn = D ./ vecnorm(D); res = y; sel = []; target = s2 * numel(y);
for it = 1:300
    [~, j] = max(abs(Dn' * res)); sel(end+1) = j; %#ok<AGROW>
    c = D(:, sel) \ y; res = y - D(:, sel) * c;
    if real(res'*res) <= target, break; end
end
Hh = zeros(M);
for ii = 1:numel(sel)
    k = floor((sel(ii)-1)/nt) + 1; ti = sel(ii) - (k-1)*nt;
    Hh = Hh + diag(ones(M-abs(hs(k)),1), -hs(k)) * diag(c(ii) * E(:, ti));
end
end
function Z = shr(X, h)
Z = zeros(size(X)); M = size(X, 1);
if h >= 0, Z(1+h:M, :) = X(1:M-h, :); else, Z(1:M+h, :) = X(1-h:M, :); end
end
function [m, h] = tci(x)
x = x(~isnan(x)); n = numel(x); nu = n - 1; tq = sqrt(nu*(1/betaincinv(0.05, nu/2, 0.5) - 1));
m = mean(x); h = tq*std(x)/sqrt(n);
end
function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
function W = noise(M, n, s2), W = sqrt(s2/2) * (randn(M,n) + 1j*randn(M,n)); end
