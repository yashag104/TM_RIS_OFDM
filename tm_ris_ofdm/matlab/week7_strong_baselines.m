%WEEK7_STRONG_BASELINES (points 2b, 2c) The strongest honest baselines.
%   2b) CE-BEM with ALL harmonics (|h| <= 255: first 4 trials only, ~8700 unknowns). Our earlier BEM kept |h| <= 12 (25 harmonics), so a
%       reviewer could say it failed only because it was truncated. Here BEM keeps
%       |h| <= 12, 40, 80 and ALL (|h| <= M-1 = 255), with the same LMMSE/DPSS prior.
%       Many harmonics need many pilots, so each BEM is also given Np = 2, 8, 32.
%   2c) ICI equalisers from the high-mobility OFDM literature, all fed the SAME
%       estimate (ours, and the true H as an oracle), compared with our full MMSE:
%         banded MMSE  (Schniter 2004 / Rugini 2005 style: per subcarrier, window +/-Q)
%         LSQR / CG    (Hrycak et al. 2010 style: iterative, matrix-free in principle)
%         one-tap
%   Oracle: MMSE with the true H throughout.
%   Physical channel, S1, Q = 16, K = 8, M = 256, 20 trials, 20 / 30 dB.

clear; rng(3800);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
L2 = ceil(M*df*0.6e-6) + 12; [F2, lam2] = dpss_basis(bins, df, 0.6e-6, L2);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
snrs = [20 30]; ntrial = 20;
Hmaxes = [12 40 128 M-1]; Nps = [2 8 32];   % ALL (M-1) only on the first 4 trials (cost)
eqnames = {'onetap', 'banded4', 'banded16', 'banded32', 'lsqr10', 'lsqr30', 'full_mmse'};
rb = zeros(ntrial, numel(Hmaxes), numel(Nps), 2); re = zeros(ntrial, numel(eqnames), 2, 2);  % est x [ours, oracle]
rours = zeros(ntrial, 2); rperf = zeros(ntrial, 2);
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
    Pc = mean(sum(abs(H(:, evalb)).^2, 1));
    for is = 1:2
        s2 = Pc / 10^(snrs(is)/10);
        rperf(t, is) = mmse_rate(H, H, s2, evalb);
        S2 = qpsk(M, 2); Y2 = H*S2 + noise(M, 2, s2);
        Ho = est_general(Y2, S2, T, ones(Q,1), F, lamF, s2, 8, 20);
        rours(t, is) = mmse_rate(Ho, H, s2, evalb);
        % ---- 2b: BEM with more harmonics and more pilots
        for ip = 1:numel(Nps)
            S = qpsk(M, Nps(ip)); Y = H*S + noise(M, Nps(ip), s2);
            for ih = 1:numel(Hmaxes)
                if Hmaxes(ih) == M-1 && t > 4, rb(t, ih, ip, is) = NaN; continue; end
                Hb = est_bem(Y, S, F2, lam2, s2, -Hmaxes(ih):Hmaxes(ih));
                rb(t, ih, ip, is) = mmse_rate(Hb, H, s2, evalb);
            end
        end
        % ---- 2c: literature ICI equalisers on our estimate and on the oracle
        for src = 1:2
            if src == 1, He = Ho; else, He = H; end
            for ie = 1:numel(eqnames)
                re(t, ie, is, src) = eq_rate(eqnames{ie}, He, H, s2, evalb);
            end
        end
    end
    fprintf('trial %d done\n', t);
end
R = struct('snr_db', snrs, 'ntrial', ntrial, 'bem_hmax', Hmaxes, 'bem_Np', Nps, ...
    'bem_rate', squeeze(mean(rb, 1, 'omitnan')), 'ours_Np2', mean(rours, 1), 'perfect', mean(rperf, 1), ...
    'eq_names', {eqnames}, 'eq_rate_ours_oracle', squeeze(mean(re, 1)));
for is = 1:2
    fprintf('SNR %d dB: ours (Np=2) %.2f | perfect %.2f\n', snrs(is), R.ours_Np2(is), R.perfect(is));
    for ih = 1:numel(Hmaxes)
        fprintf('   BEM |h|<=%3d: Np=2/8/32 -> %s\n', Hmaxes(ih), sprintf('%.2f ', squeeze(R.bem_rate(ih, :, is))));
    end
    for ie = 1:numel(eqnames)
        fprintf('   %-9s on ours %.2f | on oracle %.2f\n', eqnames{ie}, R.eq_rate_ours_oracle(ie, is, 1), R.eq_rate_ours_oracle(ie, is, 2));
    end
end
fid = fopen(fullfile(out, 'week7_strong_baselines.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

%% ---------------------------------------------------------------------
function r = eq_rate(name, He, H, s2, evalb)
% Linear equaliser W built from the estimate He, scored on the true H (exact SINR).
M = size(H, 1);
switch name
    case 'onetap'
        W = diag(1 ./ diag(He));
    case 'full_mmse'
        W = (He'*He + s2*eye(M)) \ He';
    otherwise
        if startsWith(name, 'banded')
            Wd = str2double(name(7:end)); W = zeros(M);
            for m = 1:M
                rows = max(1, m-Wd):min(M, m+Wd); cols = max(1, m-2*Wd):min(M, m+2*Wd);
                Hw = He(rows, cols); k = find(cols == m);
                w = (Hw*Hw' + s2*eye(numel(rows))) \ Hw(:, k);
                W(m, rows) = w';
            end
        else                                                   % lsqrN: N CG iterations
            nit = str2double(name(5:end));
            W = cg_matrix(He, s2, nit);
        end
end
G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end

function W = cg_matrix(He, s2, nit)
% The linear map y -> x_hat produced by nit CG iterations on (He'He + s2 I) x = He' y
% (CG with a zero start is linear in y for a fixed iteration count; built column-wise).
M = size(He, 1); A = He'*He + s2*eye(M); B = He';
X = zeros(M); Rr = B; Pp = Rr; rs = sum(conj(Rr).*Rr, 1);
for it = 1:nit
    Ap = A*Pp; al = rs ./ sum(conj(Pp).*Ap, 1);
    X = X + al.*Pp; Rr = Rr - al.*Ap; rn = sum(conj(Rr).*Rr, 1);
    Pp = Rr + (rn./rs).*Pp; rs = rn;
end
W = X;
end

function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
function W = noise(M, n, s2), W = sqrt(s2/2) * (randn(M,n) + 1j*randn(M,n)); end

function Hh = est_bem(Y, S, F2, lam, s2, hs)
[M, Np] = size(S); L = size(F2, 2); Nh = numel(hs); y = Y(:);
Phi = zeros(M*Np, Nh*L);
for k = 1:Nh, for p = 1:Np
    Phi((p-1)*M+(1:M), (k-1)*L+(1:L)) = shr(S(:,p) .* F2, hs(k));
end, end
w = repmat(lam(:), Nh, 1);
sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
Dw = sqrt(sc*w);
c = Dw .* (((Dw .* (Phi'*Phi)) .* Dw.' + s2*eye(numel(w))) \ (Dw .* (Phi'*y)));
Hh = zeros(M);
for k = 1:Nh, Hh = Hh + diag(ones(M-abs(hs(k)),1), -hs(k)) * diag(F2 * c((k-1)*L+(1:L))); end
end
function Z = shr(X, h)
Z = zeros(size(X)); M = size(X, 1);
if h >= 0, Z(1+h:M, :) = X(1:M-h, :); else, Z(1:M+h, :) = X(1-h:M, :); end
end
function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
