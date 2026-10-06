%WEEK6_BASELINES_PILOTS (#8 baselines, #11 pilot design)
%   (#8) estimators at Np = 2 full pilot symbols: ours (timing-only); CE-BEM (Q=12,
%        LMMSE); OMP over a harmonic x delay dictionary (compressed sensing, stops at
%        the noise level); plus unstructured banded LS with 32 orthogonal pilot symbols
%        (16x our pilot overhead) as a brute-force reference.
%   (#11) pilot layouts with the same number of pilot resource elements (512):
%        full    2 symbols, every subcarrier
%        comb2   4 symbols, every 2nd subcarrier, other REs empty (NR DMRS-like)
%        comb4   8 symbols, every 4th subcarrier, other REs empty
%        comb2d  4 symbols, every 2nd subcarrier, DATA on the other REs (unknown to the
%                estimator: it leaks into pilot bins through the switching), plus
%        comb2dd the same with one decision-directed pass (data decided and reused)
%   S1, D = 0.75, K = 8, Q = 16, M = 256, physical channel, 10 trials, 20/30 dB.

clear; rng(3200);
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
snrs = [20 30]; ntrial = 10;
bnames = {'ours', 'bem', 'omp', 'ls32', 'perfect'};
pnames = {'full', 'comb2', 'comb4', 'comb2d', 'comb2dd'};
rb = zeros(ntrial, numel(bnames), 2); rpil = zeros(ntrial, numel(pnames), 2);
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
    Pc = mean(sum(abs(H(:, evalb)).^2, 1));
    for is = 1:2
        s2 = Pc / 10^(snrs(is)/10);
        rp = mmse_rate(H, H, s2, evalb);
        % ---- #8 baselines, 2 full pilot symbols
        S = qpsk(M, 2); Y = H*S + noise(M, 2, s2);
        rb(t, 1, is) = mmse_rate(est_general(Y, S, T, ones(Q,1), F, lamF, s2, 8, 20), H, s2, evalb);
        rb(t, 2, is) = mmse_rate(est_bem(Y, S, F2, lam2, s2, -12:12), H, s2, evalb);
        rb(t, 3, is) = mmse_rate(est_omp(Y, S, s2, -16:16, 0:1/(2*M*df):0.6e-6, bins, df), H, s2, evalb);
        S32 = qpsk(M, 32); Y32 = H*S32 + noise(M, 32, s2);
        rb(t, 4, is) = mmse_rate(est_banded(Y32, S32, -15:15), H, s2, evalb);
        rb(t, 5, is) = rp;
        % ---- #11 pilot layouts, 512 pilot REs each
        for ip = 1:numel(pnames)
            switch pnames{ip}
                case 'full',   Sp = qpsk(M, 2); Xd = zeros(M, 2);
                case 'comb2',  Sp = comb(M, 4, 2); Xd = zeros(M, 4);
                case 'comb4',  Sp = comb(M, 8, 4); Xd = zeros(M, 8);
                otherwise,     Sp = comb(M, 4, 2); Xd = qam16(M, 4) .* (Sp == 0);
            end
            Yp = H*(Sp + Xd) + noise(M, size(Sp, 2), s2);
            Hh = est_general(Yp, Sp, T, ones(Q,1), F, lamF, s2, 8, 20);
            if strcmp(pnames{ip}, 'comb2dd')
                Xh = qamdec((Hh'*Hh + s2*eye(M)) \ (Hh'*Yp)) .* (Sp == 0);
                Hh = est_general(Yp, Sp + Xh, T, ones(Q,1), F, lamF, s2, 8, 20);
            end
            rpil(t, ip, is) = mmse_rate(Hh, H, s2, evalb);
        end
        fprintf('trial %d snr %d done\n', t, snrs(is));
    end
end
R = struct('snr_db', snrs, 'baseline_names', {bnames}, 'baseline_rate', squeeze(mean(rb, 1)), ...
    'baseline_ci95', squeeze(1.96*std(rb, 0, 1)/sqrt(ntrial)), 'pilot_names', {pnames}, ...
    'pilot_rate', squeeze(mean(rpil, 1)), 'pilot_ci95', squeeze(1.96*std(rpil, 0, 1)/sqrt(ntrial)));
for is = 1:2
    fprintf('#8  %d dB: %s\n', snrs(is), strjoin(arrayfun(@(i) sprintf('%s %.2f', bnames{i}, R.baseline_rate(i, is)), 1:numel(bnames), 'uni', 0), ' | '));
    fprintf('#11 %d dB: %s\n', snrs(is), strjoin(arrayfun(@(i) sprintf('%s %.2f', pnames{i}, R.pilot_rate(i, is)), 1:numel(pnames), 'uni', 0), ' | '));
end
fid = fopen(fullfile(out, 'week6_baselines_pilots.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

%% ---------------------------------------------------------------------
function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
function W = noise(M, n, s2), W = sqrt(s2/2) * (randn(M,n) + 1j*randn(M,n)); end
function S = comb(M, n, step)
% n symbols; symbol i carries QPSK on subcarriers offset mod(i-1, step) + 1 : step : M.
S = zeros(M, n); q = qpsk(M, n);
for i = 1:n, idx = mod(i-1, step) + 1 : step : M; S(idx, i) = q(idx, i); end          % equal energy per pilot RE (no NR power boost)
end
function X = qam16(M, n)
lv = [-3 -1 1 3] / sqrt(10);
X = lv(randi(4, M, n)) + 1j*lv(randi(4, M, n));
end
function X = qamdec(Z)
lv = [-3 -1 1 3] / sqrt(10);
q = @(v) lv(min(max(round((v*sqrt(10) + 3)/2), 0), 3) + 1);
X = reshape(q(real(Z)) + 1j*q(imag(Z)), size(Z));
end

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

function Hh = est_omp(Y, S, s2, hs, taus, bins, df)
% OMP over atoms (h, tau): H = sum S_h diag(c e^{-j2pi m df tau}). Stops when the
% residual reaches the noise energy, or after 300 atoms.
[M, Np] = size(S); y = Y(:);
E = exp(-1j*2*pi*bins(:)*df*taus);                         % M x Nt
nh = numel(hs); nt = numel(taus);
D = zeros(M*Np, nh*nt);
for k = 1:nh, for p = 1:Np
    D((p-1)*M+(1:M), (k-1)*nt+(1:nt)) = shr(S(:,p) .* E, hs(k));
end, end
nrm = vecnorm(D); Dn = D ./ nrm;
res = y; sel = []; target = s2 * numel(y);
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

function Hh = est_banded(Y, S, hs)
% Unstructured LS of every band entry (needs Np >= numel(hs) pilots).
M = size(S, 1); Hh = zeros(M);
for r = 1:M
    cols = r - hs; ok = cols >= 1 & cols <= M; cols = cols(ok);
    Hh(r, cols) = (S(cols, :).' \ Y(r, :).').';
end
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
