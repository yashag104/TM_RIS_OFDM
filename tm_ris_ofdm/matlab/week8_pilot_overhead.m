%WEEK8_PILOT_OVERHEAD Pilot count vs refresh interval vs accuracy: data for net throughput.
%   A user moves at 3, 30 or 120 km/h (RIS -> user paths/taps with Doppler, intra-symbol ICI
%   included, as week5_mobility). Np pilot symbols are sent at symbols 0..Np-1 (each sees
%   the channel of its own symbol); the estimate is then used for data at later symbols n.
%   The achievable rate is recorded as a function of n on a grid up to 320 symbols, so the
%   net throughput of any refresh interval N follows afterwards (figures/make_figures.py):
%       net(N, Np) = (N - Np)/N * mean_{n = Np..N-1} rate(n).
%   Channels:
%     phys     6-path Rician model used throughout (30 ns rms), 8 trials
%     TDLA/C/D 3GPP TR 38.901 TDL-A, C, D at 100 ns delay spread, per-tap Doppler, 4 trials
%   Receivers:
%     ours     est_general (timing only), Np = 1, 2, 4, 8
%              phys: window 0.3 us, rank 8; TDL: window 1 us (10 x DS), rank 16 (week6_channels_fix)
%     bem      CE-BEM (DPSS x harmonics, LMMSE), Np = 2, 4, 8, 16, 32
%              phys: |h| <= 12, 40, 128, window 0.6 us; TDL: |h| <= 12, 40, window 2.34 us
%              (|h| <= 128 with the TDL window needs ~7700 unknowns, ~1 GB per estimate)
%              all kept; the analysis takes BEM's best per case, in its favour
%     onetap   per-subcarrier LS from the pilots, Np = 1, 2, 4, 8
%     perfect  true channel at symbol n, no pilots (upper bound)
%   Two rate metrics per point:
%     cpe      MMSE SINR with |diag(WH)|^2 as signal: assumes the receiver tracks the
%              common phase/gain of each subcarrier (as with 5G phase-tracking RS)
%     nocpe    the receiver trusts its own effective gain diag(W Hhat); the mismatch
%              |diag(WH) - diag(W Hhat)|^2 counts as interference (no tracking)
%   S1, Q = 16, K = 8, M = 256, 20/30 dB.

clear; rng(4900);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3; Tsym = 1/df + 2.34e-6; Tcp = 2.34e-6;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
speeds = [3 30 120]; snrs = [20 30];
grid = [1:8 10 12 15 20 25 30 40 50 60 80 100 130 160 200 250 320];
NpO = [1 2 4 8]; NpB = [2 4 8 16 32];
chans = {'phys', 'TDLA', 'TDLC', 'TDLD'};
R = struct('speeds_kmh', speeds, 'snr_db', snrs, 'grid', grid, 'Tsym_s', Tsym, ...
    'Np_ours', NpO, 'Np_bem', NpB, 'metrics', {{'cpe', 'nocpe'}}, ...
    'dims_ours', 'Np x speed x snr x grid x metric', 'dims_bem', 'Np x hmax x speed x snr x grid x metric', ...
    'dims_perfect', 'speed x snr x grid');
nG = numel(grid);
for ic = 1:numel(chans)
    if strcmp(chans{ic}, 'phys')
        ntrial = 8; r = 8; [F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
        L2 = ceil(M*df*0.6e-6) + 12; [F2, lam2] = dpss_basis(bins, df, 0.6e-6, L2); hB = [12 40 128];
    else
        ntrial = 4; r = 16; ds = 100e-9; tw = min(Tcp, 10*ds);
        [F, lamF] = dpss_basis(bins, df, tw, ceil(M*df*tw) + 10);
        L2 = ceil(M*df*Tcp) + 12; [F2, lam2] = dpss_basis(bins, df, Tcp, L2); hB = [12 40];
    end
    rO = nan(ntrial, numel(NpO), numel(speeds), 2, nG, 2); rT = rO;
    rB = nan(ntrial, numel(NpB), numel(hB), numel(speeds), 2, nG, 2);
    rP = nan(ntrial, numel(speeds), 2, nG);
    for t = 1:ntrial
        for iv = 1:numel(speeds)
            fD = speeds(iv)/3.6 * fc / 3e8;
            if strcmp(chans{ic}, 'phys'), [As, Adp, nu] = moving_channel(Q, bins, df, fc, popt, fD);
            else, [As, Adp, nu] = tdl38901_paths(Q, bins, df, fc, chans{ic}(4), ds, fD); end
            nAll = unique([0:max(NpB)-1, grid]); Hc = containers.Map('KeyType', 'double', 'ValueType', 'any');
            for n = nAll, Hc(n) = channel_at(n, bins, fc, df, seqs, As, Adp, nu, P, M, Tsym); end
            H0 = Hc(0); Pc = mean(sum(abs(H0(:, evalb)).^2, 1));
            for is = 1:2
                s2 = Pc / 10^(snrs(is)/10);
                S = qpsk(M, max(NpB)); Y = zeros(M, max(NpB));
                for p = 1:max(NpB), Hp = Hc(p-1); Y(:, p) = Hp*S(:, p) + noise(M, 1, s2); end
                for ig = 1:nG, Hn = Hc(grid(ig)); rP(t, iv, is, ig) = rates(Hn, Hn, s2, evalb) * [1; 0]; end
                for ip = 1:numel(NpO)
                    Np = NpO(ip);
                    Ho = est_general(Y(:,1:Np), S(:,1:Np), T, ones(Q,1), F, lamF, s2, r, 20);
                    Ht = diag(mean(Y(:,1:Np) ./ S(:,1:Np), 2));
                    for ig = find(grid >= Np)
                        Hn = Hc(grid(ig));
                        rO(t, ip, iv, is, ig, :) = rates(Ho, Hn, s2, evalb);
                        rT(t, ip, iv, is, ig, :) = rates(Ht, Hn, s2, evalb);
                    end
                end
                for ip = 1:numel(NpB)
                    Np = NpB(ip);
                    for ih = 1:numel(hB)
                        Hb = est_bem(Y(:,1:Np), S(:,1:Np), F2, lam2, s2, -hB(ih):hB(ih));
                        for ig = find(grid >= Np)
                            rB(t, ip, ih, iv, is, ig, :) = rates(Hb, Hc(grid(ig)), s2, evalb);
                        end
                    end
                end
            end
            fprintf('%s trial %d, %d km/h done\n', chans{ic}, t, speeds(iv));
        end
    end
    R.(chans{ic}) = struct('ntrial', ntrial, 'rank', r, 'hmax_bem', hB, ...
        'ours', mean(rO, 1, 'omitnan'), 'onetap', mean(rT, 1, 'omitnan'), ...
        'bem', mean(rB, 1, 'omitnan'), 'perfect', mean(rP, 1, 'omitnan'));
    fid = fopen(fullfile(out, 'week8_pilot_overhead.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);
    ig = find(grid == 20); C = R.(chans{ic});
    for iv = 1:numel(speeds)
        fprintf('%s %3d km/h, rate at symbol 20 (cpe) @20 dB: ours Np=2 %.2f Np=8 %.2f | bem Np=32 %.2f | onetap %.2f | perfect %.2f\n', ...
            chans{ic}, speeds(iv), C.ours(1, 2, iv, 1, ig, 1), C.ours(1, 4, iv, 1, ig, 1), ...
            max(C.bem(1, 5, :, iv, 1, ig, 1)), C.onetap(1, 2, iv, 1, ig, 1), C.perfect(1, iv, 1, ig));
    end
end

%% ---------------------------------------------------------------------
function r = rates(Hh, H, s2, evalb)
% [cpe, nocpe] rates of MMSE with the estimate Hh on the true channel H
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; g = diag(G); gh = diag(W*Hh);
it = sum(abs(G).^2, 2) - abs(g).^2 + s2*sum(abs(W).^2, 2);
s1 = abs(g).^2 ./ it;
s0 = abs(gh).^2 ./ (it + abs(g - gh).^2);
r = [mean(log2(1 + s1(evalb))), mean(log2(1 + s0(evalb)))];
end
function [As, Adp, nu] = moving_channel(Q, bins, df, fc, popt, fD)
lam = 3e8 / fc; k = 2*pi/lam; n = round(sqrt(Q));
[gx, gz] = meshgrid(0:n-1, 0:n-1); r = [gx(:), zeros(Q,1), gz(:)] * lam/2;
[As, ~] = phys_channel(Q, bins, df, fc, popt);
P = popt.P; Kr = 10^(popt.K_db/10);
tau = [0, sort(-popt.rms*log(rand(1, P)))];
pw = [Kr/(Kr+1), (1/(Kr+1)) * exp(-tau(2:end)/popt.rms) / sum(exp(-tau(2:end)/popt.rms))];
g = sqrt(pw) .* exp(1j*2*pi*rand(1, P+1));
u = randn(P+1, 3); u(:,2) = abs(u(:,2)); u = u ./ vecnorm(u, 2, 2);
ph = exp(-1j*k*(r*u.'));
Adp = zeros(Q, numel(bins), P+1);
for l = 1:P+1, Adp(:, :, l) = (ph(:, l) * g(l)) * exp(-1j*2*pi*tau(l)*(bins(:).'*df)); end
nu = fD * cos(2*pi*rand(1, P+1));
end
function H = channel_at(n, bins, fc, df, seqs, As, Adp, nu, P, M, Tsym)
H = zeros(M);
for l = 1:size(Adp, 3)
    Ad = Adp(:, :, l) * exp(1j*2*pi*nu(l)*n*Tsym);
    H = H + ici_matrix(M, nu(l)/df) * coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
end
end
function C = ici_matrix(M, eps)
k = (0:M-1).' - (0:M-1); x = eps - k;
C = exp(1j*pi*x*(M-1)/M) .* snc(x) ./ snc(x/M);
end
function y = snc(x), y = ones(size(x)); nz = abs(x) > 1e-12; y(nz) = sin(pi*x(nz)) ./ (pi*x(nz)); end
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
