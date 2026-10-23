%WEEK11_OVERHEAD_HELDOUT_PART Same as WEEK11_OVERHEAD_HELDOUT (identical seeds per channel index), run for one channel
%   index icsel (set before calling) and written to week11_overhead_heldout_<chan>.json. Used to parallelise.
%   Review E4 (pre-registered 2026-10-22): net throughput after pilot overhead with
%   the refresh interval N and pilot count Np chosen on a SELECTION set of channels and reported on a
%   disjoint HELD-OUT set (20 + 20 channels per case), with t-intervals, 20 and 30 dB.
%   Channels (CP-valid per d1_cp_check): shared-path 30 ns, TDL-A/C/D 100 ns. Speeds 3/30/120 km/h.
%   Receivers: ours (Np 1/2/4/8); CE-BEM |h| <= 128 and ALL harmonics (Np 8/16/32, matrix-free LMMSE);
%   one-tap (Np 1/2/4/8); perfect. ALL estimated receivers use the robust MMSE (error power from the pilot
%   residual). Three detector metrics: per-subcarrier genie, common-phase tracking (main), none.
%   The channel at symbol n is sum_l exp(j 2 pi nu_l n Ts) C_l H_l with per-tap coupling matrices H_l and
%   Doppler ICI matrices C_l precomputed (exact, by linearity). Selection is done in figures/ (Python).

if ~exist('icsel', 'var'), error('set icsel'); end
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; Tsym = 1/df + 2.34e-6; Tcp = 2.34e-6;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16; K = 8;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
speeds = [3 30 120]; snrs = [20 30]; nset = 20;
grid = [1:8 10 12 15 20 25 30 40 50 60 80 100 130 160 200 250 320];
NpO = [1 2 4 8]; NpB = [8 16 32]; hB = {-128:128, -(M-1):(M-1)}; nG = numel(grid);
chans = {'phys30', 'TDLA100', 'TDLC100', 'TDLD100'};
R = struct('speeds_kmh', speeds, 'snr_db', snrs, 'grid', grid, 'Np_ours', NpO, 'Np_bem', NpB, 'hmax_bem', [128 M-1], ...
           'metrics', {{'sub', 'cpe', 'none'}}, 'nset', nset, 'sets', {{'selection', 'heldout'}});
for ic = icsel
    if ic == 1, tw = 0.3e-6; L = 16; r = 8; tw2 = 0.6e-6;
    else, tw = 1e-6; L = ceil(M*df*tw) + 10; r = 16; tw2 = 2e-6; end
    [F, lamF] = dpss_basis(bins, df, tw, L);
    L2 = ceil(M*df*tw2) + 12; [F2, lam2] = dpss_basis(bins, df, tw2, L2);
    C = struct();
    for iset = 1:2
        rng(11400 + 100*ic + 10*iset);                       % disjoint seeds for selection / held-out
        rO = nan(nset, numel(NpO), numel(speeds), 2, nG, 3); rT = rO;
        rB = nan(nset, numel(NpB), numel(hB), numel(speeds), 2, nG, 3); rP = nan(nset, numel(speeds), 2, nG);
        for t = 1:nset
            for iv = 1:numel(speeds)
                fD = speeds(iv)/3.6 * fc / 3e8;
                if ic == 1, [As, Adp, nu] = moving_phys(Q, bins, df, fc, popt, fD);
                else, [As, Adp, nu] = tdl38901_paths(Q, bins, df, fc, chans{ic}(4), 100e-9, fD); end
                nl = size(Adp, 3); CH = cell(1, nl);
                for l = 1:nl, CH{l} = ici_matrix(M, nu(l)/df) * coupling_matrix(bins, fc, df, seqs, As, Adp(:, :, l), P, M-1); end
                Hat = @(n) hsum(CH, nu, n, Tsym);
                H0 = Hat(0); Pc = mean(sum(abs(H0(:, evalb)).^2, 1));
                Hg = arrayfun(@(n) Hat(n), grid, 'UniformOutput', false);
                for is = 1:2
                    s2 = Pc / 10^(snrs(is)/10);
                    S = qpsk(M, max(NpB)); Y = zeros(M, max(NpB));
                    for p = 1:max(NpB), Y(:, p) = Hat(p-1)*S(:, p) + noise(M, 1, s2); end
                    for ig = 1:nG, rP(t, iv, is, ig) = rates(Hg{ig}, Hg{ig}, s2, evalb, 0) * [1; 0; 0]; end
                    for ip = 1:numel(NpO)
                        Np = NpO(ip); Sp = S(:, 1:Np); Yp = Y(:, 1:Np);
                        Ho = est_general(Yp, Sp, T, ones(Q,1), F, lamF, s2, r, 20); seO = resid_se(Yp, Sp, Ho, s2);
                        Ht = diag(mean(Yp ./ Sp, 2)); seT = resid_se(Yp, Sp, Ht, s2);
                        for ig = find(grid >= Np)
                            rO(t, ip, iv, is, ig, :) = rates(Ho, Hg{ig}, s2, evalb, seO);
                            rT(t, ip, iv, is, ig, :) = rates(Ht, Hg{ig}, s2, evalb, seT);
                        end
                    end
                    for ip = 1:numel(NpB)
                        Np = NpB(ip); Sp = S(:, 1:Np); Yp = Y(:, 1:Np);
                        for ih = 1:numel(hB)
                            Hb = est_bem_cg(Yp, Sp, F2, lam2, s2, hB{ih}); seB = resid_se(Yp, Sp, Hb, s2);
                            for ig = find(grid >= Np), rB(t, ip, ih, iv, is, ig, :) = rates(Hb, Hg{ig}, s2, evalb, seB); end
                        end
                    end
                end
            end
            fprintf('%s set %d channel %2d done\n', chans{ic}, iset, t);
        end
        C.(R.sets{iset}) = struct('ours', rO, 'onetap', rT, 'bem', rB, 'perfect', rP);
    end
    R.(chans{ic}) = C;
    fid = fopen(fullfile(out, sprintf('week11_overhead_heldout_%s.json', chans{ic})), 'w'); fprintf(fid, '%s', jsonencode(R)); fclose(fid);
end

%% ---------------------------------------------------------------------
function H = hsum(CH, nu, n, Ts)
H = zeros(size(CH{1}));
for l = 1:numel(CH), H = H + exp(1j*2*pi*nu(l)*n*Ts) * CH{l}; end
end
function se = resid_se(Y, S, Hh, s2)
se = max(norm(Y - Hh*S, 'fro')^2/numel(Y) - s2, 0);
end
function r = rates(Hh, H, s2, evalb, se2)
% [sub, cpe, none] rates of the robust MMSE built from Hh, on the true H
M = size(H, 1);
W = (Hh'*Hh + (s2 + se2)*eye(M)) \ Hh'; G = W*H; g = diag(G); gh = diag(W*Hh);
it = sum(abs(G).^2, 2) - abs(g).^2 + s2*sum(abs(W).^2, 2);
s1 = abs(g).^2 ./ it;
a = (gh(evalb)' * g(evalb)) / (gh(evalb)' * gh(evalb));
sc = abs(a*gh).^2 ./ (it + abs(g - a*gh).^2);
s0 = abs(gh).^2 ./ (it + abs(g - gh).^2);
r = [mean(log2(1 + s1(evalb))), mean(log2(1 + sc(evalb))), mean(log2(1 + s0(evalb)))];
end
function [As, Adp, nu] = moving_phys(Q, bins, df, fc, popt, fD)
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
function C = ici_matrix(M, eps)
k = (0:M-1).' - (0:M-1); x = eps - k;
C = exp(1j*pi*x*(M-1)/M) .* snc(x) ./ snc(x/M);
end
function y = snc(x), y = ones(size(x)); nz = abs(x) > 1e-12; y(nz) = sin(pi*x(nz)) ./ (pi*x(nz)); end
function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
function W = noise(M, n, s2), W = sqrt(s2/2) * (randn(M,n) + 1j*randn(M,n)); end
