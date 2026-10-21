%WEEK8_TRACKING_120 (L3) Track a fast user with the static input side fixed.
%   Fixed BS -> RIS hop; the user moves at 120 km/h (paths rotate with their Doppler,
%   and the intra-symbol Doppler ICI is included in the true channel, as week5_mobility).
%   Pilots on symbols 0 and 1, then data symbols 2..41 (16-QAM at 20 dB, 64-QAM at 30 dB).
%   The input side W is learned beforehand from Nf = 10 frames of other users
%   (est_twoscale 'learn', 2 pilots each).
%     reuse_single  week5 baseline: est_general on the pilots, reused for all symbols
%     reuse_2scale  U-only estimate on the pilots with learned W, reused
%     track         as reuse_2scale, then after each data symbol: MMSE equalise, decide,
%                   update U by LMMSE with prior mean = previous U and prior variance
%                   P0w * sclU * lam, P0w = 2 (2 pi fD Tsym)^2 (set by physics, not tuned)
%   Gap to perfect (rate of MMSE with the estimate vs the true channel at that symbol).
%   Bar (pre-registered): track, gap <= 0.25 bit/s/Hz averaged over symbols 1-40 after
%   the pilots at 20 dB, and <= 0.5 averaged over 1-20 at 30 dB.
%   S1, Q = 16, K = 8, M = 256, physical channel (6 paths, 30 ns), 8 trials.

clear; rng(4500);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3; Tsym = 1/df + 2.34e-6;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8; Nf = 10; Nd = 40;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
fD = 120/3.6 * fc / 3e8; P0w = 2*(2*pi*fD*Tsym)^2;
snrs = [20 30]; qam = [16 64]; ntrial = 8;
names = {'reuse_single', 'reuse_2scale', 'track'};
gap = zeros(ntrial, numel(names), Nd, 2);
for t = 1:ntrial
    [As, geo] = static_hop(Q, bins, df, fc, popt);
    [Adp, nu] = user_paths(geo, bins, df, popt, fD);
    Hat = @(n) channel_at(n, bins, fc, df, seqs, As, Adp, nu, P, M, Tsym);
    Htr = cell(1, Nf);
    for k = 1:Nf, [Ak, ~] = user_paths(geo, bins, df, popt, 0); Htr{k} = coupling_matrix(bins, fc, df, seqs, As, sum(Ak, 3), P, M-1); end
    H0 = Hat(0); Pc = mean(sum(abs(H0(:, evalb)).^2, 1));
    Hn = arrayfun(Hat, 2:Nd+1, 'UniformOutput', false);
    for is = 1:2
        s2 = Pc / 10^(snrs(is)/10); Mq = qam(is);
        Ys = cell(1, Nf); Ss = Ys;
        for k = 1:Nf, Ss{k} = qpsk(M, 2); Ys{k} = Htr{k}*Ss{k} + noise(M, 2, s2); end
        st = est_twoscale('learn', Ys, Ss, T, F, lamF, s2, Q, 8);
        S = qpsk(M, 2); Y = [H0*S(:,1), Hat(1)*S(:,2)] + noise(M, 2, s2);
        Hs = est_general(Y, S, T, ones(Q,1), F, lamF, s2, 8, 20);
        o = est_twoscale('frame', Y, S, T, F, lamF, s2, st.W, st.sclU); H2 = o{1}; U = o{2};
        Htr_ = H2;
        for n = 1:Nd
            H = Hn{n}; rp = mmse_rate(H, H, s2, evalb);
            gap(t, 1, n, is) = rp - mmse_rate(Hs, H, s2, evalb);
            gap(t, 2, n, is) = rp - mmse_rate(H2, H, s2, evalb);
            gap(t, 3, n, is) = rp - mmse_rate(Htr_, H, s2, evalb);
            % receive data symbol n, equalise with the current estimate, decide, update U
            X = reshape(qammod_(randi(Mq, M, 1) - 1, Mq), M, 1); y = H*X + noise(M, 1, s2);
            Xh = qamdec_(((Htr_'*Htr_ + s2*eye(M)) \ Htr_') * y, Mq);
            o = est_twoscale('frame', y, Xh, T, F, lamF, s2, st.W, st.sclU, U, P0w); Htr_ = o{1}; U = o{2};
        end
    end
    fprintf('trial %d done\n', t);
end
g = squeeze(mean(gap, 1));                                   % names x Nd x snr
R = struct('speed_kmh', 120, 'fD_hz', fD, 'P0w', P0w, 'names', {names}, 'snr_db', snrs, 'qam', qam, ...
    'gap_bits', g, 'avg_1_40_20dB', mean(g(:, 1:40, 1), 2).', 'avg_1_20_30dB', mean(g(:, 1:20, 2), 2).');
for e = 1:numel(names)
    fprintf('%-13s avg gap 1-40 @20 %.2f | 1-20 @30 %.2f | at symbols 1/10/20/40 @20: %s\n', names{e}, ...
        R.avg_1_40_20dB(e), R.avg_1_20_30dB(e), sprintf('%.2f ', g(e, [1 10 20 40], 1)));
end
fid = fopen(fullfile(out, 'week8_tracking_120.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

%% ---------------------------------------------------------------------
function [As, geo] = static_hop(Q, bins, df, fc, popt)
[As, ~] = phys_channel(Q, bins, df, fc, popt);
lam = 3e8 / fc; n = round(sqrt(Q));
[gx, gz] = meshgrid(0:n-1, 0:n-1);
geo = struct('r', [gx(:), zeros(Q,1), gz(:)] * lam/2, 'k', 2*pi/lam, 'Q', Q);
end
function [Adp, nu] = user_paths(geo, bins, df, popt, fD)
P = popt.P; Kr = 10^(popt.K_db/10);
tau = [0, sort(-popt.rms*log(rand(1, P)))];
pw = [Kr/(Kr+1), (1/(Kr+1)) * exp(-tau(2:end)/popt.rms) / sum(exp(-tau(2:end)/popt.rms))];
g = sqrt(pw) .* exp(1j*2*pi*rand(1, P+1));
u = randn(P+1, 3); u(:,2) = abs(u(:,2)); u = u ./ vecnorm(u, 2, 2);
ph = exp(-1j*geo.k*(geo.r*u.'));
Adp = zeros(geo.Q, numel(bins), P+1);
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
function x = qammod_(d, Mq)
m = sqrt(Mq); lv = -(m-1):2:(m-1);
x = (lv(mod(d, m) + 1) + 1j*lv(floor(d / m) + 1)) / sqrt(2*(Mq-1)/3);
end
function x = qamdec_(z, Mq)
m = sqrt(Mq); sc = sqrt(2*(Mq-1)/3); lv = -(m-1):2:(m-1);
q = @(v) lv(min(max(round((v + m - 1)/2), 0), m-1) + 1);
x = reshape((q(real(z)*sc) + 1j*q(imag(z)*sc)) / sc, size(z));
end
function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
