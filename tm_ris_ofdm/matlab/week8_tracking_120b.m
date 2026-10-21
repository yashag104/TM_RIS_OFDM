%WEEK8_TRACKING_120B (L3b) Kalman tracking of U at 120 km/h, input side learned.
%   Setup as week8_tracking_120 (same seed, channels, symbols). Differences:
%     - W and the U prior Ru from learn_input_side (stacked-SVD init, learned prior)
%     - a proper Kalman filter on u = vec(U):
%         prior at the pilots: u ~ CN(0, C0), C0 = kron(Ru, diag(lam))
%         predict:  P <- P + Qp,  Qp = P0w C0,  P0w = 2 (2 pi fD Tsym)^2 (physics, not tuned)
%         update:   y = A(x_hat) u + n,  noise var = max(s2, |y - A u_pred|^2 / M)
%   Receivers: reuse_single (week5 baseline), reuse_2scale (learned W, pilots only), kalman.
%   Bar (pre-registered, as L3): kalman, gap <= 0.25 averaged over symbols 1-40 at 20 dB
%   and <= 0.5 over 1-20 at 30 dB.

clear; rng(4500);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3; Tsym = 1/df + 2.34e-6;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8; Nf = 10; Nd = 40;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16); L = size(F, 2);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
fD = 120/3.6 * fc / 3e8; P0w = 2*(2*pi*fD*Tsym)^2;
snrs = [20 30]; qam = [16 64]; ntrial = 8;
names = {'reuse_single', 'reuse_2scale', 'kalman'};
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
        st = learn_input_side(Ys, Ss, T, F, lamF, s2, Q, 8);
        S = qpsk(M, 2); Y = [H0*S(:,1), Hat(1)*S(:,2)] + noise(M, 2, s2);
        Hs = est_general(Y, S, T, ones(Q,1), F, lamF, s2, 8, 20);
        r = size(st.W{1}, 2);
        C0 = kron(st.Ru, diag(lamF)); C0 = (C0 + C0')/2; Qp = P0w*C0;
        A = amat(S, T, F, st.W);                                   % pilots: Kalman init
        [u, Pk] = kupdate(zeros(L*r, 1), C0, A, Y(:), s2);
        H2 = hmat(u, F, T, st.W);
        for n = 1:Nd
            H = Hn{n}; rp = mmse_rate(H, H, s2, evalb);
            Pk = Pk + Qp; Hk = hmat(u, F, T, st.W);               % predict
            gap(t, 1, n, is) = rp - mmse_rate(Hs, H, s2, evalb);
            gap(t, 2, n, is) = rp - mmse_rate(H2, H, s2, evalb);
            gap(t, 3, n, is) = rp - mmse_rate(Hk, H, s2, evalb);
            X = reshape(qammod_(randi(Mq, M, 1) - 1, Mq), M, 1); y = H*X + noise(M, 1, s2);
            Xh = qamdec_(((Hk'*Hk + s2*eye(M)) \ Hk') * y, Mq);
            A = amat(Xh, T, F, st.W);
            sn = max(s2, mean(abs(y - A*u).^2));
            [u, Pk] = kupdate(u, Pk, A, y, sn);                    % update
        end
    end
    fprintf('trial %d done\n', t);
end
g = squeeze(mean(gap, 1));
R = struct('speed_kmh', 120, 'fD_hz', fD, 'P0w', P0w, 'names', {names}, 'snr_db', snrs, 'qam', qam, ...
    'gap_bits', g, 'avg_1_40_20dB', mean(g(:, 1:40, 1), 2).', 'avg_1_20_30dB', mean(g(:, 1:20, 2), 2).');
for e = 1:numel(names)
    fprintf('%-13s avg gap 1-40 @20 %.2f | 1-20 @30 %.2f | at symbols 1/10/20/40 @20: %s| @30: %s\n', names{e}, ...
        R.avg_1_40_20dB(e), R.avg_1_20_30dB(e), sprintf('%.2f ', g(e, [1 10 20 40], 1)), sprintf('%.2f ', g(e, [1 10 20 40], 2)));
end
fid = fopen(fullfile(out, 'week8_tracking_120b.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

%% ---------------------------------------------------------------------
function A = amat(S, T, F, W)
[M, Np] = size(S); L = size(F, 2); r = size(W{1}, 2); A = zeros(M*Np, L*r);
for p = 1:Np
    g = T{1} * (S(:,p) .* W{1}) + T{2} * (S(:,p) .* W{2});
    for a = 1:r, A((p-1)*M+(1:M), (a-1)*L+(1:L)) = F .* g(:, a); end
end
end
function [u, P] = kupdate(u, P, A, y, sn)
S = A*P*A' + sn*eye(size(A, 1)); S = (S + S')/2;
Kg = (P*A') / S;
u = u + Kg*(y - A*u); P = P - Kg*A*P; P = (P + P')/2;
end
function H = hmat(u, F, T, W)
L = size(F, 2); U = reshape(u, L, []); FU = F*U;
H = T{1} .* (FU*W{1}.') + T{2} .* (FU*W{2}.');
end
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
