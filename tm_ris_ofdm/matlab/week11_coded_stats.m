%WEEK11_CODED_STATS Review F2: coded BLER with >= 2000 blocks per point near the waterfall (pre-registered).
%   As week6_coded, with 200 channels x 10 blocks per SNR, per-channel error counts saved for a clustered
%   bootstrap of the 10%-BLER crossing, and CE-BEM also given the robust equaliser (C1).
%   (Original header follows.)
%WEEK6_CODED (#6) Coded performance: block error rate and throughput.
%   Code: rate-1/2 convolutional, K = 7, generators (133, 171) octal, soft-decision
%   Viterbi, random interleaver, 16-QAM Gray, max-log LLRs. Each receiver computes its
%   own LLR reliability from its own channel estimate (mismatch is not hidden):
%     MMSE receivers: unbiased output x = z / mu, noise variance (1 - mu) / mu, with
%       mu = diag(W Hh) from the ESTIMATE; one-tap: noise variance s2 / |g|^2 (it does
%       not know the harmonic interference, as a conventional receiver would not).
%   Packet: the central 128 subcarriers x 2 OFDM symbols x 4 bits = 1024 coded bits
%   (506 information bits + 6 tail bits). Physical channel, S1, Q = 16, M = 256,
%   Np = 2 pilot symbols per channel, 30 channels x 4 packets per SNR point.

clear; rng(11200);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8; Np = 2; Nd = 2;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
L2 = ceil(M*df*0.6e-6) + 12; [F2, lam2] = dpss_basis(bins, df, 0.6e-6, L2);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
trel = make_trellis();
ninfo = 506; ncoded = 2*(ninfo + 6);                   % 1024 = 128 x 2 x 4
snrs = [6 7 8 9 10 11 12 14]; nchan = 200; npk = 10;
names = {'onetap', 'bem12', 'bem40', 'bem40_robust', 'ours', 'perfect'};
blk = zeros(numel(names), numel(snrs)); npkt = 0; blkc = zeros(numel(names), numel(snrs), nchan);
for c = 1:nchan
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
    Pc = mean(sum(abs(H(:, evalb)).^2, 1));
    for is = 1:numel(snrs)
        s2 = Pc / 10^(snrs(is)/10);
        S = qpsk(M, Np); Yp = H*S + noise(M, Np, s2);
        Hb = est_bem_cg(Yp, S, F2, lam2, s2, -12:12); Hb40 = est_bem_cg(Yp, S, F2, lam2, s2, -40:40);
        se40 = max(norm(Yp - Hb40*S, 'fro')^2/(M*Np) - s2, 0);
        Ho = est_general(Yp, S, T, ones(Q,1), F, lamF, s2, 8, 20);
        g1 = mean(Yp ./ S, 2);
        for k = 1:npk
            info = randi([0 1], ninfo, 1);
            cw = conv_enc([info; zeros(6,1)], trel);
            perm = randperm(ncoded); cwi = cw(perm);
            X = zeros(M, Nd); X(evalb, :) = reshape(qam16_map(cwi), numel(evalb), Nd);
            X(setdiff(1:M, evalb), :) = reshape(qam16_map(randi([0 1], 4*(M-numel(evalb))*Nd, 1)), [], Nd);
            Yd = H*X + noise(M, Nd, s2);
            for e = 1:numel(names)
                switch names{e}
                    case 'onetap', xh = Yd ./ g1; nv = s2 ./ abs(g1).^2 .* ones(1, Nd);
                    case 'bem12',  [xh, nv] = mmse_out(Hb, Yd, s2);
                    case 'bem40',  [xh, nv] = mmse_out(Hb40, Yd, s2);
                    case 'bem40_robust', [xh, nv] = mmse_out(Hb40, Yd, s2 + se40);
                    case 'ours',   [xh, nv] = mmse_out(Ho, Yd, s2);
                    case 'perfect',[xh, nv] = mmse_out(H, Yd, s2);
                end
                llr = qam16_llr(xh(evalb, :), nv(evalb, :));
                dec = viterbi_soft(deint(llr, perm), trel, ninfo + 6);
                err = any(dec(1:ninfo) ~= info); blk(e, is) = blk(e, is) + err; blkc(e, is, c) = blkc(e, is, c) + err;
            end
        end
    end
    npkt = npkt + npk;
    if mod(c, 10) == 0, fprintf('channel %d done\n', c); end
end
bler = blk / npkt;
R = struct('snr_db', snrs, 'names', {names}, 'bler', bler, 'packets_per_point', npkt, 'blocks_per_channel', npk, 'errors_per_channel', blkc, ...
    'code', 'conv K=7 (133,171) rate 1/2, soft Viterbi', 'info_bits', ninfo);
R.throughput_bits_per_subcarrier = (1 - bler) * (ninfo / ncoded) * 4;  % 16-QAM, rate ~1/2
for e = 1:numel(names), fprintf('%-8s BLER %s\n', names{e}, sprintf('%.3f ', bler(e,:))); end
fid = fopen(fullfile(out, 'week11_coded_stats.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

%% ---------------------------------------------------------------------
function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
function W = noise(M, n, s2), W = sqrt(s2/2) * (randn(M,n) + 1j*randn(M,n)); end

function [xh, nv] = mmse_out(Hh, Y, s2)
% Unbiased MMSE output and its noise variance, both from the (estimated) Hh.
M = size(Hh, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh';
mu = real(diag(W*Hh)); mu = min(max(mu, 1e-6), 1 - 1e-6);
xh = (W*Y) ./ mu; nv = ((1 - mu) ./ mu) .* ones(1, size(Y, 2));
end

function x = qam16_map(b)
% Gray 16-QAM, unit average energy; b: 4 bits per symbol (I1 I2 Q1 Q2).
b = reshape(b, 4, []);
gi = 2*b(1,:) + b(2,:); gq = 2*b(3,:) + b(4,:);
lv = [-3 -1 3 1];                                     % Gray: 00->-3 01->-1 11->1 10->3
x = (lv(gi + 1) + 1j*lv(gq + 1)).' / sqrt(10);
end

function llr = qam16_llr(xh, nv)
% Max-log LLRs (positive favours bit 0) for the Gray map above, per axis.
xh = xh(:); nv = nv(:); a = 1/sqrt(10);
pts = [-3 -1 1 3] * a; bits1 = [0 0 1 1]; bits2 = [0 1 1 0];
llr = zeros(4, numel(xh));
for ax = 1:2
    if ax == 1, v = real(xh); else, v = imag(xh); end
    d = (v - pts).^2 ./ nv;                           % N x 4
    l1 = min(d(:, bits1 == 1), [], 2) - min(d(:, bits1 == 0), [], 2);
    l2 = min(d(:, bits2 == 1), [], 2) - min(d(:, bits2 == 0), [], 2);
    llr(2*ax-1, :) = l1.'; llr(2*ax, :) = l2.';
end
llr = llr(:);
end

function y = deint(llr, perm), y = zeros(size(llr)); y(perm) = llr; end

function t = make_trellis()
g = [bin2vec(oct2dec_('133')); bin2vec(oct2dec_('171'))];   % 2 x 7, MSB = current bit
ns = 64; t.next = zeros(ns, 2); t.out = zeros(ns, 2, 2);
for s = 0:ns-1
    for u = 0:1
        reg = [u, bitget(s, 6:-1:1)];                         % current input + 6 memory bits
        t.out(s+1, u+1, :) = mod(g * reg.', 2);
        t.next(s+1, u+1) = bitshift(s, -1) + u*32;
    end
end
end
function v = bin2vec(d), v = bitget(d, 7:-1:1); end
function d = oct2dec_(s), d = 0; for ch = s, d = d*8 + (ch - '0'); end, end

function c = conv_enc(u, t)
s = 0; c = zeros(2*numel(u), 1);
for i = 1:numel(u)
    c(2*i-1:2*i) = squeeze(t.out(s+1, u(i)+1, :));
    s = t.next(s+1, u(i)+1);
end
end

function u = viterbi_soft(llr, t, n)
% Vectorised soft Viterbi: maximise sum (1 - 2c) * llr; start and end in state 0.
% Next state ns = (s >> 1) + 32 u, so ns has predecessors 2*mod(ns,32) and +1, input floor(ns/32).
ns = 64; nsv = (0:ns-1).';
pA = 2*mod(nsv, 32) + 1; pB = pA + 1; uin = floor(nsv/32) + 1;     % 1-based
oA = [t.out(sub2ind([ns 2 2], pA, uin, ones(ns,1))), t.out(sub2ind([ns 2 2], pA, uin, 2*ones(ns,1)))];
oB = [t.out(sub2ind([ns 2 2], pB, uin, ones(ns,1))), t.out(sub2ind([ns 2 2], pB, uin, 2*ones(ns,1)))];
sA = 1 - 2*oA; sB = 1 - 2*oB;
pm = -inf(ns, 1); pm(1) = 0; choice = false(ns, n);
for i = 1:n
    l = llr(2*i-1:2*i);
    mA = pm(pA) + sA * l; mB = pm(pB) + sB * l;
    choice(:, i) = mB > mA; pm = max(mA, mB);
end
u = zeros(n, 1); st = 1;
for i = n:-1:1
    u(i) = uin(st) - 1;
    if choice(st, i), st = pB(st); else, st = pA(st); end
end
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
function Z = shr(X, h)
Z = zeros(size(X)); M = size(X, 1);
if h >= 0, Z(1+h:M, :) = X(1:M-h, :); else, Z(1:M+h, :) = X(1-h:M, :); end
end
