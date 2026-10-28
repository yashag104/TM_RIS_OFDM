%WEEK13_FULLBAND Review round 2 (pre-registered 2026-10-28, Step 13b): full band, M = 3276 (98 MHz), where the element
%   is dispersive (phase drift ~47 deg across the band). Matrix-free throughout.
%   Receivers (CG-MMSE on the estimate; SER on the central half; NMSE on 48 sampled central columns):
%     onetap, true      references (one-tap with the true diagonal; CG with the true operator)
%     ours              timing-only low-rank estimator (Algorithm 1, no fallback; est_timing_lowrank as week11_fullband30)
%     flatcal           ideal-element LMMSE with per-state flat complex scalars Gamma_s(fc) instead of +-1
%     flatest           the same with Gamma_2/Gamma_1 estimated from the pilots (element-free linear stage)
%     genie             LMMSE given the true dispersive Gamma_s(f) (linear, full L x L coefficients)
%     ideal             LMMSE with the +-1 element
%     bem               CE-BEM, |h| <= 8, two-hop DPSS prior (as week11_fullband30)
%   Np = 1, 2, 4 pilot symbols on the same channels; 16-QAM at 20 dB and 64-QAM at 28 dB; 30 new channels (rng 13500+).
%   Shared-path channel (LoS + 5 paths, K = 10 dB, 30 ns rms). Saves after every channel.
%   Set trials (and outname) to run a subset of channels in a separate process; draws are identical per channel.

if ~exist('ntrial', 'var'), ntrial = 30; end
if ~exist('outname', 'var'), outname = 'week13_fullband'; end
if ~exist('trials', 'var'), trials = 1:ntrial; end              % channel subset (seeds are per channel)
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3;
M = 3276; bins = -M/2:M/2-1; evalb = round(M/4)+1:round(3*M/4); f = fc + bins(:)*df;
Q = 16; K = 8; D = 0.75; Nd = 6; nit_cg = 25; r = 8; Qb = 8; Nps = [1 2 4];
base = [ones(1, round(D*K)), 2*ones(1, K - round(D*K))]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
L = ceil(M*df*0.2e-6) + 10; [F, lamF] = dpss_basis(bins, df, 0.2e-6, L);
L2 = ceil(M*df*0.4e-6) + 12; [F2, lam2] = dpss_basis(bins, df, 0.4e-6, L2);
Tst = {fast_coupling(bins, fc, df, base, ones(1,M), ones(1,M), 'state1'), ...
       fast_coupling(bins, fc, df, base, ones(1,M), ones(1,M), 'state2')};
gam_true = {element_response(f, 1, P), element_response(f, 2, P)};
gam_flat = {element_response(fc, 1, P), element_response(fc, 2, P)};
cfg = {20, 16; 28, 64};
names = {'onetap', 'true', 'ours', 'flatcal', 'flatest', 'genie', 'ideal', 'bem'};
cols = evalb(round(linspace(1, numel(evalb), 48))); E = zeros(M, numel(cols)); E(sub2ind(size(E), cols, 1:numel(cols))) = 1;
fp = fullfile(out, [outname '.json']);
R = struct('M', M, 'L', L, 'L2', L2, 'rank', r, 'Np', Nps, 'names', {names}, 'cfg', {cfg});
ser = nan(numel(trials), numel(names), numel(Nps), size(cfg, 1)); nm = ser;
for it = 1:numel(trials)
    t = trials(it); rng(13500 + t);
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    op = fast_coupling(bins, fc, df, seqs, As, Ad, P); HE = op.H(E);
    Pc = mean(sum(abs(HE).^2, 1));
    S4 = (sign(randn(M,4)) + 1j*sign(randn(M,4))) / sqrt(2); W4 = randn(M,4) + 1j*randn(M,4);
    Wd = randn(M,Nd) + 1j*randn(M,Nd); dd = randi(64, M, Nd) - 1;
    for ic = 1:size(cfg, 1)
        s2 = Pc / 10^(cfg{ic,1}/10); Mq = cfg{ic,2}; d = mod(dd, Mq); X = qammod_(d, Mq);
        Yd = op.H(X) + sqrt(s2/2)*Wd; Y4 = op.H(S4) + sqrt(s2/2)*W4;
        for ip = 1:numel(Nps)
            Np = Nps(ip); S = S4(:, 1:Np); Yp = Y4(:, 1:Np);
            ops = cell(1, numel(names));
            ops{2} = op;
            [U, V] = est_timing_lowrank(Yp, S, Tst, F, lamF, s2, r, [], []);
            ops{3} = timing_op(bins, fc, df, base, F, U, V);
            ops{4} = lin_op(Yp, S, Tst, gam_flat, F, lamF, s2, bins, fc, df, base);
            [C1, C2] = lin2(Yp, S, Tst, F, lamF, s2);
            rho = (C1(:)' * C2(:)) / max(real(C1(:)' * C1(:)), eps);
            ops{5} = lin_op(Yp, S, Tst, {1, rho}, F, lamF, s2, bins, fc, df, base);
            ops{6} = lin_op(Yp, S, Tst, gam_true, F, lamF, s2, bins, fc, df, base);
            ops{7} = lin_op(Yp, S, Tst, {1, -1}, F, lamF, s2, bins, fc, df, base);
            c = est_bem(Yp, S, F2, lam2, s2, -Qb:Qb); ops{8} = bem_op(F2, c, -Qb:Qb);
            for e = 1:numel(names)
                if e == 1, z = Yd ./ onetap_diag(op, M); else, z = cg_mmse(ops{e}, Yd, s2, nit_cg); end
                ser(it, e, ip, ic) = mean(qamdemod_(z(evalb,:), Mq) ~= d(evalb,:), 'all');
                if e >= 2, nm(it, e, ip, ic) = norm(ops{e}.H(E) - HE, 'fro')^2 / norm(HE, 'fro')^2; end
            end
            fprintf('ch %2d | %d dB %d-QAM Np %d | SER %s\n', t, cfg{ic,1}, Mq, Np, sprintf('%.3g ', ser(it, :, ip, ic)));
        end
    end
    R.trials = trials(1:it); R.ser = ser(1:it,:,:,:); R.nmse = nm(1:it,:,:,:);
    fid = fopen(fp, 'w'); fprintf(fid, '%s', jsonencode(R)); fclose(fid);
end
fprintf('done\n');

%% ---------------------------------------------------------------------
function op = lin_op(Y, S, Tst, gam, F, lam, s2, bins, fc, df, base)
% LMMSE of H = sum_s T_s o (F C F^T) diag(gam_s) (gam_s scalar or M x 1), then its matrix-free operator
[M, Np] = size(S); L = size(F, 2); y = Y(:);
Phi = zeros(M*Np, L*L);
for p = 1:Np
    Gs = Tst{1}.H(gam{1} .* (S(:,p) .* F)) + Tst{2}.H(gam{2} .* (S(:,p) .* F));
    for i = 1:L, Phi((p-1)*M+(1:M), (i-1)*L+(1:L)) = F(:, i) .* Gs; end
end
w = kron(lam(:), lam(:));
sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
Dw = sqrt(sc*w);
c = Dw .* (((Dw .* (Phi'*Phi)) .* Dw.' + s2*eye(numel(w))) \ (Dw .* (Phi'*y)));
C = reshape(c, L, L).';
[Uc, Sc, Vc] = svd(C); rk = L; sq = sqrt(Sc(1:rk,1:rk));
FU = F*Uc(:,1:rk)*sq; FV = F*conj(Vc(:,1:rk))*sq; seqv = repmat(base, rk, 1);
o1 = fast_coupling(bins, fc, df, seqv, (gam{1} .* FV).', FU.', 'state1');
o2 = fast_coupling(bins, fc, df, seqv, (gam{2} .* FV).', FU.', 'state2');
op.H = @(x) o1.H(x) + o2.H(x);
op.Ht = @(y) o1.Ht(y) + o2.Ht(y);
end
function [C1, C2] = lin2(Y, S, Tst, F, lam, s2)
% element-free linear stage (two coefficient matrices), used only for the ratio of flatest
[M, Np] = size(S); L = size(F, 2); y = Y(:);
Phi = zeros(M*Np, 2*L*L);
for s = 1:2, for p = 1:Np
    Gs = Tst{s}.H(S(:,p) .* F);
    for i = 1:L, Phi((p-1)*M+(1:M), (s-1)*L*L+(i-1)*L+(1:L)) = F(:, i) .* Gs; end
end, end
pl = kron(lam(:), lam(:)); w = [pl; pl];
sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
Dw = sqrt(sc*w);
c = Dw .* (((Dw .* (Phi'*Phi)) .* Dw.' + s2*eye(numel(w))) \ (Dw .* (Phi'*y)));
C1 = reshape(c(1:L*L), L, L).'; C2 = reshape(c(L*L+1:end), L, L).';
end
function [U, V] = est_timing_lowrank(Y, S, Tst, F, lam, s2, r, U, V)
% H = sum_s T_s o ((F U)(F V{s}).'). Init (if U empty) from the per-state linear
% LMMSE (2 L^2 unknowns), then shared-output-factor ALS with DPSS priors.
[M, Np] = size(S); L = size(F, 2); y = Y(:);
G = cell(2, Np);
for s = 1:2, for p = 1:Np, G{s,p} = Tst{s}.H(S(:,p) .* F); end, end
if isempty(U)
    pl = kron(lam(:), lam(:)); Phi = zeros(M*Np, 2*L*L);
    for s = 1:2
        for p = 1:Np
            for i = 1:L, Phi((p-1)*M+(1:M), (s-1)*L*L+(i-1)*L+(1:L)) = F(:, i) .* G{s,p}; end
        end
    end
    w = [pl; pl];
    sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
    Dw = sqrt(sc*w);
    c = Dw .* (((Dw .* (Phi'*Phi)) .* Dw.' + s2*eye(numel(w))) \ (Dw .* (Phi'*y)));
    C1 = reshape(c(1:L*L), L, L).'; C2 = reshape(c(L*L+1:end), L, L).';
    [Uc, Sg, Wc] = svd([C1 C2]);                       % shared output side
    U = Uc(:, 1:r) * sqrt(Sg(1:r,1:r));
    V = {conj(Wc(1:L, 1:r)) * sqrt(Sg(1:r,1:r)), conj(Wc(L+1:end, 1:r)) * sqrt(Sg(1:r,1:r))};
end
wl = repmat(lam(:), r, 1);
scl = mean(abs([U(:); V{1}(:); V{2}(:)]).^2 ./ [wl; wl; wl]);
Du = sqrt(scl*wl); Dv = sqrt(scl*[wl; wl]);
for it = 1:15
    A = zeros(M*Np, L*r);                               % U-step
    for p = 1:Np
        acc = G{1,p}*V{1} + G{2,p}*V{2};                % M x r
        for k = 1:r, A((p-1)*M+(1:M), (k-1)*L+(1:L)) = F .* acc(:, k); end
    end
    U = reshape(Du .* (((Du .* (A'*A)) .* Du.' + s2*eye(L*r)) \ (Du .* (A'*y))), L, r);
    FU = F*U; A = zeros(M*Np, 2*L*r);                  % V-step (both states jointly)
    for s = 1:2
        for p = 1:Np
            for k = 1:r, A((p-1)*M+(1:M), (s-1)*L*r+(k-1)*L+(1:L)) = FU(:, k) .* G{s,p}; end
        end
    end
    v = Dv .* (((Dv .* (A'*A)) .* Dv.' + s2*eye(2*L*r)) \ (Dv .* (A'*y)));
    V = {reshape(v(1:L*r), L, r), reshape(v(L*r+1:end), L, r)};
end
end

function op = timing_op(bins, fc, df, base, F, U, V)
% Operator of the timing-only estimate: per state, r virtual elements.
r = size(U, 2); seqv = repmat(base, r, 1);
o1 = fast_coupling(bins, fc, df, seqv, (F*V{1}).', (F*U).', 'state1');
o2 = fast_coupling(bins, fc, df, seqv, (F*V{2}).', (F*U).', 'state2');
op.H = @(x) o1.H(x) + o2.H(x);
op.Ht = @(y) o1.Ht(y) + o2.Ht(y);
end

function c = est_bem(Y, S, F2, lam, s2, hs)
[M, Np] = size(S); L = size(F2, 2); Nh = numel(hs); y = Y(:);
Phi = zeros(M*Np, Nh*L);
for k = 1:Nh, for p = 1:Np
    Phi((p-1)*M+(1:M), (k-1)*L+(1:L)) = shrow(S(:,p) .* F2, hs(k));
end, end
w = repmat(lam(:), Nh, 1);
sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
Dw = sqrt(sc*w);
c = Dw .* (((Dw .* (Phi'*Phi)) .* Dw.' + s2*eye(numel(w))) \ (Dw .* (Phi'*y)));
c = reshape(c, L, Nh);
end

function op = bem_op(F2, c, hs)
g = F2 * c;                                             % M x Nh, per-harmonic gains
op.H = @(x) bem_apply(g, hs, x, false);
op.Ht = @(y) bem_apply(g, hs, y, true);
end

function out = bem_apply(g, hs, x, adj)
out = zeros(size(x));
for k = 1:numel(hs)
    if ~adj, out = out + shrow(g(:,k) .* x, hs(k));
    else, out = out + conj(g(:,k)) .* shrow(x, -hs(k)); end
end
end

function Z = shrow(X, h)
Z = zeros(size(X)); M = size(X, 1);
if h >= 0, Z(1+h:M, :) = X(1:M-h, :); else, Z(1:M+h, :) = X(1-h:M, :); end
end

function dg = onetap_diag(op, M)
dg = zeros(M, 1); blk = 256;
for s = 1:blk:M
    idx = s:min(M, s+blk-1);
    E = zeros(M, numel(idx)); E(sub2ind(size(E), idx, 1:numel(idx))) = 1;
    Hc = op.H(E); dg(idx) = Hc(sub2ind(size(Hc), idx, 1:numel(idx)));
end
end

function z = cg_mmse(op, Y, s2, nit)
b = op.Ht(Y); z = zeros(size(b)); r = b; p = r; rs = sum(conj(r).*r, 1);
for it = 1:nit
    Ap = op.Ht(op.H(p)) + s2*p;
    a = rs ./ sum(conj(p).*Ap, 1);
    z = z + a.*p; r = r - a.*Ap;
    rn = sum(conj(r).*r, 1);
    p = r + (rn./rs).*p; rs = rn;
end
end

function x = qammod_(d, Mq)
m = sqrt(Mq); lv = -(m-1):2:(m-1);
x = (lv(mod(d, m) + 1) + 1j*lv(floor(d / m) + 1)) / sqrt(2*(Mq-1)/3);
end
function d = qamdemod_(z, Mq)
m = sqrt(Mq); sc = sqrt(2*(Mq-1)/3);
ir = min(max(round((real(z)*sc + m - 1)/2), 0), m-1);
iq = min(max(round((imag(z)*sc + m - 1)/2), 0), m-1);
d = ir + m*iq;
end
function x = qamdec_(z, Mq)
x = qammod_(qamdemod_(z, Mq), Mq);
end
