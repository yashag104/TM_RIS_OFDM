%WEEK4_FULLBAND_TRIALS Full n78 band, 10 trials, 95% CIs, timing-only vs CE-BEM.
%   M = 3276, physical channel, S1 (identical sequences), matrix-free throughout.
%   Receivers (SER on the central half of the band, Np = 2 pilots, Nd = 6 data):
%     onetap        conventional one-tap with the true diagonal (favourable to it)
%     true          CG-MMSE with the true operator
%     timing_p      timing-only low-rank estimator, pilots only
%     timing_dd     + 1 decision-directed pass
%     bem_p, bem_dd CE-BEM (Q = 8, two-hop DPSS, LMMSE prior), pilots only / + 1 DD pass
%   Timing-only model (from Lemma 1): the output hop is common to both states,
%     H = sum_s T_s o ((F U)(F V_s).'),  U, V_1, V_2 in C^{L x r}, r = 8.

clear; rng(2200);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3;
M = 3276; bins = -M/2:M/2-1; evalb = round(M/4)+1:round(3*M/4);
Q = 16; K = 8; D = 0.75; Np = 2; Nd = 6; nit_cg = 25; r = 8; Qb = 8;
base = [ones(1, round(D*K)), 2*ones(1, K - round(D*K))];
seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
L = ceil(M*df*0.2e-6) + 10; [F, lamF] = dpss_basis(bins, df, 0.2e-6, L);
L2 = ceil(M*df*0.4e-6) + 12; [F2, lam2] = dpss_basis(bins, df, 0.4e-6, L2);
Tst = {fast_coupling(bins, fc, df, base, ones(1,M), ones(1,M), 'state1'), ...
       fast_coupling(bins, fc, df, base, ones(1,M), ones(1,M), 'state2')};
cfg = {20, 16; 28, 64}; ntrial = 10;
names = {'onetap', 'true', 'timing_p', 'timing_dd', 'bem_p', 'bem_dd'};
R = struct('M', M, 'L', L, 'L2', L2, 'rank', r, 'bemQ', Qb, 'ntrial', ntrial, 'names', {names});
for ic = 1:size(cfg, 1)
    snr_db = cfg{ic,1}; Mq = cfg{ic,2};
    ser = zeros(ntrial, numel(names));
    for t = 1:ntrial
        [As, Ad] = phys_channel(Q, bins, df, fc, popt);
        op = fast_coupling(bins, fc, df, seqs, As, Ad, P);
        cols = evalb(round(linspace(1, numel(evalb), 24)));
        E = zeros(M, numel(cols)); E(sub2ind(size(E), cols, 1:numel(cols))) = 1;
        Pc = mean(sum(abs(op.H(E)).^2, 1)); s2 = Pc / 10^(snr_db/10);
        S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
        Yp = op.H(S) + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
        d = randi(Mq, M, Nd) - 1; X = qammod_(d, Mq);
        Yd = op.H(X) + sqrt(s2/2)*(randn(M,Nd) + 1j*randn(M,Nd));
        z = cell(1, numel(names));
        z{1} = Yd ./ onetap_diag(op, M);
        z{2} = cg_mmse(op, Yd, s2, nit_cg);
        % timing-only
        [U, V] = est_timing_lowrank(Yp, S, Tst, F, lamF, s2, r, [], []);
        opt = timing_op(bins, fc, df, base, F, U, V);
        z{3} = cg_mmse(opt, Yd, s2, nit_cg);
        Xh = qamdec_(z{3}, Mq);
        [U, V] = est_timing_lowrank([Yp Yd], [S Xh], Tst, F, lamF, s2, r, U, V);
        z{4} = cg_mmse(timing_op(bins, fc, df, base, F, U, V), Yd, s2, nit_cg);
        % CE-BEM
        c = est_bem(Yp, S, F2, lam2, s2, -Qb:Qb);
        opb = bem_op(F2, c, -Qb:Qb);
        z{5} = cg_mmse(opb, Yd, s2, nit_cg);
        Xh = qamdec_(z{5}, Mq);
        c = est_bem([Yp Yd], [S Xh], F2, lam2, s2, -Qb:Qb);
        z{6} = cg_mmse(bem_op(F2, c, -Qb:Qb), Yd, s2, nit_cg);
        for e = 1:numel(names)
            ser(t, e) = mean(qamdemod_(z{e}(evalb,:), Mq) ~= d(evalb,:), 'all');
        end
        fprintf('  snr %d trial %d: %s\n', snr_db, t, sprintf('%.3g ', ser(t,:)));
    end
    key = sprintf('snr%d_qam%d', snr_db, Mq);
    R.(key) = struct('ser_mean', mean(ser, 1), 'ser_ci95', 1.96*std(ser, 0, 1)/sqrt(ntrial), 'ser_trials', ser);
    fprintf('%s SER mean (95%% CI):\n', key);
    for e = 1:numel(names), fprintf('  %-10s %.4f +- %.4f\n', names{e}, mean(ser(:,e)), 1.96*std(ser(:,e))/sqrt(ntrial)); end
end
fid = fopen(fullfile(out, 'week4_fullband_trials.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

%% ---------------------------------------------------------------------
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
