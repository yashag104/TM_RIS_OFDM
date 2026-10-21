%WEEK8_FULLBAND_WARM (L4) Measured cost of warm-start updates at the full band.
%   M = 3276, user at 30 km/h (RIS -> user paths rotate with their Doppler; the
%   intra-symbol Doppler ICI, -47 dB at this speed, is neglected). Refresh every 40
%   symbols from Np = 2 pilots, 6 blocks; 6 data symbols in the middle of each block,
%   equalised by CG-MMSE with the estimated operator.
%     cold   full estimate every block (pooled init + 15 ALS iterations), as week4
%     warm   block 0 cold, then 2 ALS iterations from the previous U, V
%     true   CG-MMSE with the true operator (reference)
%   Reported: SER on the central half, wall-clock per update (blocks 1-5).
%   Bar (pre-registered): warm >= 5x faster than cold, SER within 1.2x of cold.
%   S1, Q = 16, r = 8, 16-QAM, 20 dB, 4 trials.

clear; rng(4400);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3; Tsym = 1/df + 2.34e-6;
M = 3276; bins = -M/2:M/2-1; evalb = round(M/4)+1:round(3*M/4);
Q = 16; K = 8; Np = 2; Nd = 6; nit_cg = 25; r = 8; Mq = 16; snr_db = 20;
Nb = 40; B = 6; ntrial = 4;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
L = ceil(M*df*0.2e-6) + 10; [F, lamF] = dpss_basis(bins, df, 0.2e-6, L);
Tst = {fast_coupling(bins, fc, df, base, ones(1,M), ones(1,M), 'state1'), ...
       fast_coupling(bins, fc, df, base, ones(1,M), ones(1,M), 'state2')};
fD = 30/3.6 * fc / 3e8;
err = zeros(ntrial, 3); cnt = 0; tm = zeros(ntrial, 2);
for t = 1:ntrial
    [As, Adp, nu] = moving_channel(Q, bins, df, fc, popt, fD);
    opAt = @(n) fast_coupling(bins, fc, df, seqs, As, sum(Adp .* reshape(exp(1j*2*pi*nu*n*Tsym), 1, 1, []), 3), P);
    op0 = opAt(0);
    cols = evalb(round(linspace(1, numel(evalb), 24)));
    E = zeros(M, numel(cols)); E(sub2ind(size(E), cols, 1:numel(cols))) = 1;
    Pc = mean(sum(abs(op0.H(E)).^2, 1)); s2 = Pc / 10^(snr_db/10);
    Uw = []; Vw = [];
    for b = 0:B-1
        n0 = b*Nb;
        opp = {opAt(n0), opAt(n0 + 1)}; opd = opAt(n0 + Nb/2);
        S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
        Yp = [opp{1}.H(S(:,1)), opp{2}.H(S(:,2))] + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
        d = randi(Mq, M, Nd) - 1; X = qammod_(d, Mq);
        Yd = opd.H(X) + sqrt(s2/2)*(randn(M,Nd) + 1j*randn(M,Nd));
        t0 = tic; [Uc, Vc] = est_timing_lowrank(Yp, S, Tst, F, lamF, s2, r, [], [], 15); tc = toc(t0);
        if b == 0, Uw = Uc; Vw = Vc; tw = tc;
        else, t0 = tic; [Uw, Vw] = est_timing_lowrank(Yp, S, Tst, F, lamF, s2, r, Uw, Vw, 2); tw = toc(t0); end
        if b > 0
            z = {cg_mmse(timing_op(bins, fc, df, base, F, Uc, Vc), Yd, s2, nit_cg), ...
                 cg_mmse(timing_op(bins, fc, df, base, F, Uw, Vw), Yd, s2, nit_cg), ...
                 cg_mmse(opd, Yd, s2, nit_cg)};
            for e = 1:3, err(t, e) = err(t, e) + sum(qamdemod_(z{e}(evalb,:), Mq) ~= d(evalb,:), 'all'); end
            tm(t, :) = tm(t, :) + [tc tw] / (B-1);
            if t == 1, cnt = cnt + numel(evalb)*Nd; end
        end
        fprintf('  trial %d block %d: cold %.1f s, warm %.1f s\n', t, b, tc, tw);
    end
end
ser = sum(err, 1) / (cnt*ntrial);
R = struct('M', M, 'speed_kmh', 30, 'refresh_symbols', Nb, 'snr_db', snr_db, 'qam', Mq, 'ntrial', ntrial, ...
    'names', {{'cold', 'warm', 'true'}}, 'ser', ser, 'ser_trials', err / cnt, ...
    'time_cold_s', mean(tm(:,1)), 'time_warm_s', mean(tm(:,2)), 'speedup', mean(tm(:,1))/mean(tm(:,2)));
fprintf('SER cold %.4f warm %.4f true %.4f | time cold %.1f s warm %.1f s (%.1fx)\n', ser, R.time_cold_s, R.time_warm_s, R.speedup);
fid = fopen(fullfile(out, 'week8_fullband_warm.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

%% ---------------------------------------------------------------------
function [As, Adp, nu] = moving_channel(Q, bins, df, fc, popt, fD)
% As static (phys_channel); RIS -> user hop as separate paths with Doppler nu_l.
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

function [U, V] = est_timing_lowrank(Y, S, Tst, F, lam, s2, r, U, V, nit)
% As week4_fullband_trials, with the iteration count as an argument (warm start if U given).
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
    [Uc, Sg, Wc] = svd([C1 C2]);
    U = Uc(:, 1:r) * sqrt(Sg(1:r,1:r));
    V = {conj(Wc(1:L, 1:r)) * sqrt(Sg(1:r,1:r)), conj(Wc(L+1:end, 1:r)) * sqrt(Sg(1:r,1:r))};
end
wl = repmat(lam(:), r, 1);
scl = mean(abs([U(:); V{1}(:); V{2}(:)]).^2 ./ [wl; wl; wl]);
Du = sqrt(scl*wl); Dv = sqrt(scl*[wl; wl]);
for it = 1:nit
    A = zeros(M*Np, L*r);
    for p = 1:Np
        acc = G{1,p}*V{1} + G{2,p}*V{2};
        for k = 1:r, A((p-1)*M+(1:M), (k-1)*L+(1:L)) = F .* acc(:, k); end
    end
    U = reshape(Du .* (((Du .* (A'*A)) .* Du.' + s2*eye(L*r)) \ (Du .* (A'*y))), L, r);
    FU = F*U; A = zeros(M*Np, 2*L*r);
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
r = size(U, 2); seqv = repmat(base, r, 1);
o1 = fast_coupling(bins, fc, df, seqv, (F*V{1}).', (F*U).', 'state1');
o2 = fast_coupling(bins, fc, df, seqv, (F*V{2}).', (F*U).', 'state2');
op.H = @(x) o1.H(x) + o2.H(x);
op.Ht = @(y) o1.Ht(y) + o2.Ht(y);
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
