%WEEK4_FULLBAND_S2_SHIFT Per-element cyclic shifts at full band via the shift structure.
%   For a cyclic shift k of the base sequence, t_{q,s}[h] = t_{base,s}[h] e^{-j2pi h k/K};
%   with h = mbar - m this is diag(e^{-j2pi mbar k/K}) T_s diag(e^{+j2pi m k/K}). So
%     H = sum_s T_s o sum_k (D_k F U) Lam_k (D_k F V_s)^H-conj  ... i.e.
%     H x = sum_s sum_k D_k [ (F U) Lam_k (F V_s).' o T_s ] D_k^* x
%   with U (output, shared), V_1, V_2 (input per state), Lam_k (r x r) per shift value.
%   Timing-only (state indicators), no element model. Given (U, V) the Lam step is a
%   convex linear solve over all shifts jointly.
%   Bar (same as week3_fullband_s2, fixed before running): operator NMSE <= -30 dB at
%   28 dB, and SER(est) <= 1.5 x SER(true) at 20 dB.

clear; rng(2300);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3;
M = 3276; bins = -M/2:M/2-1; evalb = round(M/4)+1:round(3*M/4);
Q = 16; K = 8; D = 0.75; Np = 2; Nd = 6; nit_cg = 25; r = 6;
base = [ones(1, round(D*K)), 2*ones(1, K - round(D*K))];
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
L = ceil(M*df*0.2e-6) + 10; [F, lamF] = dpss_basis(bins, df, 0.2e-6, L);
Tst = {fast_coupling(bins, fc, df, base, ones(1,M), ones(1,M), 'state1'), ...
       fast_coupling(bins, fc, df, base, ones(1,M), ones(1,M), 'state2')};
cfg = {20, 16; 28, 64}; ntrial = 5;
R = struct('M', M, 'L', L, 'rank', r, 'ntrial', ntrial);
for ic = 1:size(cfg, 1)
    snr_db = cfg{ic,1}; Mq = cfg{ic,2};
    ser = zeros(ntrial, 3); nm = zeros(ntrial, 1);
    for t = 1:ntrial
        shifts = randi(K, Q, 1) - 1;
        seqs = zeros(Q, K); for q = 1:Q, seqs(q,:) = circshift(base, shifts(q)); end
        [As, Ad] = phys_channel(Q, bins, df, fc, popt);
        op = fast_coupling(bins, fc, df, seqs, As, Ad, P);
        xt = randn(M, 6) + 1j*randn(M, 6); Hx = op.H(xt);
        s2 = mean(abs(Hx(:)).^2) / 2 / 10^(snr_db/10);
        S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
        Yp = op.H(S) + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
        d = randi(Mq, M, Nd) - 1; X = qammod_(d, Mq);
        Yd = op.H(X) + sqrt(s2/2)*(randn(M,Nd) + 1j*randn(M,Nd));
        ks = unique(shifts).';
        [U, V, Lam] = est_shift(Yp, S, Tst, F, lamF, s2, r, ks, K, bins);
        ope = shift_op(Tst, F, U, V, Lam, ks, K, bins);
        e = ope.H(xt) - Hx; nm(t) = sum(abs(e(:)).^2) / sum(abs(Hx(:)).^2);
        z = {Yd ./ onetap_diag(op, M), cg_mmse(op, Yd, s2, nit_cg), cg_mmse(ope, Yd, s2, nit_cg)};
        for j = 1:3, ser(t, j) = mean(qamdemod_(z{j}(evalb,:), Mq) ~= d(evalb,:), 'all'); end
        fprintf('  snr %d trial %d: op NMSE %.1f dB | SER one-tap %.3g true %.3g est %.3g\n', ...
            snr_db, t, 10*log10(nm(t)), ser(t,:));
    end
    key = sprintf('snr%d_qam%d', snr_db, Mq);
    R.(key) = struct('op_nmse_db_mean', 10*log10(mean(nm)), 'op_nmse_db_trials', 10*log10(nm).', ...
        'ser_mean', mean(ser,1), 'ser_ci95', 1.96*std(ser,0,1)/sqrt(ntrial), 'names', {{'onetap','true','est'}});
    fprintf('%s: op NMSE %.1f dB | SER one-tap %.3g | true %.3g | est %.3g (ratio %.2f)\n', key, ...
        10*log10(mean(nm)), mean(ser(:,1)), mean(ser(:,2)), mean(ser(:,3)), mean(ser(:,3))/mean(ser(:,2)));
end
fid = fopen(fullfile(out, 'week4_fullband_s2_shift.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

%% ---------------------------------------------------------------------
function D = dmod(bins, k, K)
D = exp(-1j*2*pi*bins(:)*k/K);                 % diag of D_k (M x 1)
end

function [U, V, Lam] = est_shift(Y, S, Tst, F, lam, s2, r, ks, K, bins)
[M, Np] = size(S); L = size(F, 2); y = Y(:); nk = numel(ks);
% G{s,p,k} = T_s (D_k^* S_p .* F)   (M x L): input-side regressors per shift.
% The output-side modulation D_k lives in FD{k} = D_k F below (applied once).
G = cell(2, Np, nk);
for ik = 1:nk
    Dk = dmod(bins, ks(ik), K);
    for s = 1:2, for p = 1:Np
        G{s,p,ik} = Tst{s}.H(conj(Dk) .* S(:,p) .* F);
    end, end
end
FD = cell(1, nk); for ik = 1:nk, FD{ik} = dmod(bins, ks(ik), K) .* F; end   % D_k F (output side)
% init: per-shift linear model with small common C_s (rank-limited by SVD), shared U
pl = kron(lam(:), lam(:)); Phi = zeros(M*Np, 2*L*L);
for s = 1:2, for p = 1:Np
    acc = zeros(M, L*L);
    for ik = 1:nk
        for i = 1:L, acc(:, (i-1)*L+(1:L)) = acc(:, (i-1)*L+(1:L)) + FD{ik}(:, i) .* G{s,p,ik}; end
    end
    Phi((p-1)*M+(1:M), (s-1)*L*L+(1:L*L)) = acc;
end, end
w = [pl; pl];
sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
Dw = sqrt(sc*w);
c = Dw .* (((Dw .* (Phi'*Phi)) .* Dw.' + s2*eye(numel(w))) \ (Dw .* (Phi'*y)));
[Uc, Sg, Wc] = svd([reshape(c(1:L*L), L, L).' reshape(c(L*L+1:end), L, L).']);
U = Uc(:, 1:r) * sqrt(Sg(1:r,1:r));
V = {conj(Wc(1:L, 1:r)) * sqrt(Sg(1:r,1:r)), conj(Wc(L+1:end, 1:r)) * sqrt(Sg(1:r,1:r))};
Lam = repmat({eye(r) / nk}, 1, nk);
wl = repmat(lam(:), r, 1); rr = 1e-3;
for it = 1:25
    % Lam-step (convex, all shifts jointly): column (ik, a, b) = (FD_k U(:,a)) .* (G_{s,k} V_s(:,b)) summed over s
    A = zeros(M*Np, nk*r*r);
    for ik = 1:nk
        FU = FD{ik} * U;
        for p = 1:Np
            GV = G{1,p,ik}*V{1} + G{2,p,ik}*V{2};          % M x r
            for a = 1:r
                A((p-1)*M+(1:M), (ik-1)*r*r+(a-1)*r+(1:r)) = FU(:, a) .* GV;
            end
        end
    end
    lv = rsolve(A, y, rr);
    for ik = 1:nk, Lam{ik} = reshape(lv((ik-1)*r*r+(1:r*r)), r, r).'; end
    % U-step: column (i, a) = sum_k FD_k(:,i) .* (G_k V Lam_k^T)(:, a)
    A = zeros(M*Np, L*r);
    for p = 1:Np
        for ik = 1:nk
            GV = (G{1,p,ik}*V{1} + G{2,p,ik}*V{2}) * Lam{ik}.';   % M x r
            for a = 1:r
                A((p-1)*M+(1:M), (a-1)*L+(1:L)) = A((p-1)*M+(1:M), (a-1)*L+(1:L)) + FD{ik} .* GV(:, a);
            end
        end
    end
    U = reshape(rsolve(A, y, rr), L, r);
    % V-step (both states): column (s, j, b) = sum_k (FD_k U Lam_k)(:, b) .* G_{s,k}(:, j)
    A = zeros(M*Np, 2*L*r);
    for p = 1:Np
        for ik = 1:nk
            FUL = FD{ik} * U * Lam{ik};                       % M x r
            for s = 1:2
                for b = 1:r
                    idx = (s-1)*L*r + (b-1)*L + (1:L);
                    A((p-1)*M+(1:M), idx) = A((p-1)*M+(1:M), idx) + FUL(:, b) .* G{s,p,ik};
                end
            end
        end
    end
    v = rsolve(A, y, rr);
    V = {reshape(v(1:L*r), L, r), reshape(v(L*r+1:end), L, r)};
end
end

function op = shift_op(Tst, F, U, V, Lam, ks, K, bins)
op.H = @(x) apply(x, false); op.Ht = @(y) apply(y, true);
    function out = apply(x, adj)
        out = zeros(size(x));
        for ik = 1:numel(ks)
            Dk = dmod(bins, ks(ik), K); FU = Dk .* (F*U*Lam{ik});   % M x r  output factors
            for s = 1:2
                FV = F*V{s};                                        % M x r  input factors (before D_k^*)
                for a = 1:size(FU, 2)
                    if ~adj
                        out = out + FU(:, a) .* Tst{s}.H(conj(Dk) .* FV(:, a) .* x);
                    else
                        out = out + Dk .* conj(FV(:, a)) .* Tst{s}.Ht(conj(FU(:, a)) .* x);
                    end
                end
            end
        end
    end
end

function x = rsolve(A, y, rel)
Gm = A'*A;
x = (Gm + rel*real(trace(Gm))/size(Gm,1)*eye(size(Gm))) \ (A'*y);
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
