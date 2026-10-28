% Full band, one channel, very high SNR, 4 pilots: NMSE of the linear estimators. An exact model should approach the
% delay-basis floor; a mis-specified one should floor at its model error. Prints only.
P = params_openris(); fc = 3.594e9; df = 30e3; M = 3276; bins = -M/2:M/2-1; Q = 16; f = fc + bins(:)*df;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1); evalb = round(M/4)+1:round(3*M/4);
L = ceil(M*df*0.2e-6) + 10; [F, lamF] = dpss_basis(bins, df, 0.2e-6, L);
Tst = {fast_coupling(bins, fc, df, base, ones(1,M), ones(1,M), 'state1'), fast_coupling(bins, fc, df, base, ones(1,M), ones(1,M), 'state2')};
rng(13501); [As, Ad] = phys_channel(Q, bins, df, fc, struct('P', 5, 'K_db', 10, 'rms', 30e-9));
op = fast_coupling(bins, fc, df, seqs, As, Ad, P);
cols = evalb(round(linspace(1, numel(evalb), 48))); E = zeros(M, 48); E(sub2ind(size(E), cols, 1:48)) = 1; HE = op.H(E);
S = (sign(randn(M,4)) + 1j*sign(randn(M,4))) / sqrt(2);
gt = {element_response(f, 1, P), element_response(f, 2, P)}; gf = {element_response(fc, 1, P), element_response(fc, 2, P)};
for snr = [20 40 60]
    s2 = mean(sum(abs(HE).^2, 1)) / 10^(snr/10); Y = op.H(S) + sqrt(s2/2)*(randn(M,4) + 1j*randn(M,4));
    o = {lin_op(Y, S, Tst, gt, F, lamF, s2, bins, fc, df, base), lin_op(Y, S, Tst, gf, F, lamF, s2, bins, fc, df, base), lin_op(Y, S, Tst, {1, -1}, F, lamF, s2, bins, fc, df, base)};
    fprintf('%d dB: NMSE genie %.1f | flatcal %.1f | ideal %.1f dB\n', snr, cellfun(@(x) 10*log10(norm(x.H(E) - HE, 'fro')^2/norm(HE, 'fro')^2), o));
end
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
