%WEEK3_FULLBAND_S2_DIAG (diagnostic: floor vs SNR, and restarts) Per-element sequences at full band: shared path-subspace estimator.
%   All elements see the same propagation paths, so their hop responses share a
%   low-dimensional subspace:  a_q = U*alpha_q,  s_q = V*beta_q  (U, V: L x r shared)
%       H = sum_q B_q o ( (F U alpha_q) (F V beta_q).' )
%   Unknowns 2*L*r + 2*Q*r (vs 2*Q*L for independent rank-1 per element).
%   Block ALS over (U, alpha, V, beta), relative ridge per block.
%
%   Kill criterion (fixed before running): with Np = 2 at M = 3276,
%     operator NMSE <= -30 dB at 28 dB, and SER(est) <= 1.5 x SER(true H) at 20 dB.
%   Compared against the previous independent rank-1 estimator (from week3_fullband).

clear; rng(1000);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3;
M = 3276; bins = -M/2:M/2-1; evalb = round(M/4)+1:round(3*M/4);
Q = 16; K = 8; D = 0.75; Np = 2; Nd = 4; nit_cg = 25;
base = [ones(1, round(D*K)), 2*ones(1, K - round(D*K))];
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
tau = 0.2e-6; L = ceil(M*df*tau) + 10;
[F, lamF] = dpss_basis(bins, df, tau, L);
cfg = {28, 64; 60, 64}; ranks = [6]; ntrial = 2;
R = struct('M', M, 'L', L, 'ranks', ranks);
for ic = 1:size(cfg, 1)
    snr_db = cfg{ic,1}; Mq = cfg{ic,2};
    nm = zeros(1, numel(ranks)); err = zeros(1, 2 + numel(ranks)); nsym = 0;
    for t = 1:ntrial
        seqs = zeros(Q,K); for q = 1:Q, seqs(q,:) = circshift(base, randi(K)-1); end
        [As, Ad] = phys_channel(Q, bins, df, fc, popt);
        op = fast_coupling(bins, fc, df, seqs, As, Ad, P);
        Bop = cell(1,Q); for q = 1:Q, Bop{q} = fast_coupling(bins, fc, df, seqs(q,:), ones(1,M), ones(1,M), P); end
        xt = randn(M, 6) + 1j*randn(M, 6); Hx = op.H(xt);
        s2 = mean(abs(Hx(:)).^2) / 2 / 10^(snr_db/10) * 2 / 2;   % unit-power symbols
        S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
        Yp = op.H(S) + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
        d = randi(Mq, M, Nd) - 1; X = qammod_(d, Mq);
        Yd = op.H(X) + sqrt(s2/2)*(randn(M,Nd) + 1j*randn(M,Nd));
        G = cell(Q, Np);
        for q = 1:Q, for p = 1:Np, G{q,p} = Bop{q}.H(S(:,p) .* F); end, end
        z = {cg_mmse(op, Yd, s2, nit_cg)};
        for ir = 1:numel(ranks)
            best = inf;
            for rs = 0:3
                [Ae1, Se1] = est_subspace(Yp, G, F, lamF, s2, ranks(ir), Q, rs);
                op1 = fast_coupling(bins, fc, df, seqs, Se1, Ae1, P);
                res = norm(op1.H(S) - Yp, 'fro')^2;                   % pilot residual (no truth used)
                fprintf('    restart %d: pilot residual %.3e, op NMSE %.1f dB\n', rs, res, ...
                    10*log10(sum(abs(op1.H(xt) - Hx).^2, 'all') / sum(abs(Hx).^2, 'all')));
                if res < best, best = res; Ae = Ae1; Se = Se1; end
            end
            ope = fast_coupling(bins, fc, df, seqs, Se, Ae, P);
            e = ope.H(xt) - Hx;
            nm(ir) = nm(ir) + sum(abs(e(:)).^2) / sum(abs(Hx(:)).^2);
            z{end+1} = cg_mmse(ope, Yd, s2, nit_cg); %#ok<AGROW>
        end
        dg = onetap_diag(op, M);
        z = [{Yd ./ dg}, z];
        for e = 1:numel(z)
            err(e) = err(e) + sum(qamdemod_(z{e}(evalb,:), Mq) ~= d(evalb,:), 'all');
        end
        nsym = nsym + numel(evalb)*Nd;
    end
    key = sprintf('snr%d_qam%d', snr_db, Mq);
    R.(key) = struct('ser_onetap', err(1)/nsym, 'ser_true', err(2)/nsym, ...
        'ser_est', err(3:end)/nsym, 'op_nmse_db', 10*log10(nm/ntrial));
    fprintf('%s: SER one-tap %.3g | true %.3g | est r=%s: %s | op NMSE %s dB\n', key, err(1)/nsym, err(2)/nsym, ...
        mat2str(ranks), sprintf('%.3g ', err(3:end)/nsym), sprintf('%.1f ', 10*log10(nm/ntrial)));
end
fid = fopen(fullfile(out, 'week3_fullband_s2_diag.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function [Ae, Se] = est_subspace(Yp, G, F, lam, s2, r, Q, restart)
% Block ALS for H = sum_q B_q o ((F U a_q)(F V b_q).'), with G{q,p} = B_q (S_p .* F).
[M, L] = size(F); Np = size(G, 2); y = Yp(:);
% init: common C from mean regressor (as if all B_q were equal), top-r SVD
pl = kron(lam(:), lam(:)); Phi = zeros(M*Np, L*L);
for p = 1:Np
    Gm = zeros(M, L); for q = 1:Q, Gm = Gm + G{q,p}; end
    for i = 1:L, Phi((p-1)*M+(1:M), (i-1)*L+(1:L)) = F(:, i) .* Gm; end
end
sc0 = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * pl)/numel(y), eps);
Dw = sqrt(sc0*pl);
c = Dw .* (((Dw .* (Phi'*Phi)) .* Dw.' + s2*eye(L*L)) \ (Dw .* (Phi'*y)));
[Uc, Sg, Vc] = svd(reshape(c, L, L).');
U = Uc(:, 1:r) * sqrt(Sg(1:r,1:r)); V = conj(Vc(:, 1:r)) * sqrt(Sg(1:r,1:r));
al = ones(r, Q); be = ones(r, Q);
if restart > 0                                    % random init on the same scale
    sU = norm(U, 'fro'); sV = norm(V, 'fro');
    U = (randn(L, r) + 1j*randn(L, r)) .* sqrt(lam(:)); U = U / norm(U, 'fro') * sU;
    V = (randn(L, r) + 1j*randn(L, r)) .* sqrt(lam(:)); V = V / norm(V, 'fro') * sV;
    al = randn(r, Q) + 1j*randn(r, Q); be = randn(r, Q) + 1j*randn(r, Q);
end
rr = 1e-3;
for it = 1:30
    FV = F*V; FU = F*U;
    % g{q,p} = B_q (S_p .* (F V b_q)) = G{q,p} * (V b_q)
    % U-step: column (i,k) = F(:,i) .* sum_q al(k,q) g{q,p}
    A = zeros(M*Np, L*r);
    for p = 1:Np
        acc = zeros(M, r);
        for q = 1:Q, acc = acc + (G{q,p} * (V*be(:,q))) * al(:,q).'; end
        for k = 1:r, A((p-1)*M+(1:M), (k-1)*L+(1:L)) = F .* acc(:, k); end
    end
    U = reshape(rsolve(A, y, rr), L, r); FU = F*U;
    % alpha-step: column (q,k) = (F U(:,k)) .* g{q,p}
    A = zeros(M*Np, r*Q);
    for p = 1:Np
        for q = 1:Q
            g = G{q,p} * (V*be(:,q));
            A((p-1)*M+(1:M), (q-1)*r+(1:r)) = FU .* g;
        end
    end
    al = reshape(rsolve(A, y, rr), r, Q);
    % V-step: column (j,k) = sum_q be(k,q) (F U a_q) .* G{q,p}(:,j)
    A = zeros(M*Np, L*r);
    for p = 1:Np
        for k = 1:r
            acc = zeros(M, L);
            for q = 1:Q, acc = acc + be(k,q) * ((FU*al(:,q)) .* G{q,p}); end
            A((p-1)*M+(1:M), (k-1)*L+(1:L)) = acc;
        end
    end
    V = reshape(rsolve(A, y, rr), L, r); FV = F*V;
    % beta-step: column (q,k) = (F U a_q) .* (G{q,p} V(:,k))
    A = zeros(M*Np, r*Q);
    for p = 1:Np
        for q = 1:Q
            A((p-1)*M+(1:M), (q-1)*r+(1:r)) = (FU*al(:,q)) .* (G{q,p} * V);
        end
    end
    be = reshape(rsolve(A, y, rr), r, Q);
end
Ae = (F*U*al).'; Se = (F*V*be).';
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
