%WEEK11_THEORY Review G1 + G2 (pre-registered 2026-10-22), identical patterns (G = 1), M = 256.
%   Linear first stage y = Phi c + w, c ~ CN(0, P) (vec of the two L x L delay matrices, DPSS prior, scale
%   by moment matching), exact model for G = 1.
%   G1  effective degrees of freedom dof = tr(hat) = sum_i e_i/(e_i + s2), e_i eigenvalues of P^1/2 Phi^H Phi P^1/2,
%       versus observations M Np, for Np in {1, 2, 4, 8}; plus unknown counts of the low-rank model for G in {1,4,16}.
%   G2  Bayesian MSE bound (posterior covariance) of the coupling matrix,
%       BMSE_H = tr(Psi Sigma_post Psi^H), Psi = d vec(H) / d c, computed through the Toeplitz structure
%       (Psi^H Psi from conj(T_s) o T_s' and Khatri-Rao products of F), normalised by ||H||^2, versus the
%       NMSE achieved by the same LMMSE on the true channel. A gap means the DPSS prior is not the true
%       channel distribution (the bound holds only under the assumed prior).
%   Shared-path channel, 30 channels, 20/30 dB.
clear; rng(11600);
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; K = 8;
M = 256; bins = -M/2:M/2-1; Q = 16; evalb = 65:192;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lam] = dpss_basis(bins, df, 0.3e-6, 16); L = size(F, 2);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
% Psi^H Psi through the structure: <psi_(s,i,j), psi_(s',i',j')> = sum K_ss'(mb,m) conj(F(mb,i))F(mb,i') conj(F(m,j))F(m,j')
KR = zeros(M, L*L); for i = 1:L, for ip = 1:L, KR(:, (i-1)*L + ip) = conj(F(:, i)) .* F(:, ip); end, end
Gram = zeros(2*L*L);
for s = 1:2
    for sp = 1:2
        Kss = conj(T{s}) .* T{sp}; Z = Kss * KR;          % M x L^2 over (j, j')
        Gb = KR.' * Z;                                    % (i,i') x (j,j')
        blk = zeros(L*L);
        for i = 1:L, for j = 1:L, for ip = 1:L, for jp = 1:L
            blk((i-1)*L + j, (ip-1)*L + jp) = Gb((i-1)*L + ip, (j-1)*L + jp);
        end, end, end, end
        Gram((s-1)*L*L + (1:L*L), (sp-1)*L*L + (1:L*L)) = blk;
    end
end
Gram = (Gram + Gram')/2;
snrs = [20 30]; Nps = [1 2 4 8]; ntrial = 30;
dof = nan(ntrial, numel(Nps), 2); bnd = dof; ach = dof;
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1); Pc = mean(sum(abs(H(:, evalb)).^2, 1)); nH = norm(H, 'fro')^2;
    for is = 1:2
        s2 = Pc / 10^(snrs(is)/10); Sall = qpsk(M, 8); Yall = H*Sall + noise(M, 8, s2);
        for ip = 1:numel(Nps)
            S = Sall(:, 1:Nps(ip)); y = reshape(Yall(:, 1:Nps(ip)), [], 1); Np = Nps(ip);
            Phi = zeros(M*Np, 2*L*L);
            for s = 1:2, for p = 1:Np
                Gs = T{s} * (S(:,p) .* F);
                for i = 1:L, Phi((p-1)*M+(1:M), (s-1)*L*L+(i-1)*L+(1:L)) = F(:, i) .* Gs; end
            end, end
            w = repmat(kron(lam(:), lam(:)), 2, 1);
            sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
            D = sqrt(sc*w); A = (D .* (Phi'*Phi)) .* D.';
            e = real(eig((A + A')/2)); dof(t, ip, is) = sum(e ./ (e + s2));
            Sig = D .* ((A + s2*eye(numel(w))) \ (s2*diag(ones(numel(w), 1)))) .* D.';     % posterior covariance
            bnd(t, ip, is) = real(trace(Sig * Gram)) / nH;
            c = D .* ((A + s2*eye(numel(w))) \ (D .* (Phi'*y)));
            C1 = reshape(c(1:L*L), L, L).'; C2 = reshape(c(L*L+1:end), L, L).';
            Hh = T{1} .* (F*C1*F.') + T{2} .* (F*C2*F.');
            ach(t, ip, is) = norm(Hh - H, 'fro')^2 / nH;
        end
    end
    fprintf('trial %2d: 30 dB dof %s | bound/achieved NMSE dB Np=2 %.1f / %.1f\n', t, sprintf('%.0f ', squeeze(dof(t, :, 2))), ...
        10*log10(bnd(t, 2, 2)), 10*log10(ach(t, 2, 2)));
end
R = struct('ntrial', ntrial, 'snr_db', snrs, 'Np', Nps, 'L', L, 'unknowns_linear', 2*L*L, 'observations', M*Nps, ...
    'dof_mean', squeeze(mean(dof, 1)), 'bmse_bound_db', 10*log10(squeeze(mean(bnd, 1))), 'nmse_achieved_db', 10*log10(squeeze(mean(ach, 1))), ...
    'lowrank_unknowns_G1_4_16', [3*L*8 + 64, 3*L*8 + 4*64, 3*L*8 + 16*64]);
for is = 1:2
    fprintf('%d dB: dof %s of %s obs | bound %s dB | achieved %s dB\n', snrs(is), sprintf('%.0f ', R.dof_mean(:, is)), ...
        sprintf('%d ', M*Nps), sprintf('%.1f ', R.bmse_bound_db(:, is)), sprintf('%.1f ', R.nmse_achieved_db(:, is)));
end
fid = fopen(fullfile(out, 'week11_theory.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
function W = noise(M, n, s2), W = sqrt(s2/2) * (randn(M,n) + 1j*randn(M,n)); end
