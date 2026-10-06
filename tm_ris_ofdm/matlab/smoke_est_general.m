% Smoke test for est_general: random patterns, compare with perfect H and (G=1) est_lmmse.
clear; rng(7);
P = params_openris(); fc = 3.594e9; df = 30e3; M = 256; bins = -M/2:M/2-1; evalb = 65:192;
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9); [F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
Q = 16; K = 8;
for G = [1 4 16]
  for Np = [2 4]
    seqG = randi(2, G, K); for g = 1:G, while all(seqG(g,:) == seqG(g,1)), seqG(g,:) = randi(2,1,K); end, end
    grp = mod((0:Q-1).', G) + 1; seqs = seqG(grp, :);
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
    T = cell(G,2); for g=1:G, for s=1:2, T{g,s} = coupling_matrix(bins, fc, df, seqG(g,:), ones(1,M), ones(1,M), sprintf('state%d',s), M-1); end, end
    s2 = mean(sum(abs(H(:,evalb)).^2,1)) / 100;
    S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np)))/sqrt(2); Y = H*S + sqrt(s2/2)*(randn(M,Np)+1j*randn(M,Np));
    tic; Hg = est_general(Y, S, T, grp, F, lamF, s2, 8, 20); tg = toc;
    e = @(A) 10*log10(norm(A(evalb,evalb)-H(evalb,evalb),'fro')^2/norm(H(evalb,evalb),'fro')^2);
    msg = sprintf('G=%2d Np=%d: est_general NMSE %.1f dB (%.1fs)', G, Np, e(Hg), tg);
    if G == 1, Hl = est_lmmse(Y, S, {T{1,1},T{1,2}}, F, lamF, s2, 'linear', 15); msg = [msg sprintf(' | est_lmmse %.1f dB', e(Hl))]; end
    disp(msg);
  end
end
