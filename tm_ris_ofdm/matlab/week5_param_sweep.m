%WEEK5_PARAM_SWEEP Are the results an artefact of our chosen settings?
%   Sweeps every design choice one at a time around the default
%   (Q = 16, K = 8, D = 0.75, identical sequences), physical channel, M = 256,
%   Np = 2 pilot symbols, SNR 20 and 30 dB, 12 channels per point.
%     A) number of elements Q in {4, 16, 64, 144}
%     B) slots per symbol K in {2, 4, 8, 16} (duty 0.5, and 0.75 where K allows)
%     C) duty cycle D in {1..7}/8 with K = 8
%     D) random switching patterns: G distinct random K-slot sequences assigned to
%        the 16 elements, G in {1, 2, 4, 8, 16}, Np in {2, 4}
%   Receivers: one-tap; MMSE with CE-BEM (Q=12, LMMSE prior); MMSE with ours
%   (timing-only est_general, same settings at every point); perfect H.
%   Reported: rate (bit/s/Hz) mean and 95% CI; ours' gap to perfect.

clear; rng(2600);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
L2 = ceil(M*df*0.6e-6) + 12; [F2, lam2] = dpss_basis(bins, df, 0.6e-6, L2);
snrs = [20 30]; ntrial = 12; R = struct();

% A) number of elements
for Q = [4 16 64 144]
    R.(sprintf('A_Q%d', Q)) = run_point(Q, 8, 6, 1, 2, ntrial, snrs, P, fc, df, M, bins, evalb, popt, F, lamF, F2, lam2);
end
% B) slots per symbol
for K = [2 4 8 16]
    for D = [0.5 0.75]
        n1 = round(D*K); if n1 < 1 || n1 >= K || abs(n1/K - D) > 1e-9, continue; end
        R.(sprintf('B_K%d_D%03d', K, round(1000*D))) = run_point(16, K, n1, 1, 2, ntrial, snrs, P, fc, df, M, bins, evalb, popt, F, lamF, F2, lam2);
    end
end
% C) duty cycle
for n1 = 1:7
    R.(sprintf('C_D%d_8', n1)) = run_point(16, 8, n1, 1, 2, ntrial, snrs, P, fc, df, M, bins, evalb, popt, F, lamF, F2, lam2);
end
% D) random patterns, G distinct sequences
for G = [1 2 4 8 16]
    for Np = [2 4]
        R.(sprintf('D_G%d_Np%d', G, Np)) = run_point(16, 8, -1, G, Np, ntrial, snrs, P, fc, df, M, bins, evalb, popt, F, lamF, F2, lam2);
    end
end
fid = fopen(fullfile(out, 'week5_param_sweep.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);
disp('done');

function res = run_point(Q, K, n1, G, Np, ntrial, snrs, P, fc, df, M, bins, evalb, popt, F, lamF, F2, lam2)
% n1 = slots in state 1 (identical duty pattern); n1 = -1 means random sequences.
names = {'onetap', 'bem', 'ours', 'perfect'};
rate = zeros(ntrial, numel(names), numel(snrs));
for t = 1:ntrial
    if n1 > 0
        base = [ones(1, n1), 2*ones(1, K - n1)];
        seqG = base; grp = ones(Q, 1);
    else
        seqG = zeros(G, K);
        for g = 1:G
            sq = randi(2, 1, K); while all(sq == sq(1)), sq = randi(2, 1, K); end   % not constant
            seqG(g, :) = sq;
        end
        grp = mod((0:Q-1).', G) + 1; grp = grp(randperm(Q));
    end
    seqs = seqG(grp, :);
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
    Gn = size(seqG, 1);
    T = cell(Gn, 2);
    for g = 1:Gn
        for s = 1:2
            T{g,s} = coupling_matrix(bins, fc, df, seqG(g,:), ones(1,M), ones(1,M), sprintf('state%d', s), M-1);
        end
    end
    Pc = mean(sum(abs(H(:, evalb)).^2, 1));
    for is = 1:numel(snrs)
        s2 = Pc / 10^(snrs(is)/10);
        S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
        Y = H*S + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
        g1 = mean(Y ./ S, 2); w = conj(g1) ./ abs(g1).^2; dH = diag(H);
        sinr = abs(w.*dH).^2 ./ (abs(w).^2 .* (sum(abs(H).^2, 2) - abs(dH).^2 + s2));
        rate(t, 1, is) = mean(log2(1 + sinr(evalb)));
        rate(t, 2, is) = mmse_rate(est_bem(Y, S, F2, lam2, s2, -12:12), H, s2, evalb);
        Hh = est_general(Y, S, T, grp, F, lamF, s2, 8, 20);   % one estimator for every point
        rate(t, 3, is) = mmse_rate(Hh, H, s2, evalb);
        rate(t, 4, is) = mmse_rate(H, H, s2, evalb);
    end
end
res = struct('Q', Q, 'K', K, 'n1', n1, 'G', G, 'Np', Np, 'snr_db', snrs, 'names', {names}, ...
    'rate_mean', squeeze(mean(rate, 1)), 'rate_ci95', squeeze(1.96*std(rate, 0, 1)/sqrt(ntrial)));
gap = squeeze(mean(rate(:,4,:) - rate(:,3,:), 1));
fprintf('Q=%3d K=%2d n1=%2d G=%2d Np=%d | ours %s | perfect %s | bem %s | onetap %s | gap %s\n', Q, K, n1, G, Np, ...
    sprintf('%.2f ', res.rate_mean(3,:)), sprintf('%.2f ', res.rate_mean(4,:)), ...
    sprintf('%.2f ', res.rate_mean(2,:)), sprintf('%.2f ', res.rate_mean(1,:)), sprintf('%.2f ', gap));
end

function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; Gm = W*H; sg = abs(diag(Gm)).^2;
sinr = sg ./ (sum(abs(Gm).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
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
