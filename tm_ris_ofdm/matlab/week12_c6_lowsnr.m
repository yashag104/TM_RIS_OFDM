%WEEK12_C6_LOWSNR C6 learned-prior LMMSE baseline and low-SNR rank selection (pre-registered 2026-10-25, Step 9).
%   Part 1 (low-SNR rank): ours at ranks {1,2,4,8} (shared path, Np = 2) and {2,4,8,16} (TDL-C, Np = 4), SNR 0..30,
%     on 20 selection channels and 30 held-out channels; pilot residual stored for the BIC rule. Selection in Python.
%   Part 2 (C6): PCA prior of vec(H) learned from Ntrain in {100,300,1000} channels of the same generator; LMMSE
%     from the pilots. Ours on the same 30 test channels. SNR sweep and (shared path) Np sweep at 20 dB.
%   Set parts (default 1:2), ntrial, nsel, ntrain before calling to run a subset. Saves after every channel.

if ~exist('parts', 'var'), parts = 1:2; end
if ~exist('ntrial', 'var'), ntrial = 30; end
if ~exist('nsel', 'var'), nsel = 20; end
if ~exist('ntrains', 'var'), ntrains = [100 300 1000]; end
if ~exist('outname', 'var'), outname = 'week12_c6_lowsnr'; end
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; K = 8;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
cfg.phys30 = struct('tw', 0.3e-6, 'L', 16, 'r', 8, 'ranks', [1 2 4 8], 'Np', 2);
cfg.TDLC100 = struct('tw', 1e-6, 'L', ceil(M*df*1e-6) + 10, 'r', 16, 'ranks', [2 4 8 16], 'Np', 4);
chans = {'phys30', 'TDLC100'}; snrs = 0:5:30;
fp = fullfile(out, [outname '.json']);
if exist(fp, 'file'), R = jsondecode(fileread(fp)); else, R = struct(); end
mkH = @(ic) chan(ic, Q, bins, df, fc, popt, seqs, P, M);

%% ---------------- Part 1: low-SNR rank ----------------
if ismember(1, parts)
    for ic = 1:2
        cf = cfg.(chans{ic}); Np = cf.Np; nr = numel(cf.ranks);
        [F, lamF] = dpss_basis(bins, df, cf.tw, cf.L);
        for iset = 1:2
            if iset == 1, nch = nsel; seed0 = 13100; nm = 'selection'; else, nch = ntrial; seed0 = 13200; nm = 'heldout'; end
            rt = nan(nch, nr, numel(snrs)); res = rt; rp = nan(nch, numel(snrs));
            for t = 1:nch
                rng(seed0 + 100*ic + t);
                H = mkH(ic); Pc = mean(sum(abs(H(:, evalb)).^2, 1));
                S = qpsk(M, Np); W0 = randn(M, Np) + 1j*randn(M, Np);
                for is = 1:numel(snrs)
                    s2 = Pc / 10^(snrs(is)/10); Y = H*S + sqrt(s2/2)*W0;
                    rp(t, is) = rate(H, H, s2, 0, evalb);
                    for ir = 1:nr
                        Hh = est_general(Y, S, T, ones(Q,1), F, lamF, s2, cf.ranks(ir), 20);
                        e2 = norm(Y - Hh*S, 'fro')^2; res(t, ir, is) = e2 / s2;
                        rt(t, ir, is) = rate(Hh, H, s2, max(e2/(M*Np) - s2, 0), evalb);
                    end
                end
                fprintf('P1 %s %s ch %2d | 0 dB gap by rank %s | 30 dB %s\n', chans{ic}, nm, t, ...
                    sprintf('%.3f ', rp(t,1) - rt(t,:,1)), sprintf('%.3f ', rp(t,end) - rt(t,:,end)));
                R.P1.(chans{ic}).(nm) = struct('rate', rt(1:t,:,:), 'resid_over_s2', res(1:t,:,:), 'rate_perfect', rp(1:t,:));
                R.P1.(chans{ic}).ranks = cf.ranks; R.P1.(chans{ic}).Np = Np; R.P1.(chans{ic}).L = cf.L;
                R.P1.snr_db = snrs; R.P1.M = M;
                save_json(fp, R);
            end
        end
    end
end

%% ---------------- Part 2: C6 learned-prior LMMSE ----------------
if ismember(2, parts)
    Nps = [1 2 4 8 16 32];
    for ic = 1:2
        cf = cfg.(chans{ic}); Np0 = cf.Np;
        [F, lamF] = dpss_basis(bins, df, cf.tw, cf.L);
        % training set (single precision to bound memory: M^2 x Ntrain complex)
        Nt = max(ntrains); X = zeros(M*M, Nt, 'single'); tic;
        for n = 1:Nt
            rng(13000 + 10000*ic + n); Hn = mkH(ic); X(:, n) = single(Hn(:));
        end
        fprintf('P2 %s: %d training channels in %.0f s\n', chans{ic}, Nt, toc);
        B = cell(1, numel(ntrains)); lam = B;
        for k = 1:numel(ntrains)
            Xk = double(X(:, 1:ntrains(k)));
            [V, D] = eig(Xk'*Xk); [d, o] = sort(real(diag(D)), 'descend'); V = V(:, o);
            keep = d >= 1e-6*d(1); d = d(keep); V = V(:, keep);
            B{k} = (Xk*V) ./ sqrt(d.');                       % orthonormal PCA basis, M^2 x Kk
            lam{k} = d / ntrains(k);                           % prior variances of the coefficients
            fprintf('P2 %s: Ntrain %d -> %d components\n', chans{ic}, ntrains(k), numel(d));
        end
        clear X Xk V
        ne = numel(ntrains) + 1;                               % learned (each Ntrain) + ours
        rt = nan(ntrial, ne, numel(snrs)); nmse = rt; rp = nan(ntrial, numel(snrs));
        rtn = nan(ntrial, ne, numel(Nps)); nmn = rtn; rpn = nan(ntrial, 1);
        for t = 1:ntrial
            rng(13300 + 100*ic + t);
            H = mkH(ic); Pc = mean(sum(abs(H(:, evalb)).^2, 1));
            Sp = qpsk(M, 32); W0 = randn(M, 32) + 1j*randn(M, 32);
            for is = 1:numel(snrs)
                s2 = Pc / 10^(snrs(is)/10); Yall = H*Sp + sqrt(s2/2)*W0;
                S = Sp(:, 1:Np0); Y = Yall(:, 1:Np0);
                rp(t, is) = rate(H, H, s2, 0, evalb);
                for k = 1:ne
                    if k < ne, Hh = est_learned(Y, S, B{k}, lam{k}, s2, M);
                    else, Hh = est_general(Y, S, T, ones(Q,1), F, lamF, s2, cf.r, 20); end
                    e2 = norm(Y - Hh*S, 'fro')^2;
                    nmse(t, k, is) = norm(Hh(:, evalb) - H(:, evalb), 'fro')^2 / norm(H(:, evalb), 'fro')^2;
                    rt(t, k, is) = rate(Hh, H, s2, max(e2/(M*Np0) - s2, 0), evalb);
                end
                if ic == 1 && snrs(is) == 20
                    rpn(t) = rp(t, is);
                    for ip = 1:numel(Nps)
                        Np = Nps(ip); S = Sp(:, 1:Np); Y = Yall(:, 1:Np);
                        for k = 1:ne
                            if k < ne, Hh = est_learned(Y, S, B{k}, lam{k}, s2, M);
                            else, Hh = est_general(Y, S, T, ones(Q,1), F, lamF, s2, cf.r, 20); end
                            e2 = norm(Y - Hh*S, 'fro')^2;
                            nmn(t, k, ip) = norm(Hh(:, evalb) - H(:, evalb), 'fro')^2 / norm(H(:, evalb), 'fro')^2;
                            rtn(t, k, ip) = rate(Hh, H, s2, max(e2/(M*Np) - s2, 0), evalb);
                        end
                    end
                end
            end
            fprintf('P2 %s ch %2d | 20 dB gap (learned %s; ours) %s\n', chans{ic}, t, mat2str(ntrains), ...
                sprintf('%.3f ', rp(t, 5) - rt(t, :, 5)));
            C = struct('snr_db', snrs, 'Np', Np0, 'ntrain', ntrains, 'ncomp', cellfun(@numel, lam), ...
                'names', {[arrayfun(@(n) sprintf('learned%d', n), ntrains, 'uni', 0), {'ours'}]}, ...
                'rate', rt(1:t,:,:), 'nmse', nmse(1:t,:,:), 'rate_perfect', rp(1:t,:));
            if ic == 1
                C.Np_sweep = Nps; C.rate_np = rtn(1:t,:,:); C.nmse_np = nmn(1:t,:,:); C.rate_perfect_np = rpn(1:t);
            end
            R.P2.(chans{ic}) = C; save_json(fp, R);
        end
        clear B
    end
end
fprintf('done\n');

%% ---------------------------------------------------------------------
function H = chan(ic, Q, bins, df, fc, popt, seqs, P, M)
if ic == 1, [As, Ad] = phys_channel(Q, bins, df, fc, popt);
else, [As, Ad] = tdl38901(Q, bins, df, fc, 'C', 100e-9); end
H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
end
function Hh = est_learned(Y, S, B, lam, s2, M)
% LMMSE with the learned PCA prior: vec(H) = B c, c ~ CN(0, diag(lam)); y = vec(H S) + n
Np = size(S, 2); Kc = size(B, 2);
Phi = reshape(pagemtimes(reshape(B, M, M, Kc), S), M*Np, Kc);
y = Y(:); Dl = sqrt(lam(:));
A = Phi .* Dl.';                                            % Phi diag(sqrt(lam))
c = Dl .* (A' * ((A*A' + s2*eye(M*Np)) \ y));
Hh = reshape(B*c, M, M);
end
function r = rate(Hh, H, s2, se2, evalb)
M = size(H, 1);
W = (Hh'*Hh + (s2 + se2)*eye(M)) \ Hh'; G = W*H; g = diag(G);
it = sum(abs(G).^2, 2) - abs(g).^2 + s2*sum(abs(W).^2, 2);
r = mean(log2(1 + abs(g(evalb)).^2 ./ it(evalb)));
end
function save_json(fp, R)
fid = fopen(fp, 'w'); fprintf(fid, '%s', jsonencode(R)); fclose(fid);
end
function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
