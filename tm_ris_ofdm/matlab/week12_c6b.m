%WEEK12_C6B Learned-prior LMMSE under distribution shift, and a prior self-trained on our own estimates
%   (pre-registered 2026-10-25, Step 10). Ntrain = 300 per training set, 30 test channels, SNR {0,10,20,30}.
%   Shared family: train phys30 (true), phys100K3 (true), phys30 (self: ours, Np = 8, 20 dB); tests phys30, phys100K3.
%   TDL family (100 ns): train C/A/D (true), C (self); tests C, A, D. Ours on every test channel. Saves per channel.

if ~exist('ntrial', 'var'), ntrial = 30; end
if ~exist('ntrain', 'var'), ntrain = 300; end
if ~exist('outname', 'var'), outname = 'week12_c6b'; end
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
gen.phys30    = @() physH(struct('P', 5, 'K_db', 10, 'rms', 30e-9), Q, bins, df, fc, seqs, P, M);
gen.phys100K3 = @() physH(struct('P', 5, 'K_db', 3, 'rms', 100e-9), Q, bins, df, fc, seqs, P, M);
gen.TDLC = @() tdlH('C', Q, bins, df, fc, seqs, P, M);
gen.TDLA = @() tdlH('A', Q, bins, df, fc, seqs, P, M);
gen.TDLD = @() tdlH('D', Q, bins, df, fc, seqs, P, M);
fam(1) = struct('name', 'shared', 'tw', 0.3e-6, 'L', 16, 'r', 8, 'Np', 2, 'train', {{'phys30', 'phys100K3'}}, ...
    'self', 'phys30', 'tests', {{'phys30', 'phys100K3'}}, 'priors', {{{'phys30', 'self'}, {'phys100K3', 'phys30'}}});
fam(2) = struct('name', 'tdl', 'tw', 1e-6, 'L', ceil(M*df*1e-6) + 10, 'r', 16, 'Np', 4, 'train', {{'TDLC', 'TDLA', 'TDLD'}}, ...
    'self', 'TDLC', 'tests', {{'TDLC', 'TDLA', 'TDLD'}}, 'priors', {{{'TDLC', 'self'}, {'TDLA', 'TDLC'}, {'TDLD', 'TDLC'}}});
snrs = [0 10 20 30];
fp = fullfile(out, [outname '.json']);
if exist(fp, 'file'), R = jsondecode(fileread(fp)); else, R = struct(); end

for f = 1:2
    fm = fam(f); [F, lamF] = dpss_basis(bins, df, fm.tw, fm.L);
    % ---- training priors: true-H sets and the self-trained set
    pri = struct();
    for it = 1:numel(fm.train)
        X = zeros(M*M, ntrain); g = gen.(fm.train{it});
        for n = 1:ntrain, rng(14000 + 1000*f + 100*it + 10000*n); H = g(); X(:, n) = H(:); end
        [pri.(fm.train{it}).B, pri.(fm.train{it}).lam] = pca_prior(X);
        fprintf('%s: prior %s (true H) %d components\n', fm.name, fm.train{it}, numel(pri.(fm.train{it}).lam));
    end
    X = zeros(M*M, ntrain); g = gen.(fm.self);
    for n = 1:ntrain
        rng(15000 + 1000*f + 10000*n); H = g(); Pc = mean(sum(abs(H(:, evalb)).^2, 1)); s2 = Pc / 10^(20/10);
        S = qpsk(M, 8); Y = H*S + sqrt(s2/2)*(randn(M, 8) + 1j*randn(M, 8));
        Hh = est_general(Y, S, T, ones(Q,1), F, lamF, s2, fm.r, 20); X(:, n) = Hh(:);
    end
    [pri.self.B, pri.self.lam] = pca_prior(X); clear X
    fprintf('%s: self-trained prior %d components\n', fm.name, numel(pri.self.lam));
    % ---- tests
    for ie = 1:numel(fm.tests)
        tn = fm.tests{ie}; pl = fm.priors{ie}; ne = numel(pl) + 1; g = gen.(tn);
        rt = nan(ntrial, ne, numel(snrs)); nm = rt; rp = nan(ntrial, numel(snrs));
        for t = 1:ntrial
            rng(16000 + 1000*f + 100*ie + t); H = g(); Pc = mean(sum(abs(H(:, evalb)).^2, 1));
            S = qpsk(M, fm.Np); W0 = randn(M, fm.Np) + 1j*randn(M, fm.Np);
            for is = 1:numel(snrs)
                s2 = Pc / 10^(snrs(is)/10); Y = H*S + sqrt(s2/2)*W0;
                rp(t, is) = rate(H, H, s2, 0, evalb);
                for k = 1:ne
                    if k < ne, Hh = est_learned(Y, S, pri.(pl{k}).B, pri.(pl{k}).lam, s2, M);
                    else, Hh = est_general(Y, S, T, ones(Q,1), F, lamF, s2, fm.r, 20); end
                    e2 = norm(Y - Hh*S, 'fro')^2;
                    nm(t, k, is) = norm(Hh(:, evalb) - H(:, evalb), 'fro')^2 / norm(H(:, evalb), 'fro')^2;
                    rt(t, k, is) = rate(Hh, H, s2, max(e2/(M*fm.Np) - s2, 0), evalb);
                end
            end
            fprintf('%s test %s ch %2d | 20 dB gap (%s, ours) %s\n', fm.name, tn, t, strjoin(pl, ', '), ...
                sprintf('%.3f ', rp(t, 3) - rt(t, :, 3)));
            R.(tn) = struct('family', fm.name, 'snr_db', snrs, 'Np', fm.Np, 'ntrain', ntrain, ...
                'names', {[cellfun(@(p) ['learned_' p], pl, 'uni', 0), {'ours'}]}, ...
                'ncomp', cellfun(@(p) numel(pri.(p).lam), pl), ...
                'rate', rt(1:t,:,:), 'nmse', nm(1:t,:,:), 'rate_perfect', rp(1:t,:));
            save_json(fp, R);
        end
    end
    clear pri
end
fprintf('done\n');

%% ---------------------------------------------------------------------
function H = physH(opts, Q, bins, df, fc, seqs, P, M)
[As, Ad] = phys_channel(Q, bins, df, fc, opts); H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
end
function H = tdlH(m, Q, bins, df, fc, seqs, P, M)
[As, Ad] = tdl38901(Q, bins, df, fc, m, 100e-9); H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
end
function [B, lam] = pca_prior(X)
[V, D] = eig(X'*X); [d, o] = sort(real(diag(D)), 'descend'); V = V(:, o);
keep = d >= 1e-6*d(1); d = d(keep); V = V(:, keep);
B = (X*V) ./ sqrt(d.'); lam = d / size(X, 2);
end
function Hh = est_learned(Y, S, B, lam, s2, M)
Np = size(S, 2); Kc = size(B, 2);
Phi = reshape(pagemtimes(reshape(B, M, M, Kc), S), M*Np, Kc);
y = Y(:); Dl = sqrt(lam(:)); A = Phi .* Dl.';
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
