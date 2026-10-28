%WEEK13_REVIEW_BASELINES Review round 2 (pre-registered 2026-10-28, Step 13a), M = 256.
%   S1: flat-calibrated ideal element (per-state complex scalars Gamma_s(fc) instead of +-1), and the same with the
%       ratio Gamma_2/Gamma_1 estimated from the pilots (from our element-free linear stage).
%   S2: structure-aware baselines and ablations: linear stage alone; ours without the DPSS basis (delay-grid basis,
%       flat prior); exact-element genie WITH the low-rank structure; per-state static pilots fed into Prop. 1
%       (H ~ sum_s T_s diag(h_s)); OMP over a delay x delay x state dictionary that uses the known kernels T.
%       'ours' is est_general as in all earlier runs (with the pilot-fit fallback); 'ours_nofb' is the ALS output
%       (Algorithm 1 as written in the paper). Variants of ours use the no-fallback output.
%   S3: sensitivity of ours to tau_max (x0.5, x2) and rank r (shared 4/16, TDL 8/24); effective degrees of freedom of
%       the linear stage (trace of its LMMSE hat matrix) per channel and SNR.
%   Shared path Np = 2, TDL-C 100 ns Np = 4; SNR 10/20/30 dB; 30 new channels (rng 13300+); metrics NMSE and
%   robust-MMSE rate loss as week12. Set ntrial, icsel, outname to run a subset. Saves after every channel.

if ~exist('ntrial', 'var'), ntrial = 30; end
if ~exist('icsel', 'var'), icsel = 1:2; end
if ~exist('outname', 'var'), outname = 'week13_review_baselines'; end
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; K = 8;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16; fin = fc + bins*df;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
cfg.phys30 = struct('tw', 0.3e-6, 'L', 16, 'r', 8, 'rlo', 4, 'rhi', 16, 'tw2', 0.6e-6);
cfg.TDLC100 = struct('tw', 1e-6, 'L', ceil(M*df*1e-6) + 10, 'r', 16, 'rlo', 8, 'rhi', 24, 'tw2', 2e-6);
G1 = element_response(fin, 1, P); G2 = element_response(fin, 2, P);            % true responses (genie only)
g1c = element_response(fc, 1, P); g2c = element_response(fc, 2, P);            % flat calibration at the carrier
snrs = [10 20 30];
ests = {'ours', 'ours_nofb', 'linear', 'ours_dft', 'ours_tw05', 'ours_tw2', 'ours_rlo', 'ours_rhi', ...
        'genie', 'genie_lr', 'ideal', 'flatcal', 'flatest', 'static', 'omp_T', 'bemall_Np32'};
chans = {'phys30', 'TDLC100'}; npc = [2 4];
fp = fullfile(out, [outname '.json']);
if exist(fp, 'file'), R = jsondecode(fileread(fp)); else, R = struct(); end
for ic = icsel
    cf = cfg.(chans{ic}); Np = npc(ic);
    [F, lamF] = dpss_basis(bins, df, cf.tw, cf.L);
    Lh = ceil(M*df*cf.tw/2) + (cf.L - ceil(M*df*cf.tw)); Ld = ceil(M*df*cf.tw*2) + (cf.L - ceil(M*df*cf.tw));
    [Fh, lamh] = dpss_basis(bins, df, cf.tw/2, max(Lh, 4));                    % tau_max x0.5 (same margin rule)
    [Fd, lamd] = dpss_basis(bins, df, cf.tw*2, Ld);                            % tau_max x2
    Fg = exp(-1j*2*pi*bins(:)*df*((0:cf.L-1)/cf.L*cf.tw)); Fg = Fg ./ vecnorm(Fg); % delay grid, no DPSS
    L2 = ceil(M*df*cf.tw2) + 12; [F2, lam2] = dpss_basis(bins, df, cf.tw2, L2);
    taus = 0:1/(4*M*df):cf.tw;                                                 % OMP delay grid
    nm = nan(ntrial, numel(ests), numel(snrs)); rt = nm; rp = nan(ntrial, numel(snrs)); dof = rp; fbu = rp;
    for t = 1:ntrial
        rng(13300 + 100*ic + t);
        if ic == 1, [As, Ad] = phys_channel(Q, bins, df, fc, popt);
        else, [As, Ad] = tdl38901(Q, bins, df, fc, 'C', 100e-9); end
        H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
        hst = {sum(Ad .* As, 1) .* G1, sum(Ad .* As, 1) .* G2};               % static RIS, state s: diagonal channel
        Pc = mean(sum(abs(H(:, evalb)).^2, 1));
        Sp = qpsk(M, 32); W0 = randn(M, 32) + 1j*randn(M, 32); Ws = randn(M, Np) + 1j*randn(M, Np);
        for is = 1:numel(snrs)
            s2 = Pc / 10^(snrs(is)/10);
            Yall = H*Sp + sqrt(s2/2)*W0; S = Sp(:, 1:Np); Y = Yall(:, 1:Np);
            rp(t, is) = metrics(H, H, s2, 0, evalb);
            [Hfb, info] = est_general_ab(Y, S, T, ones(Q,1), F, lamF, s2, cf.r, 20, struct('fallback', true));
            fbu(t, is) = info.used_fallback;
            [C1, C2, dof(t, is)] = lin_stage(Y, S, T, F, lamF, s2);
            for ie = 1:numel(ests)
                switch ests{ie}
                    case 'ours',      Hh = Hfb;                    % identical to est_general (instrumented copy)
                    case 'ours_nofb', Hh = info.Hals;
                    case 'linear',    Hh = info.H0;
                    case 'ours_dft',  Hh = nofb(Y, S, T, Fg, ones(cf.L,1)/cf.L, s2, cf.r);
                    case 'ours_tw05', Hh = nofb(Y, S, T, Fh, lamh, s2, cf.r);
                    case 'ours_tw2',  Hh = nofb(Y, S, T, Fd, lamd, s2, cf.r);
                    case 'ours_rlo',  Hh = nofb(Y, S, T, F, lamF, s2, cf.rlo);
                    case 'ours_rhi',  Hh = nofb(Y, S, T, F, lamF, s2, cf.rhi);
                    case 'genie',     Hh = est_ideal(Y, S, {T{1} .* G1, -T{2} .* G2}, F, lamF, s2);
                    case 'genie_lr',  Kg = T{1} .* G1 + T{2} .* G2; Hh = nofb(Y, S, {Kg, zeros(M)}, F, lamF, s2, cf.r);
                    case 'ideal',     Hh = est_ideal(Y, S, T, F, lamF, s2);
                    case 'flatcal',   Hh = est_ideal(Y, S, {g1c * T{1}, -g2c * T{2}}, F, lamF, s2);
                    case 'flatest',   rho = (C1(:)' * C2(:)) / max(real(C1(:)' * C1(:)), eps);
                                      Hh = est_ideal(Y, S, {T{1}, -rho * T{2}}, F, lamF, s2);
                    case 'static',    Hh = est_static(S, hst, Ws, s2, T, F2, lam2);
                    case 'omp_T',     Hh = est_omp_T(Y, S, T, s2, taus, bins, df);
                    case 'bemall_Np32', Hh = est_bem_cg(Yall, Sp, F2, lam2, s2, -(M-1):(M-1));
                end
                se_r = max(norm(Y - Hh*S, 'fro')^2/(M*Np) - s2, 0);
                if strcmp(ests{ie}, 'bemall_Np32'), se_r = max(norm(Yall - Hh*Sp, 'fro')^2/(M*32) - s2, 0); end
                if strcmp(ests{ie}, 'static'), se_r = 0; end                    % its pilots are not Y; oracle-free 0
                nm(t, ie, is) = norm(Hh(:, evalb) - H(:, evalb), 'fro')^2 / norm(H(:, evalb), 'fro')^2;
                rt(t, ie, is) = metrics(Hh, H, s2, se_r, evalb);
            end
        end
        fprintf('%s ch %2d | 20 dB loss %s | dof %.1f\n', chans{ic}, t, sprintf('%.2f ', rp(t, 2) - rt(t, :, 2)), dof(t, 2));
        R.(chans{ic}) = struct('snr_db', snrs, 'Np', Np, 'estimators', {ests}, 'nmse', nm(1:t,:,:), ...
            'rate', rt(1:t,:,:), 'rate_perfect', rp(1:t,:), 'dof_linear', dof(1:t,:), 'fallback_used', fbu(1:t,:), ...
            'L', cf.L, 'r', cf.r, 'n_unknowns_lowrank', 3*cf.L*cf.r + cf.r^2, 'n_obs', M*Np);
        save_json(fp, R);
    end
end
fprintf('done\n');

%% ---------------------------------------------------------------------
function Hh = nofb(Y, S, T, F, lam, s2, r)
[~, info] = est_general_ab(Y, S, T, ones(16,1), F, lam, s2, r, 20, struct('fallback', false));
Hh = info.Hals;
end
function [C1, C2, dof] = lin_stage(Y, S, T, F, lam, s2)
% element-free linear stage of Algorithm 1 (H ~ sum_s T_s o (F C_s F^T)); returns C_s and the trace of the hat matrix
[M, Np] = size(S); L = size(F, 2); y = Y(:);
Phi = zeros(M*Np, 2*L*L);
for s = 1:2, for p = 1:Np
    Gs = T{s} * (S(:,p) .* F);
    for i = 1:L, Phi((p-1)*M+(1:M), (s-1)*L*L+(i-1)*L+(1:L)) = F(:, i) .* Gs; end
end, end
pl = kron(lam(:), lam(:)); w = [pl; pl];
sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
Dw = sqrt(sc*w); A = (Dw .* (Phi'*Phi)) .* Dw.'; A = (A + A')/2;
c = Dw .* ((A + s2*eye(numel(w))) \ (Dw .* (Phi'*y)));
C1 = reshape(c(1:L*L), L, L).'; C2 = reshape(c(L*L+1:end), L, L).';
ev = max(real(eig(A)), 0); dof = sum(ev ./ (ev + s2));
end
function Hh = est_static(S, hst, Ws, s2, T, F2, lam2)
% per-state static pilots: the RIS holds state s during Np/2 pilot symbols, giving y = h_s .* x + n; LMMSE of h_s in
% the delay basis, then Prop. 1 with the harmonic span assumed flat: H(mb, m) ~ sum_s T_s(mb, m) h_s(m)
[M, Np] = size(S); nps = max(1, floor(Np/2)); L = size(F2, 2); Hh = zeros(M);
for s = 1:2
    idx = (s-1)*nps + (1:nps); idx = mod(idx-1, Np) + 1;
    X = S(:, idx); y = reshape(hst{s}(:) .* X + sqrt(s2/2)*Ws(:, idx), [], 1);
    Phi = zeros(M*nps, L);
    for p = 1:nps, Phi((p-1)*M+(1:M), :) = X(:, p) .* F2; end
    sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * lam2(:))/numel(y), eps);
    Dw = sqrt(sc*lam2(:));
    c = Dw .* (((Dw .* (Phi'*Phi)) .* Dw.' + s2*eye(L)) \ (Dw .* (Phi'*y)));
    Hh = Hh + T{s} .* (F2*c).';
end
end
function Hh = est_omp_T(Y, S, T, s2, taus, bins, df)
% sparse estimator that uses the known kernels: atoms T_s o (e(tau_a) e(tau_b)^T), s = 1, 2
[M, Np] = size(S); y = Y(:); E = exp(-1j*2*pi*bins(:)*df*taus); nt = numel(taus);
D = zeros(M*Np, 2*nt*nt);
for s = 1:2, for p = 1:Np
    Gs = T{s} * (S(:,p) .* E);                                               % M x nt (input delay b)
    for a = 1:nt, D((p-1)*M+(1:M), (s-1)*nt*nt + (a-1)*nt + (1:nt)) = E(:, a) .* Gs; end
end, end
Dn = D ./ vecnorm(D); res = y; sel = []; target = s2 * numel(y);
for it = 1:min(300, floor(numel(y)/2))
    [~, j] = max(abs(Dn' * res)); sel(end+1) = j; %#ok<AGROW>
    c = D(:, sel) \ y; res = y - D(:, sel) * c;
    if real(res'*res) <= target, break; end
end
Hh = zeros(M);
for ii = 1:numel(sel)
    j = sel(ii) - 1; s = floor(j/(nt*nt)) + 1; j = j - (s-1)*nt*nt; a = floor(j/nt) + 1; b = j - (a-1)*nt + 1;
    Hh = Hh + c(ii) * (T{s} .* (E(:, a) * E(:, b).'));
end
end
function r = metrics(Hh, H, s2, se2, evalb)
M = size(H, 1);
W = (Hh'*Hh + (s2 + se2)*eye(M)) \ Hh'; G = W*H; g = diag(G);
it = sum(abs(G).^2, 2) - abs(g).^2 + s2*sum(abs(W).^2, 2);
r = mean(log2(1 + abs(g(evalb)).^2 ./ it(evalb)));
end
function Hh = est_ideal(Y, S, T, F, lam, s2)
% model-aware LMMSE: H = (T{1} - T{2}) o (F C F^T); T{s} may carry the element response as column scaling
[M, Np] = size(S); L = size(F, 2); y = Y(:); Td = T{1} - T{2};
Phi = zeros(M*Np, L*L);
for p = 1:Np
    Gs = Td * (S(:,p) .* F);
    for i = 1:L, Phi((p-1)*M+(1:M), (i-1)*L+(1:L)) = F(:, i) .* Gs; end
end
w = kron(lam(:), lam(:));
sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
Dw = sqrt(sc*w);
c = Dw .* (((Dw .* (Phi'*Phi)) .* Dw.' + s2*eye(numel(w))) \ (Dw .* (Phi'*y)));
Hh = Td .* (F*reshape(c, L, L).'*F.');
end
function save_json(fp, R)
fid = fopen(fp, 'w'); fprintf(fid, '%s', jsonencode(R)); fclose(fid);
end
function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
