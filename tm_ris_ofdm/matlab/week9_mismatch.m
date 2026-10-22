%WEEK9_MISMATCH Review step 2, workstream A3 (pre-registered 2026-10-22).
%   The ground truth is built slot by slot, H = sum_k T_k o (Aout_k Ain_k^T), where T_k is the Toeplitz
%   kernel of the indicator of slot k (sum_k T_k = I), Aout_k(:,q) the output response of element q in
%   slot k and Ain_k(:,q) its input response times its effective reflection in slot k. This lets each
%   mismatch be placed where the physics puts it, outside the estimator's model:
%     m1  state-dependent output response   a_d,q -> a_d,q (1 + eps eta_{q,s}),  eps = 0.05/0.1/0.2
%     m2  mutual coupling depending on BOTH elements' states (4-neighbour grid):
%         reflection of q in slot k += kappa sum_n eta_{qn}(s_q[k], s_n[k]) Gamma_{s_n[k]} a_s,n,
%         kappa = 0.02/0.05/0.08, identical patterns (G = 1) and 4 random patterns (G = 4)
%     m3  frequency-dependent array phase (beam squint), k = 2 pi (fc + f)/c, Q = 16 and Q = 256
%     m4  per-element resonance spread, sigma = 10 MHz per element and state
%   Estimator unchanged (est_general, tau_max 0.3 us, L 16, r 8, 20 it, Np 2, pilots only).
%   Matched reference: same channel without the mismatch, same pilots and noise. 30 channels, 20/30 dB.

clear; rng(9700);
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; c0 = 3e8; K = 8;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; f = fc + bins(:)*df; Np = 2;
base = [ones(1, 6), 2*ones(1, 2)];
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
Tk = cell(1, K);                                   % slot-indicator kernels
for k = 1:K
    stv = 2*ones(1, K); stv(k) = 1;
    Tk{k} = coupling_edges(bins, fc, df, (0:K)/K, stv, ones(1,M), ones(1,M), 'state1', M-1);
end
G0 = [element_response(f, 1, P), element_response(f, 2, P)];
kern = @(pat) {coupling_matrix(bins, fc, df, pat, ones(1,M), ones(1,M), 'state1', M-1), ...
               coupling_matrix(bins, fc, df, pat, ones(1,M), ones(1,M), 'state2', M-1)};
T1 = kern(base);
rng(9701); pats4 = zeros(4, K);
for g = 1:4, while true, p = 1 + (rand(1, K) > 0.5); if any(p == 1) && any(p == 2), break; end, end, pats4(g, :) = p; end
T4 = cell(4, 2); for g = 1:4, kk = kern(pats4(g, :)); T4{g,1} = kk{1}; T4{g,2} = kk{2}; end
rng(9700);
conds = {'m1_eps0.05', 'm1_eps0.10', 'm1_eps0.20', 'm2_G1_k0.02', 'm2_G1_k0.05', 'm2_G1_k0.08', ...
         'm2_G4_k0.02', 'm2_G4_k0.05', 'm2_G4_k0.08', 'm3_Q16', 'm3_Q256', 'm4_spread10MHz'};
snrs = [20 30]; ntrial = 30;
gm = nan(ntrial, numel(conds), 2); g0 = nan(ntrial, numel(conds), 2); relmis = nan(ntrial, numel(conds)); valerr = nan(ntrial, 1);
for t = 1:ntrial
    for ic = 1:numel(conds)
        cn = conds{ic};
        if contains(cn, 'Q256'), Q = 256; else, Q = 16; end
        if contains(cn, 'G4'), G = 4; grp = mod(0:Q-1, G).' + 1; seqs = pats4(grp, :); T = T4;
        else, G = 1; grp = ones(Q, 1); seqs = repmat(base, Q, 1); T = T1; end
        [pth, r] = paths(Q, fc);
        As = build(pth.s, r, bins, df, fc, c0, false); Ad = build(pth.d, r, bins, df, fc, c0, false);
        Gam = repmat(reshape(G0, M, 1, 2), 1, Q, 1);          % M x Q x 2 per-element state responses
        H0 = truth(Tk, seqs, As, Ad, Gam, [], [], 0);
        if ic == 1 && t <= 3
            Href = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
            valerr(t) = norm(H0 - Href, 'fro')/norm(Href, 'fro');
        end
        switch cn(1:2)
            case 'm1'
                eps_ = sscanf(cn, 'm1_eps%f'); eta = (randn(Q, 2) + 1j*randn(Q, 2))/sqrt(2);
                Hm = truth(Tk, seqs, As, Ad, Gam, 1 + eps_*eta, [], 0);
            case 'm2'
                kap = sscanf(cn(regexp(cn, 'k0'):end), 'k%f');
                C = coupling(Q, kap);
                Hm = truth(Tk, seqs, As, Ad, Gam, [], C, 1);
            case 'm3'
                As2 = build(pth.s, r, bins, df, fc, c0, true); Ad2 = build(pth.d, r, bins, df, fc, c0, true);
                Hm = truth(Tk, seqs, As2, Ad2, Gam, [], [], 0);
            case 'm4'
                Gs = Gam;
                for q = 1:Q
                    for s = 1:2
                        Pq = P; Pq.f0(s) = P.f0(s) + 10e6*randn; Gs(:, q, s) = element_response(f, s, Pq);
                    end
                end
                Hm = truth(Tk, seqs, As, Ad, Gs, [], [], 0);
        end
        relmis(t, ic) = norm(Hm - H0, 'fro')/norm(H0, 'fro');
        Pc = mean(sum(abs(H0(:, evalb)).^2, 1));
        for is = 1:2
            s2 = Pc / 10^(snrs(is)/10); S = qpsk(M, Np); W = noise(M, Np, s2);
            g0(t, ic, is) = mmse_rate(H0, H0, s2, evalb) - mmse_rate(est_general(H0*S + W, S, T, grp, F, lamF, s2, 8, 20), H0, s2, evalb);
            gm(t, ic, is) = mmse_rate(Hm, Hm, s2, evalb) - mmse_rate(est_general(Hm*S + W, S, T, grp, F, lamF, s2, 8, 20), Hm, s2, evalb);
        end
    end
    fprintf('trial %2d: 30 dB loss vs matched %s\n', t, sprintf('%+.2f ', squeeze(gm(t, :, 2) - g0(t, :, 2))));
end
R = struct('ntrial', ntrial, 'snr_db', snrs, 'conditions', {conds}, 'construction_check', max(valerr), ...
    'rel_mismatch_mean', mean(relmis, 1), 'gap_matched', squeeze(mean(g0, 1)), 'gap_mismatch', squeeze(mean(gm, 1)));
fprintf('construction check (slot-wise truth vs coupling_matrix): %.1e\n', max(valerr));
for ic = 1:numel(conds)
    for is = 1:2
        [m, h] = tci(gm(:, ic, is) - g0(:, ic, is));
        R.([strrep(strrep(conds{ic}, '.', 'p'), '-', '_') sprintf('_%ddB', snrs(is))]) = [m h];
    end
    fprintf('%-15s rel. mismatch %.3f | loss vs matched @20 %+.3f +-%.3f | @30 %+.3f +-%.3f\n', conds{ic}, mean(relmis(:, ic)), ...
        R.([strrep(strrep(conds{ic}, '.', 'p'), '-', '_') '_20dB']), R.([strrep(strrep(conds{ic}, '.', 'p'), '-', '_') '_30dB']));
end
R.gap_matched_trials = g0; R.gap_mismatch_trials = gm;
fid = fopen(fullfile(out, 'week9_mismatch.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

%% ---------------------------------------------------------------------
function H = truth(Tk, seqs, As, Ad, Gam, outfac, C, cmode)
% H = sum_k T_k o (Aout_k Ain_k^T). Gam: M x Q x 2. outfac: Q x 2 output factors per state (or []).
% C: Q x Q x 2 x 2 joint-state coupling (C(q,n,s_q,s_n)), used when cmode = 1.
[M, Q] = size(As.'); K = size(seqs, 2); H = zeros(M);
for k = 1:K
    sk = seqs(:, k);
    In = zeros(M, Q); Out = Ad.';
    for q = 1:Q, In(:, q) = Gam(:, q, sk(q)) .* As(q, :).'; end
    if cmode == 1
        In2 = In;
        for q = 1:Q
            for n = find(C(q, :, 1, 1) ~= 0 | C(q, :, 2, 2) ~= 0)
                In2(:, q) = In2(:, q) + C(q, n, sk(q), sk(n)) * In(:, n);
            end
        end
        In = In2;
    end
    if ~isempty(outfac)
        for q = 1:Q, Out(:, q) = Out(:, q) * outfac(q, sk(q)); end
    end
    H = H + Tk{k} .* (Out * In.');
end
end
function C = coupling(Q, kap)
% joint-state coupling between 4-neighbours on the sqrt(Q) x sqrt(Q) grid
n = round(sqrt(Q)); [gx, gz] = meshgrid(1:n, 1:n); gx = gx(:); gz = gz(:);
C = zeros(Q, Q, 2, 2);
for q = 1:Q
    nb = find(abs(gx - gx(q)) + abs(gz - gz(q)) == 1);
    for m = nb.', C(q, m, :, :) = kap * (randn(2, 2) + 1j*randn(2, 2))/sqrt(2); end
end
end
function [pth, r] = paths(Q, fc)
lam = 3e8/fc; n = round(sqrt(Q)); [gx, gz] = meshgrid(0:n-1, 0:n-1);
r = [gx(:) - (n-1)/2, zeros(Q,1), gz(:) - (n-1)/2] * lam/2;
pth = struct('s', hop(), 'd', hop());
    function h = hop()
        P_ = 5; Kr = 10; rms = 30e-9;
        tau = [0, sort(-rms*log(rand(1, P_)))];
        pw = [Kr/(Kr+1), (1/(Kr+1)) * exp(-tau(2:end)/rms) / sum(exp(-tau(2:end)/rms))];
        g = sqrt(pw) .* exp(1j*2*pi*rand(1, P_+1));
        u = randn(P_+1, 3); u(:,2) = abs(u(:,2)); u = u ./ vecnorm(u, 2, 2);
        h = struct('tau', tau, 'g', g, 'u', u);
    end
end
function A = build(h, r, bins, df, fc, c0, squint)
% Q x M per-element response; squint: array phase uses the actual frequency fc + f
Q = size(r, 1); M = numel(bins); A = zeros(Q, M); fr = bins(:).'*df;
for l = 1:numel(h.tau)
    proj = r*h.u(l, :).';                                  % Q x 1 path-length differences
    if squint, ph = exp(-1j*2*pi*(fc + fr).*proj/c0); else, ph = exp(-1j*2*pi*fc*proj/c0) * ones(1, M); end
    A = A + h.g(l) * ph .* exp(-1j*2*pi*h.tau(l)*fr);
end
end
function [m, h] = tci(x)
n = numel(x); nu = n - 1; tq = sqrt(nu*(1/betaincinv(0.05, nu/2, 0.5) - 1));
m = mean(x); h = tq*std(x)/sqrt(n);
end
function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
function W = noise(M, n, s2), W = sqrt(s2/2) * (randn(M,n) + 1j*randn(M,n)); end
function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
