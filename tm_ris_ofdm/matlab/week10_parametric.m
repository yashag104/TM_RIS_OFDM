%WEEK10_PARAMETRIC Review step 4: C3 parametric element-fit baseline (pre-registered 2026-10-22).
%   The receiver knows the single-pole family but starts from parameters off by N(0,10 MHz) (resonance) and
%   N(0,10%) (decay rates). It fits the 6 parameters by maximising the evidence of the model-aware linear
%   LMMSE, H = B(theta) o (F C F^T), B = sum_s T_s o (1 Gamma_s(theta)^T) (fminsearch, 300 iterations), then
%   estimates with the fitted model. Compared with timing-only (est_general), the oracle model and the
%   unfitted wrong model. Truth (t1): exact single-pole cell; (t2) per-element resonance spread 10 MHz.
%   Shared-path channel, Np = 2, 30 channels, 20/30 dB.
clear; rng(10300);
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; K = 8;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16; Np = 2; f = fc + bins(:)*df;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
Bof = @(PP) T{1} .* (ones(M,1) * element_response(f, 1, PP).') + T{2} .* (ones(M,1) * element_response(f, 2, PP).');
snrs = [20 30]; ntrial = 30; names = {'timing_only', 'oracle_model', 'wrong_model', 'fitted_model'};
gap = nan(ntrial, 4, 2, 2); perr = nan(ntrial, 2, 2);     % ... x snr x truth
opts = optimset('MaxIter', 300, 'MaxFunEvals', 600, 'Display', 'off', 'TolX', 1e-3, 'TolFun', 1e-3);
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    P0 = P; P0.f0 = P.f0 + 10e6*randn(1, 2); P0.xr = P.xr .* (1 + 0.1*randn(1, 2)); P0.xi = P.xi .* (1 + 0.1*randn(1, 2));
    Plist = cell(Q, 1); for q = 1:Q, Pq = P; Pq.f0 = P.f0 + 10e6*randn(1, 2); Plist{q} = Pq; end
    Hs = {coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1), coupling_spread(bins, fc, df, seqs, As, Ad, Plist, M)};
    for it = 1:2
        H = Hs{it}; Pc = mean(sum(abs(H(:, evalb)).^2, 1));
        for is = 1:2
            s2 = Pc / 10^(snrs(is)/10); S = qpsk(M, Np); Y = H*S + noise(M, Np, s2); rp = mmse_rate(H, H, s2, evalb);
            th2P = @(th) setp(P0, th);
            obj = @(th) -evidence_B(Bof(th2P(th)), Y, S, F, lamF, s2);
            th = fminsearch(obj, zeros(1, 6), opts); Pf = th2P(th);
            perr(t, is, it) = abs(Pf.f0(1) - P.f0(1))/1e6;
            Hh = {est_general(Y, S, T, ones(Q,1), F, lamF, s2, 8, 20), est_lmmse(Y, S, {Bof(P)}, F, lamF, s2, 'linear', 15), ...
                  est_lmmse(Y, S, {Bof(P0)}, F, lamF, s2, 'linear', 15), est_lmmse(Y, S, {Bof(Pf)}, F, lamF, s2, 'linear', 15)};
            for e = 1:4, gap(t, e, is, it) = rp - mmse_rate(Hh{e}, H, s2, evalb); end
        end
    end
    fprintf('trial %2d: 30 dB gaps (t1) %s | (t2) %s | f0 err after fit %.1f MHz\n', t, sprintf('%.3f ', gap(t, :, 2, 1)), sprintf('%.3f ', gap(t, :, 2, 2)), perr(t, 2, 1));
end
R = struct('ntrial', ntrial, 'snr_db', snrs, 'names', {names}, 'truths', {{'single_pole_exact', 'element_spread_10MHz'}}, ...
    'gap_mean', squeeze(mean(gap, 1)), 'dims', 'estimator x snr x truth', 'f0_err_after_fit_MHz', squeeze(mean(perr, 1)));
for it = 1:2
    for is = 1:2
        [m, h] = tci(gap(:, 4, is, it) - gap(:, 1, is, it));
        R.(sprintf('fitted_minus_timing_%s_%ddB', R.truths{it}, snrs(is))) = [m h];
        fprintf('%-22s %d dB: gaps %s | fitted - timing %+.3f +-%.3f\n', R.truths{it}, snrs(is), sprintf('%.3f ', mean(gap(:, :, is, it), 1)), m, h);
    end
end
R.gap_trials = gap;
fid = fopen(fullfile(out, 'week10_parametric.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function PP = setp(P0, th)
PP = P0; PP.f0 = P0.f0 + 1e6*th(1:2); PP.xr = P0.xr .* exp(th(3:4)); PP.xi = P0.xi .* exp(th(5:6));
end
function H = coupling_spread(bins, fc, df, seqs, As, Ad, Plist, M)
H = zeros(M);
for q = 1:numel(Plist), H = H + coupling_matrix(bins, fc, df, seqs(q, :), As(q, :), Ad(q, :), Plist{q}, M-1); end
end
function le = evidence_B(B, Y, S, F, lam, s2)
[M, Np] = size(S); L = size(F, 2); y = Y(:); Phi = zeros(M*Np, L*L);
for p = 1:Np
    Gs = B * (S(:,p) .* F);
    for i = 1:L, Phi((p-1)*M+(1:M), (i-1)*L+(1:L)) = F(:, i) .* Gs; end
end
w = kron(lam(:), lam(:));
sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
Dw = sqrt(sc*w); Bm = (Dw .* (Phi'*Phi)) .* Dw.' / s2 + eye(numel(w));
Rc = chol((Bm + Bm')/2); z = Rc' \ (Dw .* (Phi'*y)) / s2;
le = -(real(y'*y)/s2 - real(z'*z)) - 2*sum(log(real(diag(Rc))));
end

function [m, h] = tci(x)
x = x(~isnan(x)); n = numel(x); nu = n - 1; tq = sqrt(nu*(1/betaincinv(0.05, nu/2, 0.5) - 1));
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
