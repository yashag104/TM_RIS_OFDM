%WEEK8_TIMING_EVIDENCE (L2b) Blind timing by Bayesian evidence instead of residual.
%   A clock offset d is exactly a delay of the output hop and an advance of the input
%   hop (T_d o (FU)(FV)^T = T o (D FU)(D^* FV)^T). A DPSS basis still represents the
%   shifted factors (projection error ~1e-4), so the pilot RESIDUAL is nearly flat in d,
%   which is why week7's residual search failed at 30 dB. The shifted factors are,
%   however, improbable under the DPSS prior. The evidence of the linear first stage,
%       y = Phi_d c + w,  c ~ CN(0, sc diag(w)),  C_d = Phi_d diag(sc w) Phi_d^H + s2 I,
%       log p(y | d) = -y^H C_d^{-1} y - log det C_d,
%   accounts for this and should peak at the true d.
%   Search: grid +-0.5 slot step 0.025, then golden section (14 steps), as week7.
%   Receivers: resid (week7 coarse stage: residual of the 8-iteration estimator),
%   evid (this), known. Final estimate always est_general (20 it.) with the chosen d.
%   Bar (pre-registered): 30 dB, evid within 1 dB of known.
%   S1, K = 8, Q = 16, M = 256, Np = 2, off-grid offsets, physical channel, 12 trials.

clear; rng(3900);                                   % same seed/trials as week7_timing_refine
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8; Np = 2; Nd = 6;
base = [ones(1, 6), 2*ones(1, 2)]; st = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
E0 = (0:K)/K;
Tnom = {coupling_edges(bins, fc, df, E0, base, ones(1,M), ones(1,M), 'state1', M-1), ...
        coupling_edges(bins, fc, df, E0, base, ones(1,M), ones(1,M), 'state2', M-1)};
snrs = [20 30]; qam = [16 64]; slope = [0.31 0.33]; ntrial = 12; grid = -0.5:0.025:0.5;
loss = zeros(ntrial, 2, 3); err = zeros(ntrial, 2, 2);
for t = 1:ntrial
    d0 = rand - 0.5;
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H = coupling_edges(bins, fc, df, repmat(((0:K) + d0)/K, Q, 1), st, As, Ad, P, M-1);
    Pc = mean(sum(abs(H(:, evalb)).^2, 1));
    for is = 1:2
        s2 = Pc / 10^(snrs(is)/10); Mq = qam(is);
        S = qpsk(M, Np); Y = H*S + noise(M, Np, s2);
        Xd = qammod_(randi(Mq, M, Nd) - 1, Mq); Yd = H*Xd + noise(M, Nd, s2); %#ok<NASGU> keep RNG as week7
        fr = @(d) resid(d, Y, S, Tnom, bins, K, F, lamF, s2, Q);
        r = arrayfun(fr, grid); [~, i0] = min(r);
        dr = golden(fr, grid(max(i0-1,1)), grid(min(i0+1,numel(grid))), 14);
        fe = @(d) -evidence(d, Y, S, Tnom, bins, K, F, lamF, s2);
        e = arrayfun(fe, grid); [~, i0] = min(e);
        de = golden(fe, grid(max(i0-1,1)), grid(min(i0+1,numel(grid))), 14);
        rp = mmse_rate(H, H, s2, evalb); ds = [dr de d0];
        for k = 1:3
            Hk = est_general(Y, S, shiftk(Tnom, bins, ds(k), K), ones(Q,1), F, lamF, s2, 8, 20);
            loss(t, is, k) = (rp - mmse_rate(Hk, H, s2, evalb)) / slope(is);
        end
        err(t, is, :) = ds(1:2) - d0;
    end
    fprintf('trial %d: d0 %+.3f | err resid %+.4f evid %+.4f (20 dB) %+.4f %+.4f (30 dB)\n', t, d0, err(t,1,:), err(t,2,:));
end
R = struct('snr_db', snrs, 'receivers', {{'resid', 'evid', 'known'}}, 'ntrial', ntrial, ...
    'loss_db', squeeze(mean(loss, 1)), 'loss_ci95', squeeze(1.96*std(loss, 0, 1)/sqrt(ntrial)), ...
    'offset_err_rms', squeeze(sqrt(mean(err.^2, 1))), 'loss_trials', loss);
for is = 1:2
    fprintf('%d dB: offset err rms resid %.4f evid %.4f slot | loss resid/evid/known %.2f/%.2f/%.2f dB\n', ...
        snrs(is), R.offset_err_rms(is, :), R.loss_db(is, :));
end
fid = fopen(fullfile(out, 'week8_timing_evidence.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function le = evidence(d, Y, S, Tnom, bins, K, F, lam, s2)
T = shiftk(Tnom, bins, d, K); [M, Np] = size(S); L = size(F, 2); y = Y(:);
Phi = zeros(M*Np, 2*L*L);
for s = 1:2
    for p = 1:Np
        Gs = T{s} * (S(:,p) .* F);
        for i = 1:L, Phi((p-1)*M+(1:M), (s-1)*L*L+(i-1)*L+(1:L)) = F(:, i) .* Gs; end
    end
end
pl = kron(lam(:), lam(:)); w = [pl; pl];
sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
Dw = sqrt(sc*w); B = (Dw .* (Phi'*Phi)) .* Dw.' / s2 + eye(numel(w));
Rc = chol((B + B')/2);
z = Rc' \ (Dw .* (Phi'*y)) / s2;
le = -(real(y'*y)/s2 - real(z'*z)) - 2*sum(log(real(diag(Rc))));   % up to a d-independent constant
end
function r = resid(d, Y, S, Tnom, bins, K, F, lamF, s2, Q)
Hc = est_general(Y, S, shiftk(Tnom, bins, d, K), ones(Q,1), F, lamF, s2, 8, 8);
r = norm(Y - Hc*S, 'fro');
end
function T = shiftk(Tnom, bins, d, K)
D = exp(-1j*2*pi*bins(:)*d/K);
T = {D .* Tnom{1} .* D', D .* Tnom{2} .* D'};
end
function x = golden(f, a, b, n)
g = (sqrt(5) - 1)/2; c = b - g*(b - a); d = a + g*(b - a); fc = f(c); fd = f(d);
for i = 1:n
    if fc < fd, b = d; d = c; fd = fc; c = b - g*(b - a); fc = f(c);
    else, a = c; c = d; fc = fd; d = a + g*(b - a); fd = f(d); end
end
x = (a + b)/2;
end
function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
function W = noise(M, n, s2), W = sqrt(s2/2) * (randn(M,n) + 1j*randn(M,n)); end
function x = qammod_(d, Mq)
m = sqrt(Mq); lv = -(m-1):2:(m-1);
x = (lv(mod(d, m) + 1) + 1j*lv(floor(d / m) + 1)) / sqrt(2*(Mq-1)/3);
end
function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
