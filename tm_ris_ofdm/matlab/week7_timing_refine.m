%WEEK7_TIMING_REFINE Precise blind timing recovery for high SNR (point 3, follow-up).
%   week6_timing_offgrid: blind offset error rms 0.0033 slot -> 3.69 dB loss at 30 dB.
%   Here, after the coarse scan + golden section (8-iteration estimator), two refinements:
%     stage2  a second golden section on +/- 0.01 slot with the FULL estimator (20 it.)
%     stage2dd the same, but the residual also uses 6 decided data symbols (DD), so the
%              timing is fitted to 8 symbols instead of 2
%   Bar (fixed before running): at 30 dB, blind within 0.5 dB of known timing.
%   S1, D = 0.75, K = 8, Q = 16, M = 256, Np = 2, off-grid offsets, 12 trials.

clear; rng(3900);
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
loss = zeros(ntrial, 2, 4); err = zeros(ntrial, 2, 3);          % [coarse stage2 stage2dd known]
for t = 1:ntrial
    d0 = rand - 0.5;
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H = coupling_edges(bins, fc, df, repmat(((0:K) + d0)/K, Q, 1), st, As, Ad, P, M-1);
    Pc = mean(sum(abs(H(:, evalb)).^2, 1));
    for is = 1:2
        s2 = Pc / 10^(snrs(is)/10); Mq = qam(is);
        S = qpsk(M, Np); Y = H*S + noise(M, Np, s2);
        Xd = qammod_(randi(Mq, M, Nd) - 1, Mq); Yd = H*Xd + noise(M, Nd, s2);
        f8 = @(d) resid(d, Y, S, Tnom, bins, K, F, lamF, s2, Q, 8);
        r = arrayfun(f8, grid); [~, i0] = min(r);
        d1 = golden(f8, grid(max(i0-1,1)), grid(min(i0+1,numel(grid))), 14);
        f20 = @(d) resid(d, Y, S, Tnom, bins, K, F, lamF, s2, Q, 20);
        d2 = golden(f20, d1 - 0.01, d1 + 0.01, 16);
        % DD: decide data with the stage-2 estimate, then refit timing on pilots + data
        H2 = est_general(Y, S, shiftk(Tnom, bins, d2, K), ones(Q,1), F, lamF, s2, 8, 20);
        Xh = qamdec_(((H2'*H2 + s2*eye(M)) \ H2') * Yd, Mq);
        fdd = @(d) resid(d, [Y Yd], [S Xh], Tnom, bins, K, F, lamF, s2, Q, 20);
        d3 = golden(fdd, d2 - 0.005, d2 + 0.005, 16);
        rp = mmse_rate(H, H, s2, evalb);
        ds = [d1 d2 d3 d0]; lr = zeros(1, 4);
        for k = 1:4
            if k == 3, Hk = est_general([Y Yd], [S Xh], shiftk(Tnom, bins, ds(k), K), ones(Q,1), F, lamF, s2, 8, 20);
            else, Hk = est_general(Y, S, shiftk(Tnom, bins, ds(k), K), ones(Q,1), F, lamF, s2, 8, 20); end
            lr(k) = (rp - mmse_rate(Hk, H, s2, evalb)) / slope(is);
        end
        loss(t, is, :) = lr; err(t, is, :) = ds(1:3) - d0;
    end
    fprintf('trial %d done\n', t);
end
R = struct('snr_db', snrs, 'stages', {{'coarse', 'stage2', 'stage2dd', 'known'}}, ...
    'loss_db', squeeze(mean(loss, 1)), 'offset_err_rms', squeeze(sqrt(mean(err.^2, 1))));
for is = 1:2
    fprintf('%d dB: offset error rms coarse/stage2/stage2dd %.4f/%.4f/%.4f slot | loss coarse/stage2/stage2dd/known %.2f/%.2f/%.2f/%.2f dB\n', ...
        snrs(is), R.offset_err_rms(is, :), R.loss_db(is, :));
end
fid = fopen(fullfile(out, 'week7_timing_refine.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function r = resid(d, Y, S, Tnom, bins, K, F, lamF, s2, Q, it)
Hc = est_general(Y, S, shiftk(Tnom, bins, d, K), ones(Q,1), F, lamF, s2, 8, it);
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
function x = qamdec_(z, Mq)
m = sqrt(Mq); sc = sqrt(2*(Mq-1)/3); lv = -(m-1):2:(m-1);
q = @(v) lv(min(max(round((v + m - 1)/2), 0), m-1) + 1);
x = reshape((q(real(z)*sc) + 1j*q(imag(z)*sc)) / sc, size(z));
end
function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
