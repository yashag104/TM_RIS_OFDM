%WEEK11_RISE_SMALL Review K2 (pre-registered 2026-10-22): realistic nanosecond switching rise times.
%   rho in {0.0005, 0.002, 0.005, 0.01, 0.02} slot (2, 8, 21, 42, 83 ns). Truth with the ramp; our receiver
%   with nominal (instant) kernels and with calibrated kernels. Shared-path channel, 30 channels, 20/30 dB.
clear; rng(11700);
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; K = 8; Nt = 32768;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16; Np = 2;
base = [ones(1, 6), 2*ones(1, 2)];
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
Tn = ramp_kernels(bins, fc, df, base, 0, M, Nt);
rhos = [0.0005 0.002 0.005 0.01 0.02]; snrs = [20 30]; ntrial = 30;
gap = nan(ntrial, numel(rhos), 2, 2);                     % [nominal, calibrated]
for ir = 1:numel(rhos)
    Tc = ramp_kernels(bins, fc, df, base, rhos(ir), M, Nt); Wq = ramp_weights(base, Q, Nt, rhos(ir));
    for t = 1:ntrial
        [As, Ad] = phys_channel(Q, bins, df, fc, popt);
        H = coupling_waveform(bins, fc, df, Wq, As, Ad, P, M-1); Pc = mean(sum(abs(H(:, evalb)).^2, 1));
        for is = 1:2
            s2 = Pc / 10^(snrs(is)/10); S = qpsk(M, Np); Y = H*S + noise(M, Np, s2); rp = mmse_rate(H, H, s2, evalb);
            gap(t, ir, 1, is) = rp - mmse_rate(est_general(Y, S, Tn, ones(Q,1), F, lamF, s2, 8, 20), H, s2, evalb);
            gap(t, ir, 2, is) = rp - mmse_rate(est_general(Y, S, Tc, ones(Q,1), F, lamF, s2, 8, 20), H, s2, evalb);
        end
    end
    fprintf('rho %.4f: gap nominal/calibrated @20 %.3f/%.3f @30 %.3f/%.3f\n', rhos(ir), squeeze(mean(gap(:, ir, :, 1), 1)), squeeze(mean(gap(:, ir, :, 2), 1)));
end
R = struct('rho', rhos, 'rise_ns', rhos*1e9/(df*K), 'snr_db', snrs, 'ntrial', ntrial, 'gap_mean', squeeze(mean(gap, 1)), 'dims', 'rho x [nominal calibrated] x snr');
R.gap_trials = gap;
fid = fopen(fullfile(out, 'week11_rise_small.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);
function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
function W = noise(M, n, s2), W = sqrt(s2/2) * (randn(M,n) + 1j*randn(M,n)); end
function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
function W = ramp_weights(base, Q, Nt, rho)
K = numel(base); t = ((0:Nt-1) + 0.5) / Nt * K;
w1 = double(base(floor(t) + 1) == 1);
if rho > 0
    for k = 0:K-1
        a = base(mod(k-1, K) + 1); b = base(k + 1);
        if a == b, continue; end
        d = t - k; d(d > K/2) = d(d > K/2) - K; d(d < -K/2) = d(d < -K/2) + K;
        in = abs(d) < rho/2; ramp = (d(in) + rho/2) / rho;
        if b == 1, w1(in) = ramp; else, w1(in) = 1 - ramp; end
    end
end
W = zeros(Q, Nt, 2); W(:, :, 1) = repmat(w1, Q, 1); W(:, :, 2) = repmat(1 - w1, Q, 1);
end

function T = ramp_kernels(bins, fc, df, base, rho, M, Nt)
% Timing kernels from the (possibly ramped) switching waveform, state indicators.
W = ramp_weights(base, 1, Nt, rho);
T = {coupling_waveform(bins, fc, df, W, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_waveform(bins, fc, df, W, ones(1,M), ones(1,M), 'state2', M-1)};
end

