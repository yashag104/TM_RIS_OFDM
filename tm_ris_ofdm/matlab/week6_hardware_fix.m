%WEEK6_HARDWARE_FIX Follow-up to week6_hardware: transients and a corrected phase-noise model.
%   (b2) Switching transients (rise time rho, fraction of a slot). The receiver builds
%        its timing kernels from a ramped waveform with rise time rho_hat:
%          nominal  rho_hat = 0 (instant switching, as before)
%          known    rho_hat = rho (rise time calibrated once for the hardware)
%          blind    rho_hat chosen on a grid by the smallest pilot residual (no truth)
%   (d2) Phase noise, CORRECTED model: one continuous Wiener phase across pilot 1,
%        pilot 2 and the data symbol (with the cyclic prefix in between), instead of a
%        fresh random phase per symbol (the bug in week6_hardware (d)).
%          raw  no common-phase correction
%          cpe  common phase error of each symbol removed (as NR phase-tracking RS do);
%               the remaining intra-symbol ICI is the true penalty
%   S1, D = 0.75, K = 8, Q = 16, M = 256, Np = 2, physical channel, 10 trials.

clear; rng(3400);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3; Tcp = 2.34e-6;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8; Np = 2; Nt = 32768;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
slope = [0.31 0.33]; snrs = [20 30]; ntrial = 10;
R = struct();

%% (b2) transients with calibrated or blindly estimated rise time
rhos = [0.05 0.1 0.2 0.3]; grid = 0:0.025:0.4;
Tg = cell(1, numel(grid));
for ig = 1:numel(grid), Tg{ig} = ramp_kernels(bins, fc, df, base, grid(ig), M, Nt); end
loss = zeros(ntrial, numel(rhos), 2, 3); rho_est = zeros(ntrial, numel(rhos));
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    for i = 1:numel(rhos)
        W = ramp_weights(base, Q, Nt, rhos(i));
        H = coupling_waveform(bins, fc, df, W, As, Ad, P, M-1);
        Pc = mean(sum(abs(H(:, evalb)).^2, 1));
        for is = 1:2
            s2 = Pc / 10^(snrs(is)/10);
            S = qpsk(M, Np); Y = H*S + noise(M, Np, s2);
            rp = mmse_rate(H, H, s2, evalb);
            Hn = est_general(Y, S, Tg{1}, ones(Q,1), F, lamF, s2, 8, 20);
            Tk = ramp_kernels(bins, fc, df, base, rhos(i), M, Nt);
            Hk = est_general(Y, S, Tk, ones(Q,1), F, lamF, s2, 8, 20);
            best = inf;
            for ig = 1:numel(grid)
                Hc = est_general(Y, S, Tg{ig}, ones(Q,1), F, lamF, s2, 8, 8);
                rs = norm(Y - Hc*S, 'fro'); if rs < best, best = rs; Hb = Hc; gb = grid(ig); end
            end
            loss(t, i, is, :) = (rp - [mmse_rate(Hn, H, s2, evalb), mmse_rate(Hk, H, s2, evalb), ...
                                       mmse_rate(Hb, H, s2, evalb)]) / slope(is);
            if is == 1, rho_est(t, i) = gb; end
        end
    end
end
R.b2_rise = rhos; R.b2_loss_db = squeeze(mean(loss, 1));          % rho x snr x [nominal known blind]
R.b2_rho_est_mean = mean(rho_est, 1); R.b2_rho_est_std = std(rho_est, 0, 1);
for i = 1:numel(rhos)
    fprintf('(b2) rise %.2f: est %.3f+/-%.3f | loss @20 nominal/known/blind %.2f/%.2f/%.2f | @30 %.2f/%.2f/%.2f\n', ...
        rhos(i), R.b2_rho_est_mean(i), R.b2_rho_est_std(i), squeeze(R.b2_loss_db(i,1,:)), squeeze(R.b2_loss_db(i,2,:)));
end

%% (d2) continuous phase noise, with and without common-phase correction
T0 = ramp_kernels(bins, fc, df, base, 0, M, Nt);
lws = [100 1e3 1e4 1e5]; loss = zeros(ntrial, numel(lws), 2, 2);
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H0 = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
    for i = 1:numel(lws)
        C = pn_three_symbols(M, lws(i), df, Tcp);                  % 3 symbols, continuous phase
        for mode = 1:2
            Hs = cell(1, 3);
            for n = 1:3
                Cn = C{n};
                if mode == 2, Cn = Cn * exp(-1j*angle(Cn(1,1)) ); end % remove common phase error
                Hs{n} = Cn * H0;
            end
            Pc = mean(sum(abs(Hs{3}(:, evalb)).^2, 1));
            for is = 1:2
                s2 = Pc / 10^(snrs(is)/10);
                S = qpsk(M, Np);
                Y = [Hs{1}*S(:,1), Hs{2}*S(:,2)] + noise(M, Np, s2);
                Hh = est_general(Y, S, T0, ones(Q,1), F, lamF, s2, 8, 20);
                loss(t, i, is, mode) = (mmse_rate(Hs{3}, Hs{3}, s2, evalb) - mmse_rate(Hh, Hs{3}, s2, evalb)) / slope(is);
            end
        end
    end
end
R.d2_linewidth_hz = lws; R.d2_loss_db = squeeze(mean(loss, 1));     % lw x snr x [raw cpe]
for i = 1:numel(lws)
    fprintf('(d2) linewidth %6.0f Hz: loss @20 raw/cpe %.2f/%.2f | @30 %.2f/%.2f\n', lws(i), ...
        squeeze(R.d2_loss_db(i,1,:)), squeeze(R.d2_loss_db(i,2,:)));
end
fid = fopen(fullfile(out, 'week6_hardware_fix.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

%% ---------------------------------------------------------------------
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

function C = pn_three_symbols(M, lw, df, Tcp)
% Continuous Wiener phase over three OFDM symbols separated by cyclic prefixes.
Ts = 1/(df*M); ncp = round(Tcp/Ts); ntot = 3*(M + ncp);
phi = cumsum(sqrt(2*pi*lw*Ts) * randn(ntot, 1));
C = cell(1, 3);
for n = 1:3
    s0 = (n-1)*(M + ncp) + ncp;                    % start of the useful part
    c = fft(exp(1j*phi(s0 + (1:M)))) / M;
    k = mod((0:M-1).' - (0:M-1), M);
    C{n} = c(k + 1);
end
end

function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
function W = noise(M, n, s2), W = sqrt(s2/2) * (randn(M,n) + 1j*randn(M,n)); end
function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
