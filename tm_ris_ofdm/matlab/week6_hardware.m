%WEEK6_HARDWARE (#5) Hardware imperfections, receiver unchanged (timing-only, nominal).
%   (a) element-to-element spread: each element's resonances shifted by N(0, sf),
%       decay rates scaled by (1 + N(0, 0.05)); sf in {0, 5, 10, 20} MHz
%   (b) switching transients: linear ramp between states lasting rho of a slot,
%       rho in {0, 0.05, 0.1, 0.2, 0.3} (breaks the instant-settling assumption)
%   (c) receiver carrier-frequency offset eps (in subcarriers) in {0, .005, .01, .02, .05}
%   (d) oscillator phase noise (Wiener), linewidth in {0, 100, 1e3, 1e4} Hz, new
%       realisation every OFDM symbol (pilots and data differ)
%   S1, D = 0.75, K = 8, Q = 16, M = 256, Np = 2, physical channel, 10 trials.
%   Loss (dB, rate-equivalent) vs perfect knowledge of the impaired channel.

clear; rng(2900);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8; Np = 2; Nt = 32768;   % 32768 points: <0.07% vs exact edges
base = [ones(1, 6), 2*ones(1, 2)];
seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
Tnom = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
        coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
slope = [0.31 0.33]; snrs = [20 30]; ntrial = 10;
R = struct();

%% (a) element spread
sfs = [0 5e6 10e6 20e6]; loss = zeros(ntrial, numel(sfs), 2);
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    for i = 1:numel(sfs)
        Pl = cell(Q, 1);
        for q = 1:Q
            Pq = P; Pq.f0 = P.f0 + sfs(i)*randn(1, 2);
            if sfs(i) > 0, Pq.xr = P.xr .* (1 + 0.05*randn(1,2)); Pq.xi = P.xi .* (1 + 0.05*randn(1,2)); end
            Pl{q} = Pq;
        end
        W = square_weights(base, Q, Nt, 0);
        H = coupling_waveform(bins, fc, df, W, As, Ad, Pl, M-1);
        loss(t, i, :) = run_est({H, H, H}, Tnom, F, lamF, Np, evalb, snrs, slope, Q);
    end
end
R.a_spread_hz = sfs; R.a_loss_db = squeeze(mean(loss, 1));
fprintf('(a) element spread (MHz) %s -> loss dB @20 %s | @30 %s\n', mat2str(sfs/1e6), sprintf('%.2f ', R.a_loss_db(:,1)), sprintf('%.2f ', R.a_loss_db(:,2)));

%% (b) switching transients
rhos = [0 0.05 0.1 0.2 0.3]; loss = zeros(ntrial, numel(rhos), 2);
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    for i = 1:numel(rhos)
        W = square_weights(base, Q, Nt, rhos(i));
        H = coupling_waveform(bins, fc, df, W, As, Ad, P, M-1);
        loss(t, i, :) = run_est({H, H, H}, Tnom, F, lamF, Np, evalb, snrs, slope, Q);
    end
end
R.b_rise_fraction_of_slot = rhos; R.b_loss_db = squeeze(mean(loss, 1));
fprintf('(b) rise time (slots) %s -> loss dB @20 %s | @30 %s\n', mat2str(rhos), sprintf('%.2f ', R.b_loss_db(:,1)), sprintf('%.2f ', R.b_loss_db(:,2)));

%% (c) receiver CFO (same for pilots and data)
epss = [0 0.005 0.01 0.02 0.05]; loss = zeros(ntrial, numel(epss), 2);
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H0 = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
    for i = 1:numel(epss)
        H = ici_cfo(M, epss(i)) * H0;
        loss(t, i, :) = run_est({H, H, H}, Tnom, F, lamF, Np, evalb, snrs, slope, Q);
    end
end
R.c_cfo_subcarriers = epss; R.c_loss_db = squeeze(mean(loss, 1));
fprintf('(c) CFO (subcarriers) %s -> loss dB @20 %s | @30 %s\n', mat2str(epss), sprintf('%.2f ', R.c_loss_db(:,1)), sprintf('%.2f ', R.c_loss_db(:,2)));

%% (d) phase noise, new realisation per OFDM symbol
lws = [0 100 1e3 1e4]; loss = zeros(ntrial, numel(lws), 2);
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H0 = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
    for i = 1:numel(lws)
        Hs = cell(1, 3);
        for n = 1:3, Hs{n} = ici_pn(M, lws(i), df) * H0; end
        loss(t, i, :) = run_est(Hs, Tnom, F, lamF, Np, evalb, snrs, slope, Q);
    end
end
R.d_linewidth_hz = lws; R.d_loss_db = squeeze(mean(loss, 1));
fprintf('(d) phase-noise linewidth (Hz) %s -> loss dB @20 %s | @30 %s\n', mat2str(lws), sprintf('%.2f ', R.d_loss_db(:,1)), sprintf('%.2f ', R.d_loss_db(:,2)));
fid = fopen(fullfile(out, 'week6_hardware.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

%% ---------------------------------------------------------------------
function W = square_weights(base, Q, Nt, rho)
% State weights over one period; each state change is a linear ramp of rho slots
% centred on the nominal edge (rho = 0: ideal instant switching).
K = numel(base); t = ((0:Nt-1) + 0.5) / Nt * K;           % time in slots
w1 = double(base(floor(t) + 1) == 1);
if rho > 0
    for k = 0:K-1
        a = base(mod(k-1, K) + 1); b = base(k + 1);
        if a == b, continue; end
        d = t - k; d(d > K/2) = d(d > K/2) - K; d(d < -K/2) = d(d < -K/2) + K;
        in = abs(d) < rho/2;
        ramp = (d(in) + rho/2) / rho;                       % 0 -> 1 across the edge
        if b == 1, w1(in) = ramp; else, w1(in) = 1 - ramp; end
    end
end
W = zeros(Q, Nt, 2); W(:, :, 1) = repmat(w1, Q, 1); W(:, :, 2) = repmat(1 - w1, Q, 1);
end

function C = ici_cfo(M, eps)
k = (0:M-1).' - (0:M-1); x = eps - k;
C = exp(1j*pi*x*(M-1)/M) .* snc(x) ./ snc(x/M);
end

function C = ici_pn(M, lw, df)
% Receive-side ICI matrix of one OFDM symbol of Wiener phase noise (M samples).
if lw == 0, C = eye(M); return; end
Ts = 1/(df*M);
phi = cumsum(sqrt(2*pi*lw*Ts) * randn(M, 1)) + 2*pi*rand;
c = fft(exp(1j*phi)) / M;                                   % c(k+1): leakage by k bins
k = mod((0:M-1).' - (0:M-1), M);
C = c(k + 1);
end

function y = snc(x)
y = ones(size(x)); nz = abs(x) > 1e-12; y(nz) = sin(pi*x(nz)) ./ (pi*x(nz));
end

function l = run_est(Hs, Tnom, F, lamF, Np, evalb, snrs, slope, Q)
M = size(Hs{1}, 1); l = zeros(1, 2); Hd = Hs{3};
Pc = mean(sum(abs(Hd(:, evalb)).^2, 1));
for is = 1:2
    s2 = Pc / 10^(snrs(is)/10);
    S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
    Y = [Hs{1}*S(:,1), Hs{2}*S(:,2)] + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
    Hh = est_general(Y, S, {Tnom{1}, Tnom{2}}, ones(Q,1), F, lamF, s2, 8, 20);
    l(is) = (mmse_rate(Hd, Hd, s2, evalb) - mmse_rate(Hh, Hd, s2, evalb)) / slope(is);
end
end

function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
