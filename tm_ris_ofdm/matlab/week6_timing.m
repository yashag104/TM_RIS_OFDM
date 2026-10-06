%WEEK6_TIMING (#2) What if the receiver's knowledge of the switching timing is imperfect?
%   Truth: switching edges shifted / jittered relative to the ideal k/K grid
%   (COUPLING_EDGES, exact Fourier integral). Receiver: timing-only estimator
%   (est_general) built from the NOMINAL timing, then MMSE.
%     (a) common offset delta (all edges, all elements), in fractions of a slot
%     (b) random jitter, independent per edge, per element and per OFDM symbol
%         (pilot 1, pilot 2 and data each see their own edges)
%     (c) blind recovery of (a): the receiver scans delta_hat on a grid, builds the
%         shifted kernels exactly (a time shift is a pair of diagonal phase ramps),
%         and keeps the delta_hat with the smallest pilot residual (no truth used).
%   S1, D = 0.75, K = 8, Q = 16, M = 256, Np = 2, physical channel, 10 trials.
%   Pre-registered bars: (a)/(b) report the tolerance at which loss reaches 1 dB at
%   20 dB; (c) must bring a 0.2-slot offset back within 0.3 dB of the known-timing case.

clear; rng(2800);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8; Np = 2;
base = [ones(1, 6), 2*ones(1, 2)];
st = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
Tnom = kernels(bins, fc, df, base, 0, M);
slope = [0.31 0.33]; snrs = [20 30]; ntrial = 10;
R = struct();

%% (a) common offset
offs = [0 0.02 0.05 0.1 0.2 0.3 0.5];
loss = zeros(ntrial, numel(offs), 2);
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    for io = 1:numel(offs)
        E = repmat(((0:K) + offs(io)) / K, Q, 1);
        H = coupling_edges(bins, fc, df, E, st, As, Ad, P, M-1);
        loss(t, io, :) = run_est(H, H, H, Tnom, F, lamF, Np, evalb, snrs, slope, Q);
    end
end
R.a_offset_slots = offs; R.a_loss_db = squeeze(mean(loss, 1)); R.a_ci = squeeze(1.96*std(loss,0,1)/sqrt(ntrial));
fprintf('(a) offset (slots) %s\n    loss dB @20: %s\n    loss dB @30: %s\n', mat2str(offs), ...
    sprintf('%.2f ', R.a_loss_db(:,1)), sprintf('%.2f ', R.a_loss_db(:,2)));

%% (b) per-edge jitter, new edges for every OFDM symbol
jit = [0 0.01 0.02 0.05 0.1];
loss = zeros(ntrial, numel(jit), 2);
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    for ij = 1:numel(jit)
        Hs = cell(1, 3);
        for n = 1:3                                  % pilot 1, pilot 2, data symbol
            E = repmat((0:K)/K, Q, 1);
            E(:, 2:K) = E(:, 2:K) + jit(ij)/K * randn(Q, K-1);   % inner edges jitter
            E = sort(E, 2);
            Hs{n} = coupling_edges(bins, fc, df, E, st, As, Ad, P, M-1);
        end
        loss(t, ij, :) = run_est(Hs{1}, Hs{2}, Hs{3}, Tnom, F, lamF, Np, evalb, snrs, slope, Q);
    end
end
R.b_jitter_slots_rms = jit; R.b_loss_db = squeeze(mean(loss, 1)); R.b_ci = squeeze(1.96*std(loss,0,1)/sqrt(ntrial));
fprintf('(b) jitter rms (slots) %s\n    loss dB @20: %s\n    loss dB @30: %s\n', mat2str(jit), ...
    sprintf('%.2f ', R.b_loss_db(:,1)), sprintf('%.2f ', R.b_loss_db(:,2)));

%% (c) blind offset recovery
grid = -0.5:0.025:0.5;
true_off = [0.05 0.2 0.4];
lossc = zeros(ntrial, numel(true_off), 2, 3);          % [nominal, blind, known]
est_off = zeros(ntrial, numel(true_off));
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    for io = 1:numel(true_off)
        E = repmat(((0:K) + true_off(io)) / K, Q, 1);
        H = coupling_edges(bins, fc, df, E, st, As, Ad, P, M-1);
        Pc = mean(sum(abs(H(:, evalb)).^2, 1));
        for is = 1:2
            s2 = Pc / 10^(snrs(is)/10);
            S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
            Y = H*S + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
            rp = mmse_rate(H, H, s2, evalb);
            H0 = est_general(Y, S, {Tnom{1}, Tnom{2}}, ones(Q,1), F, lamF, s2, 8, 20);
            % blind scan: kernels for each candidate offset (exact: diagonal phase ramps)
            best = inf;
            for dh = grid
                Tc = shift_kernels(Tnom, bins, dh, K);
                Hc = est_general(Y, S, {Tc{1}, Tc{2}}, ones(Q,1), F, lamF, s2, 8, 8);
                res = norm(Y - Hc*S, 'fro');
                if res < best, best = res; Hb = Hc; dbest = dh; end
            end
            Tk = shift_kernels(Tnom, bins, true_off(io), K);
            Hk = est_general(Y, S, {Tk{1}, Tk{2}}, ones(Q,1), F, lamF, s2, 8, 20);
            lossc(t, io, is, :) = (rp - [mmse_rate(H0, H, s2, evalb), mmse_rate(Hb, H, s2, evalb), ...
                                         mmse_rate(Hk, H, s2, evalb)]) / slope(is);
            if is == 1, est_off(t, io) = dbest; end
        end
    end
end
R.c_true_offset = true_off; R.c_loss_db = squeeze(mean(lossc, 1));   % offset x snr x [nominal blind known]
R.c_est_offset_mean = mean(est_off, 1); R.c_est_offset_std = std(est_off, 0, 1);
for io = 1:numel(true_off)
    fprintf('(c) true offset %.2f: estimated %.3f +/- %.3f | loss dB @20 nominal/blind/known %.2f/%.2f/%.2f | @30 %.2f/%.2f/%.2f\n', ...
        true_off(io), R.c_est_offset_mean(io), R.c_est_offset_std(io), squeeze(R.c_loss_db(io,1,:)), squeeze(R.c_loss_db(io,2,:)));
end
fid = fopen(fullfile(out, 'week6_timing.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

%% ---------------------------------------------------------------------
function T = kernels(bins, fc, df, base, off, M)
K = numel(base); E = ((0:K) + off) / K;
T = {coupling_edges(bins, fc, df, E, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_edges(bins, fc, df, E, base, ones(1,M), ones(1,M), 'state2', M-1)};
end

function T = shift_kernels(Tnom, bins, d, K)
% A time shift of d slots multiplies harmonic h by exp(-j 2 pi h d / K), h = mbar - m.
D = exp(-1j*2*pi*bins(:)*d/K);
T = {D .* Tnom{1} .* D', D .* Tnom{2} .* D'};
end

function l = run_est(H1, H2, Hd, Tnom, F, lamF, Np, evalb, snrs, slope, Q)
% Pilots through H1, H2 (one per pilot symbol); data through Hd. Loss vs perfect Hd.
M = size(H1, 1); l = zeros(1, 2);
Pc = mean(sum(abs(Hd(:, evalb)).^2, 1));
for is = 1:2
    s2 = Pc / 10^(snrs(is)/10);
    S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
    Y = [H1*S(:,1), H2*S(:,2)] + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
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
