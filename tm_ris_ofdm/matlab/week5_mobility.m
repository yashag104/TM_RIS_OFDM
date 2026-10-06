%WEEK5_MOBILITY Coherence time: channel ageing and intra-symbol Doppler.
%   Only the RIS -> user hop moves (BS and RIS are fixed). Each of its paths l has
%   Doppler nu_l = fD cos(phi_l), fD = v fc / c. Two effects are modelled:
%     (1) inter-symbol ageing: path l rotates by exp(j 2 pi nu_l n Tsym) over symbols n;
%     (2) intra-symbol Doppler: within a symbol, path l is a carrier offset
%         eps_l = nu_l / df, i.e. an extra ICI matrix C(eps_l) on the receive side.
%   The estimate is made once from 2 pilot symbols at n = 0, 1 (timing-only,
%   est_general), then used for data symbols n = 2 ... N. Score: MMSE rate with the
%   stale estimate on the true channel of symbol n, vs perfect knowledge of that symbol.
%   Speeds 3, 30, 120 km/h at 3.594 GHz, 30 kHz SCS (Tsym = 35.7 us with normal CP).
%   Output: rate vs n, number of symbols before a 1 dB loss, net efficiency with the
%   pilots repeated at that interval, and the intra-symbol Doppler penalty alone.

clear; rng(2700);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3; Tsym = 1/df + 2.34e-6;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8; D = 0.75; Np = 2;
base = [ones(1, round(D*K)), 2*ones(1, K - round(D*K))];
seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
speeds_kmh = [3 30 120]; snrs = [20 30]; ntrial = 10;
nvec = [2 5 10 20 40 80 150 300];
R = struct('speeds_kmh', speeds_kmh, 'snr_db', snrs, 'n', nvec, 'Tsym_s', Tsym);
for iv = 1:numel(speeds_kmh)
    fD = speeds_kmh(iv)/3.6 * fc / 3e8;
    rate_est = zeros(ntrial, numel(nvec), numel(snrs)); rate_perf = rate_est; rate_static = zeros(ntrial, numel(snrs));
    for t = 1:ntrial
        [As, Adp, nu] = moving_channel(Q, bins, df, fc, popt, fD);   % Adp: Q x M x paths
        H_at = @(n, intra) channel_at(n, intra, bins, fc, df, seqs, As, Adp, nu, P, M, Tsym);
        H0 = H_at(0, true); H1 = H_at(1, true);
        Pc = mean(sum(abs(H0(:, evalb)).^2, 1));
        for is = 1:numel(snrs)
            s2 = Pc / 10^(snrs(is)/10);
            S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
            Y = [H0*S(:,1), H1*S(:,2)] + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
            Hh = est_general(Y, S, {T{1}, T{2}}, ones(Q,1), F, lamF, s2, 8, 20);
            rate_static(t, is) = mmse_rate(Hh, H0, s2, evalb);
            for in = 1:numel(nvec)
                Hn = H_at(nvec(in), true);
                rate_est(t, in, is) = mmse_rate(Hh, Hn, s2, evalb);
                rate_perf(t, in, is) = mmse_rate(Hn, Hn, s2, evalb);
            end
        end
        % intra-symbol Doppler alone: perfect knowledge of the Doppler-free vs true channel
        if t == 1
            Hnd = H_at(0, false); Hd = H_at(0, true); s2 = Pc / 10^(30/10);
            R.(sprintf('v%d', speeds_kmh(iv))).intra_doppler_loss_bits_30dB = ...
                mmse_rate(Hd, Hd, s2, evalb) - mmse_rate(Hnd, Hd, s2, evalb);
            R.(sprintf('v%d', speeds_kmh(iv))).intra_doppler_rel_power_db = ...
                10*log10(norm(Hd - Hnd, 'fro')^2 / norm(Hd, 'fro')^2);
        end
    end
    key = sprintf('v%d', speeds_kmh(iv));
    R.(key).fD_hz = fD;
    R.(key).coherence_time_s = 0.423 / fD;
    R.(key).coherence_symbols = 0.423 / fD / Tsym;
    R.(key).rate_est = squeeze(mean(rate_est, 1));          % n x snr
    R.(key).rate_perfect = squeeze(mean(rate_perf, 1));
    R.(key).rate_at_estimation = mean(rate_static, 1);
    gap = squeeze(mean(rate_perf - rate_est, 1));           % bits, n x snr
    R.(key).gap_bits = gap;
    fprintf('%3d km/h: fD %.0f Hz, Tc %.2f ms (%.0f symbols) | intra-symbol Doppler %.1f dB rel., %.3f bit loss @30dB\n', ...
        speeds_kmh(iv), fD, 1e3*0.423/fD, 0.423/fD/Tsym, R.(key).intra_doppler_rel_power_db, R.(key).intra_doppler_loss_bits_30dB);
    for is = 1:numel(snrs)
        fprintf('   %d dB: gap (bits) at n = %s : %s\n', snrs(is), mat2str(nvec), sprintf('%.2f ', gap(:, is)));
    end
end
fid = fopen(fullfile(out, 'week5_mobility.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function [As, Adp, nu] = moving_channel(Q, bins, df, fc, popt, fD)
% Physical channel with per-path components kept for the moving (RIS -> user) hop.
lam = 3e8 / fc; k = 2*pi/lam; n = round(sqrt(Q));
[gx, gz] = meshgrid(0:n-1, 0:n-1); r = [gx(:), zeros(Q,1), gz(:)] * lam/2;
[As, ~] = phys_channel(Q, bins, df, fc, popt);                  % static BS -> RIS hop
P = popt.P; Kr = 10^(popt.K_db/10);
tau = [0, sort(-popt.rms*log(rand(1, P)))];
pw = [Kr/(Kr+1), (1/(Kr+1)) * exp(-tau(2:end)/popt.rms) / sum(exp(-tau(2:end)/popt.rms))];
g = sqrt(pw) .* exp(1j*2*pi*rand(1, P+1));
u = randn(P+1, 3); u(:,2) = abs(u(:,2)); u = u ./ vecnorm(u, 2, 2);
ph = exp(-1j*k*(r*u.'));                                       % Q x paths
Adp = zeros(Q, numel(bins), P+1);
for l = 1:P+1
    Adp(:, :, l) = (ph(:, l) * g(l)) * exp(-1j*2*pi*tau(l)*(bins(:).'*df));
end
nu = fD * cos(2*pi*rand(1, P+1));                              % per-path Doppler (Hz)
end

function H = channel_at(n, intra, bins, fc, df, seqs, As, Adp, nu, P, M, Tsym)
% Coupling at OFDM symbol n; path l rotated by its Doppler phase, and (if intra)
% followed by the receive-side ICI matrix of its carrier offset eps_l = nu_l/df.
H = zeros(M);
for l = 1:size(Adp, 3)
    Ad = Adp(:, :, l) * exp(1j*2*pi*nu(l)*n*Tsym);
    Hl = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
    if intra, Hl = ici_matrix(M, nu(l)/df) * Hl; end
    H = H + Hl;
end
end

function C = ici_matrix(M, eps)
% Receive-side ICI of a carrier offset eps (in subcarriers) over one OFDM symbol.
k = (0:M-1).' - (0:M-1);                                       % mbar - m
x = eps - k;                                                   % (m - mbar) + eps: same sign as the
                                                               % inter-symbol rotation exp(+j2pi nu t)
C = exp(1j*pi*x*(M-1)/M) .* sinc_(x) ./ sinc_(x/M);
C(abs(eps) < 1e-12 & k ~= 0) = 0;
end

function y = sinc_(x)
y = ones(size(x)); nz = abs(x) > 1e-12; y(nz) = sin(pi*x(nz)) ./ (pi*x(nz));
end

function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
