%WEEK6_PHASENOISE Fair phase-noise test (replaces (d2) of week6_hardware_fix).
%   Continuous Wiener phase over pilot 1, pilot 2 and the data symbol, common phase
%   error removed per symbol. Three receivers on the phase-noisy data symbol:
%     genie_pn   knows the data symbol's exact phase-noise ICI (unknowable; upper bound)
%     genie_ch   knows the channel H0 exactly but not the random phase noise
%                (the fair reference: what any receiver could at best know)
%     ours       timing-only estimate from the two phase-noisy pilot symbols
%   Reported: generic OFDM phase-noise penalty (genie_pn vs no phase noise) and ours vs
%   genie_ch (our estimator's extra sensitivity). Bar: ours within 0.5 dB of genie_ch
%   at 20 dB for linewidths <= 1 kHz.

clear; rng(3600);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3; Tcp = 2.34e-6;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8; Np = 2;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
slope = [0.31 0.33]; snrs = [20 30]; ntrial = 10;
lws = [10 100 1e3 1e4];
rate = zeros(ntrial, numel(lws), 2, 4);                 % [nopn genie_pn genie_ch ours]
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H0 = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
    for i = 1:numel(lws)
        C = pn3(M, lws(i), df, Tcp);
        Hs = cell(1, 3);
        for n = 1:3, Cn = C{n} * exp(-1j*angle(C{n}(1,1))); Hs{n} = Cn * H0; end
        Pc = mean(sum(abs(H0(:, evalb)).^2, 1));
        for is = 1:2
            s2 = Pc / 10^(snrs(is)/10);
            S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
            Y = [Hs{1}*S(:,1), Hs{2}*S(:,2)] + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
            Hh = est_general(Y, S, T, ones(Q,1), F, lamF, s2, 8, 20);
            rate(t, i, is, :) = [mmse_rate(H0, H0, s2, evalb), mmse_rate(Hs{3}, Hs{3}, s2, evalb), ...
                                 mmse_rate(H0, Hs{3}, s2, evalb), mmse_rate(Hh, Hs{3}, s2, evalb)];
        end
    end
end
m = squeeze(mean(rate, 1));                               % lw x snr x 4
R = struct('linewidth_hz', lws, 'snr_db', snrs, 'rate_nopn_geniepn_geniech_ours', m);
for i = 1:numel(lws)
    for is = 1:2
        fprintf('lw %6.0f Hz, %d dB: no-PN %.2f | genie-PN %.2f | genie-channel %.2f | ours %.2f  -> generic penalty %.2f dB, ours vs genie-channel %.2f dB\n', ...
            lws(i), snrs(is), squeeze(m(i, is, :)), (m(i,is,1) - m(i,is,2))/slope(is), (m(i,is,3) - m(i,is,4))/slope(is));
    end
end
fid = fopen(fullfile(out, 'week6_phasenoise.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function C = pn3(M, lw, df, Tcp)
Ts = 1/(df*M); ncp = round(Tcp/Ts); ntot = 3*(M + ncp);
phi = cumsum(sqrt(2*pi*lw*Ts) * randn(ntot, 1));
C = cell(1, 3);
for n = 1:3
    s0 = (n-1)*(M + ncp) + ncp;
    c = fft(exp(1j*phi(s0 + (1:M)))) / M;
    k = mod((0:M-1).' - (0:M-1), M);
    C{n} = c(k + 1);
end
end
function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
