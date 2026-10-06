%WEEK7_COMPLEXITY_FIX (#10) Cheaper channel updates for a moving user.
%   A user moves at 30 km/h (only the RIS -> user hop changes). The channel is
%   re-estimated every N = 20 symbols from pilots at the start of each block, for 8
%   blocks; scored on a data symbol in the middle of each block. Update variants:
%     cold      full estimate from scratch every block (20 ALS iterations, rank 8)
%     warm2     start from the previous block's estimate, 2 iterations
%     uonly     V (static BS -> RIS hop + element response) frozen from block 0;
%               only Lam and U updated, 1 iteration
%     uonly_1p  same with ONE pilot symbol instead of two
%     cold_r6   full estimate, rank 6
%   Reported: gap to perfect (bit/s/Hz), measured MATLAB time per update, and the
%   analytic operation count per update (normal equations + regressors).
%   S1, Q = 16, K = 8, M = 256, physical channel, 8 trials, 20 and 30 dB.

clear; rng(4000);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3; Tsym = 1/df + 2.34e-6;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8; N = 20; B = 8;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16); L = size(F, 2);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
fD = 30/3.6 * fc / 3e8;
snrs = [20 30]; ntrial = 8;
names = {'cold', 'warm2', 'uonly', 'uonly_1p', 'cold_r6'};
gap = zeros(ntrial, numel(names), 2); tim = zeros(ntrial, numel(names), 2);
for t = 1:ntrial
    [As, Adp, nu] = moving_channel(Q, bins, df, fc, popt, fD);
    Hat = @(n) channel_at(n, bins, fc, df, seqs, As, Adp, nu, P, M, Tsym);
    H0 = Hat(0); Pc = mean(sum(abs(H0(:, evalb)).^2, 1));
    for is = 1:2
        s2 = Pc / 10^(snrs(is)/10);
        st = cell(1, numel(names));
        for b = 0:B-1
            n0 = b*N; Hp = {Hat(n0), Hat(n0 + 1)}; Hd = Hat(n0 + N/2);
            S = qpsk(M, 2); Y = [Hp{1}*S(:,1), Hp{2}*S(:,2)] + noise(M, 2, s2);
            rp = mmse_rate(Hd, Hd, s2, evalb);
            for e = 1:numel(names)
                tic;
                switch names{e}
                    case 'cold',    [Hh, st{e}] = est_general_ws(Y, S, T, F, lamF, s2, 8, 20, [], false);
                    case 'cold_r6', [Hh, st{e}] = est_general_ws(Y, S, T, F, lamF, s2, 6, 20, [], false);
                    case 'warm2'
                        if b == 0, [Hh, st{e}] = est_general_ws(Y, S, T, F, lamF, s2, 8, 20, [], false);
                        else, [Hh, st{e}] = est_general_ws(Y, S, T, F, lamF, s2, 8, 2, st{e}, false); end
                    case {'uonly', 'uonly_1p'}
                        if b == 0
                            [Hh, st{e}] = est_general_ws(Y, S, T, F, lamF, s2, 8, 20, [], false);
                        else
                            if strcmp(names{e}, 'uonly_1p'), Yu = Y(:,1); Su = S(:,1); else, Yu = Y; Su = S; end
                            [Hh, st{e}] = est_general_ws(Yu, Su, T, F, lamF, s2, 8, 1, st{e}, true);
                        end
                end
                el = toc;
                if b > 0                                   % score the refresh blocks only
                    gap(t, e, is) = gap(t, e, is) + (rp - mmse_rate(Hh, Hd, s2, evalb)) / (B-1);
                    tim(t, e, is) = tim(t, e, is) + el / (B-1);
                end
            end
        end
    end
    fprintf('trial %d done\n', t);
end
% analytic operations per refresh (Np pilots, rank r, L basis functions, G = 1)
ops = @(Np, r, it, cold, fixV) M*Np*( (cold)*(2*L*L)^2 ) + (cold)*(2*L*L)^3/3 + ...
      it*( M*Np*(r*r)^2 + (r*r)^3/3 + M*Np*(L*r)^2 + (L*r)^3/3 + (~fixV)*(M*Np*(2*L*r)^2 + (2*L*r)^3/3) );
opsv = [ops(2,8,20,1,0), ops(2,8,2,0,0), ops(2,8,1,0,1), ops(1,8,1,0,1), ops(2,6,20,1,0)];
R = struct('names', {names}, 'snr_db', snrs, 'speed_kmh', 30, 'refresh_symbols', N, ...
    'gap_bits', squeeze(mean(gap, 1)), 'time_s', squeeze(mean(tim, 1)), 'ops_per_update', opsv, ...
    'ops_reduction_vs_cold', opsv(1) ./ opsv);
for e = 1:numel(names)
    fprintf('%-9s gap @20 %.2f @30 %.2f | time %.3f s | ops %.2e (%.0fx fewer than cold)\n', names{e}, ...
        R.gap_bits(e,1), R.gap_bits(e,2), mean(R.time_s(e,:)), opsv(e), R.ops_reduction_vs_cold(e));
end
fid = fopen(fullfile(out, 'week7_complexity_fix.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

%% ---------------------------------------------------------------------
function [As, Adp, nu] = moving_channel(Q, bins, df, fc, popt, fD)
lam = 3e8 / fc; k = 2*pi/lam; n = round(sqrt(Q));
[gx, gz] = meshgrid(0:n-1, 0:n-1); r = [gx(:), zeros(Q,1), gz(:)] * lam/2;
[As, ~] = phys_channel(Q, bins, df, fc, popt);
P = popt.P; Kr = 10^(popt.K_db/10);
tau = [0, sort(-popt.rms*log(rand(1, P)))];
pw = [Kr/(Kr+1), (1/(Kr+1)) * exp(-tau(2:end)/popt.rms) / sum(exp(-tau(2:end)/popt.rms))];
g = sqrt(pw) .* exp(1j*2*pi*rand(1, P+1));
u = randn(P+1, 3); u(:,2) = abs(u(:,2)); u = u ./ vecnorm(u, 2, 2);
ph = exp(-1j*k*(r*u.'));
Adp = zeros(Q, numel(bins), P+1);
for l = 1:P+1, Adp(:, :, l) = (ph(:, l) * g(l)) * exp(-1j*2*pi*tau(l)*(bins(:).'*df)); end
nu = fD * cos(2*pi*rand(1, P+1));
end

function H = channel_at(n, bins, fc, df, seqs, As, Adp, nu, P, M, Tsym)
H = zeros(M);
for l = 1:size(Adp, 3)
    Ad = Adp(:, :, l) * exp(1j*2*pi*nu(l)*n*Tsym);
    H = H + ici_matrix(M, nu(l)/df) * coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
end
end

function C = ici_matrix(M, eps)
k = (0:M-1).' - (0:M-1); x = eps - k;
C = exp(1j*pi*x*(M-1)/M) .* snc(x) ./ snc(x/M);
end
function y = snc(x), y = ones(size(x)); nz = abs(x) > 1e-12; y(nz) = sin(pi*x(nz)) ./ (pi*x(nz)); end
function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
function W = noise(M, n, s2), W = sqrt(s2/2) * (randn(M,n) + 1j*randn(M,n)); end
function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
