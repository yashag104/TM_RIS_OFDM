%WEEK6_DIRECTLINK (#12) Unblocked direct path, two RIS panels, co-channel interference.
%   (a) direct BS -> user path at power rd relative to the RIS path,
%       rd in {-inf, -20, -10, 0, +10} dB. Receivers: ours ignoring the direct path
%       (est_general), ours with the direct path (est_general_direct), one-tap, perfect.
%       Note one-tap benefits from a strong direct path (it is diagonal).
%   (b) two TM-RIS panels at different positions (independent path sets), each with
%       its own switching pattern (G = 2 groups); rank r = 8 vs r = 16.
%   (c) co-channel interference: an OFDM interferer through its own frequency-selective
%       channel at SIR in {inf, 20, 10, 0} dB, unknown to the receiver (treated as noise).
%   S1 pattern unless stated, Q = 16, M = 256, Np = 2, 10 trials.

clear; rng(3000);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8; Np = 2;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
[Fd, ~] = dpss_basis(bins, df, 0.6e-6, 20);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
snrs = [20 30]; ntrial = 10; R = struct();

%% (a) direct path
rds = [-inf -20 -10 0 10];
names = {'onetap', 'ours_nodirect', 'ours_direct', 'perfect'};
rate = zeros(ntrial, numel(rds), numel(names), 2);
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    Hr = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
    hd = tdl_resp(bins, df, 60e-9);                    % unit-power direct response
    Pr = mean(sum(abs(Hr(:, evalb)).^2, 1));
    for i = 1:numel(rds)
        H = Hr;
        if isfinite(rds(i)), H = Hr + diag(sqrt(Pr * 10^(rds(i)/10)) * hd); end
        Pc = mean(sum(abs(H(:, evalb)).^2, 1));
        for is = 1:2
            s2 = Pc / 10^(snrs(is)/10);
            S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
            Y = H*S + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
            g1 = mean(Y ./ S, 2); w = conj(g1) ./ abs(g1).^2; dH = diag(H);
            sinr = abs(w.*dH).^2 ./ (abs(w).^2 .* (sum(abs(H).^2, 2) - abs(dH).^2 + s2));
            rate(t, i, 1, is) = mean(log2(1 + sinr(evalb)));
            rate(t, i, 2, is) = mmse_rate(est_general(Y, S, T, ones(Q,1), F, lamF, s2, 8, 20), H, s2, evalb);
            rate(t, i, 3, is) = mmse_rate(est_general_direct(Y, S, T, ones(Q,1), F, lamF, s2, 8, 20, Fd), H, s2, evalb);
            rate(t, i, 4, is) = mmse_rate(H, H, s2, evalb);
        end
    end
end
R.a_direct_rel_db = rds; R.a_names = names; R.a_rate = squeeze(mean(rate, 1));
for i = 1:numel(rds)
    fprintf('(a) direct %+5.0f dB: rate @20 %s | @30 %s  [onetap noDirect direct perfect]\n', rds(i), ...
        sprintf('%.2f ', squeeze(R.a_rate(i,:,1))), sprintf('%.2f ', squeeze(R.a_rate(i,:,2))));
end

%% (b) two RIS panels with independent channels and different patterns
rate = zeros(ntrial, 4, 2);                             % [r8, r16, perfect, single-panel-model r8]
for t = 1:ntrial
    base2 = circshift(base, 3);
    [As1, Ad1] = phys_channel(Q, bins, df, fc, popt); [As2, Ad2] = phys_channel(Q, bins, df, fc, popt);
    H = coupling_matrix(bins, fc, df, seqs, As1, Ad1, P, M-1) + ...
        coupling_matrix(bins, fc, df, repmat(base2, Q, 1), As2, Ad2, P, M-1);
    T2 = {T{1}, T{2}; coupling_matrix(bins, fc, df, base2, ones(1,M), ones(1,M), 'state1', M-1), ...
                      coupling_matrix(bins, fc, df, base2, ones(1,M), ones(1,M), 'state2', M-1)};
    Pc = mean(sum(abs(H(:, evalb)).^2, 1));
    for is = 1:2
        s2 = Pc / 10^(snrs(is)/10);
        S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
        Y = H*S + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
        rate(t, 1, is) = mmse_rate(est_general(Y, S, T2, [], F, lamF, s2, 8, 20), H, s2, evalb);
        rate(t, 2, is) = mmse_rate(est_general(Y, S, T2, [], F, lamF, s2, 16, 20), H, s2, evalb);
        rate(t, 3, is) = mmse_rate(H, H, s2, evalb);
        rate(t, 4, is) = mmse_rate(est_general(Y, S, T, [], F, lamF, s2, 8, 20), H, s2, evalb);
    end
end
R.b_rate = squeeze(mean(rate, 1));                      % 4 x snr
fprintf('(b) two panels: rate @20 r8 %.2f r16 %.2f perfect %.2f (one-panel model %.2f) | @30 r8 %.2f r16 %.2f perfect %.2f (%.2f)\n', ...
    R.b_rate(:,1), R.b_rate(:,2));

%% (c) co-channel interference
sirs = [inf 20 10 0]; rate = zeros(ntrial, numel(sirs), 2, 2);   % [ours, perfect]
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
    Pc = mean(sum(abs(H(:, evalb)).^2, 1));
    hi = tdl_resp(bins, df, 100e-9);
    for i = 1:numel(sirs)
        for is = 1:2
            s2 = Pc / 10^(snrs(is)/10);
            pint = 0; if isfinite(sirs(i)), pint = Pc / 10^(sirs(i)/10); end
            S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
            Xi = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
            Y = H*S + sqrt(pint) * (hi(:) .* Xi) + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
            Hh = est_general(Y, S, T, ones(Q,1), F, lamF, s2 + pint, 8, 20);
            % rate with interference treated as extra white noise of the same power
            rate(t, i, 1, is) = mmse_rate(Hh, H, s2 + pint, evalb);
            rate(t, i, 2, is) = mmse_rate(H, H, s2 + pint, evalb);
        end
    end
end
R.c_sir_db = sirs; R.c_rate = squeeze(mean(rate, 1));
for i = 1:numel(sirs)
    fprintf('(c) SIR %4.0f dB: rate ours/perfect @20 %.2f/%.2f | @30 %.2f/%.2f\n', sirs(i), ...
        R.c_rate(i,1,1), R.c_rate(i,2,1), R.c_rate(i,1,2), R.c_rate(i,2,2));
end
fid = fopen(fullfile(out, 'week6_directlink.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function h = tdl_resp(bins, df, rms)
tau = (0:11) * rms / 2; pdp = exp(-tau / rms); pdp = pdp / sum(pdp);
g = sqrt(pdp/2) .* (randn(1,12) + 1j*randn(1,12));
h = (g * exp(-1j*2*pi*tau.' * (bins*df))).';
end

function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
