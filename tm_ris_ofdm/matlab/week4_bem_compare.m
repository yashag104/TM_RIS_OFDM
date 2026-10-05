%WEEK4_BEM_COMPARE Our estimator vs the literature's estimators, same pilots and seeds.
%   Estimators (all from Np pilot OFDM symbols, physical channel, S1, M = 256):
%     onetap    conventional OFDM, per-subcarrier LS + one-tap equalizer
%     toeplitz  ideal-model mixing H(mbar,m) = t(mbar - m) (as assumed in Tao et al.)
%     bemQ      CE-BEM (complex-exponential basis expansion, the standard estimator
%               for intra-symbol time variation): H = sum_{|h|<=Q} S_h diag(F2 c_h),
%               delay basis F2 covers both hops (tau = 0.6 us), LMMSE with DPSS prior
%               (a strong version of BEM, given the same prior as ours). Q = 4, 8, 12.
%     ours      timing-only: H = sum_s T_s o (F C_s F.'), LMMSE with DPSS prior
%     perfect   true H
%   Score: rate of the MMSE equalizer built from each estimate, applied to the true H.
%   20 trials, 95% CIs; Np in {1, 2, 4, 8}; SNR 20 and 30 dB.

clear; rng(2100);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8; D = 0.75;
base = [ones(1, round(D*K)), 2*ones(1, K - round(D*K))];
seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);                     % per hop (ours)
L2 = ceil(M*df*0.6e-6) + 12; [F2, lam2] = dpss_basis(bins, df, 0.6e-6, L2);   % two-hop (BEM)
Bmf = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
       coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
Qb = [4 8 12]; Nps = [1 2 4 8]; snrs = [20 30]; ntrial = 20;
names = [{'onetap', 'toeplitz'}, arrayfun(@(q) sprintf('bem%d', q), Qb, 'uni', 0), {'ours', 'perfect'}];
rate = zeros(ntrial, numel(names), numel(Nps), numel(snrs));
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
    Pc = mean(sum(abs(H(:, evalb)).^2, 1));
    for is = 1:numel(snrs)
        s2 = Pc / 10^(snrs(is)/10);
        for ip = 1:numel(Nps)
            Np = Nps(ip);
            S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
            Y = H*S + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
            for e = 1:numel(names)
                switch names{e}
                    case 'onetap'
                        g = mean(Y ./ S, 2); w = conj(g) ./ abs(g).^2; dH = diag(H);
                        sinr = abs(w.*dH).^2 ./ (abs(w).^2 .* (sum(abs(H).^2, 2) - abs(dH).^2 + s2));
                        rate(t, e, ip, is) = mean(log2(1 + sinr(evalb))); continue
                    case 'toeplitz'
                        hs = -12:12; Phi = zeros(M*Np, numel(hs));
                        for k = 1:numel(hs), Phi(:, k) = reshape(sh(S, hs(k)), [], 1); end
                        c = (Phi'*Phi + 1e-3*real(trace(Phi'*Phi))/numel(hs)*eye(numel(hs))) \ (Phi'*Y(:));
                        Hh = zeros(M); for k = 1:numel(hs), Hh = Hh + c(k)*shiftmat(M, hs(k)); end
                    case 'ours'
                        Hh = est_lmmse(Y, S, Bmf, F, lamF, s2, 'linear', 15);
                    case 'perfect'
                        Hh = H;
                    otherwise                                          % bemQ
                        q = str2double(names{e}(4:end));
                        Hh = est_bem(Y, S, F2, lam2, s2, -q:q);
                end
                rate(t, e, ip, is) = mmse_rate(Hh, H, s2, evalb);
            end
        end
    end
end
R = struct('names', {names}, 'Np', Nps, 'snr_db', snrs, 'ntrial', ntrial);
R.rate_mean = squeeze(mean(rate, 1));                % names x Np x snr
R.rate_ci95 = squeeze(1.96*std(rate, 0, 1)/sqrt(ntrial));
for is = 1:numel(snrs)
    fprintf('SNR %d dB, rate (bit/s/Hz) for Np = %s:\n', snrs(is), mat2str(Nps));
    for e = 1:numel(names)
        fprintf('  %-9s %s\n', names{e}, sprintf('%5.2f+-%4.2f  ', [R.rate_mean(e,:,is); R.rate_ci95(e,:,is)]));
    end
end
fid = fopen(fullfile(out, 'week4_bem_compare.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function Hh = est_bem(Y, S, F2, lam, s2, hs)
% CE-BEM in the frequency domain: H = sum_h S_h diag(F2 c_h), LMMSE with DPSS prior.
[M, Np] = size(S); L = size(F2, 2); Nh = numel(hs); y = Y(:);
Phi = zeros(M*Np, Nh*L);
for k = 1:Nh
    for p = 1:Np
        Phi((p-1)*M+(1:M), (k-1)*L+(1:L)) = sh(S(:,p) .* F2, hs(k));
    end
end
w = repmat(lam(:), Nh, 1);
sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
Dw = sqrt(sc*w);
c = Dw .* (((Dw .* (Phi'*Phi)) .* Dw.' + s2*eye(numel(w))) \ (Dw .* (Phi'*y)));
Hh = zeros(M);
for k = 1:Nh, Hh = Hh + shiftmat(M, hs(k)) * diag(F2 * c((k-1)*L+(1:L))); end
end

function Z = sh(X, h)
% Shift rows: output subcarrier m+h receives input m (truncated to the band).
Z = zeros(size(X)); M = size(X, 1);
if h >= 0, Z(1+h:M, :) = X(1:M-h, :); else, Z(1:M+h, :) = X(1-h:M, :); end
end

function Sm = shiftmat(M, h)
Sm = diag(ones(M - abs(h), 1), -h);
end

function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
