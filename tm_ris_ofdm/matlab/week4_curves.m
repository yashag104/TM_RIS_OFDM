%WEEK4_CURVES Standard communication-paper curves: NMSE vs SNR and BER vs SNR.
%   Physical channel, S1 (identical sequences, duty 0.75, K = 8), Q = 16, M = 256,
%   Np = 2 pilot OFDM symbols, SNR 0:2:30 dB, 30 channel realisations.
%   Estimators: ours (timing-only LMMSE), CE-BEM (Q = 12, same prior), Toeplitz model.
%   Receivers (16-QAM, Gray, uncoded): one-tap; MMSE with the Toeplitz / BEM / our
%   estimate; ours + 1 decision-directed pass; MMSE with the true channel.
%   Also stores constellation samples at 20 dB for one realisation.
%   Output: ../results/week4_curves.json

clear; rng(2400);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8; D = 0.75; Np = 2; Nd = 40; Mq = 16;
base = [ones(1, round(D*K)), 2*ones(1, K - round(D*K))];
seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
L2 = ceil(M*df*0.6e-6) + 12; [F2, lam2] = dpss_basis(bins, df, 0.6e-6, L2);
Bmf = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
       coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
snrs = 0:2:30; ntrial = 30;
est_names = {'toeplitz', 'bem', 'ours'};
rx_names = {'onetap', 'toeplitz', 'bem', 'ours', 'ours_dd', 'perfect'};
nmse = zeros(ntrial, numel(est_names), numel(snrs));
berr = zeros(ntrial, numel(rx_names), numel(snrs)); nbits = 0;
const = struct();
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
    Pc = mean(sum(abs(H(:, evalb)).^2, 1));
    d = randi(Mq, M, Nd) - 1; X = qammod_(d, Mq); bref = bits_(d(evalb,:), Mq);
    for is = 1:numel(snrs)
        s2 = Pc / 10^(snrs(is)/10);
        S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
        Yp = H*S + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
        Yd = H*X + sqrt(s2/2)*(randn(M,Nd) + 1j*randn(M,Nd));
        Hh = cell(1, 3);
        hs = -12:12; Phi = zeros(M*Np, numel(hs));
        for k = 1:numel(hs), Phi(:, k) = reshape(shr(S, hs(k)), [], 1); end
        c = (Phi'*Phi + 1e-3*real(trace(Phi'*Phi))/numel(hs)*eye(numel(hs))) \ (Phi'*Yp(:));
        Hh{1} = zeros(M); for k = 1:numel(hs), Hh{1} = Hh{1} + c(k)*diag(ones(M-abs(hs(k)),1), -hs(k)); end
        Hh{2} = est_bem(Yp, S, F2, lam2, s2, -12:12);
        Hh{3} = est_lmmse(Yp, S, Bmf, F, lamF, s2, 'linear', 15);
        for e = 1:3
            E = Hh{e}(evalb, evalb) - H(evalb, evalb);
            nmse(t, e, is) = norm(E, 'fro')^2 / norm(H(evalb, evalb), 'fro')^2;
        end
        Z = cell(1, numel(rx_names));
        g = mean(Yp ./ S, 2); Z{1} = Yd ./ g;                       % conventional: per-subcarrier LS
        for e = 1:3, Z{1+e} = mmse_eq(Hh{e}, Yd, s2); end
        Xh = qamdec_(Z{4}, Mq);
        Hdd = est_lmmse([Yp Yd], [S Xh], Bmf, F, lamF, s2, 'linear', 15);
        Z{5} = mmse_eq(Hdd, Yd, s2);
        Z{6} = mmse_eq(H, Yd, s2);
        for r = 1:numel(rx_names)
            bh = bits_(qamdemod_(Z{r}(evalb,:), Mq), Mq);
            berr(t, r, is) = sum(bh(:) ~= bref(:));
        end
        if t == 1 && snrs(is) == 20
            const.onetap = cplx(Z{1}(evalb, 1:8)); const.ours = cplx(Z{4}(evalb, 1:8));
            const.perfect = cplx(Z{6}(evalb, 1:8));
        end
    end
    nbits = nbits + numel(bref);
    fprintf('trial %d done\n', t);
end
R = struct('snr_db', snrs, 'ntrial', ntrial, 'Np', Np, 'qam', Mq, 'bits_per_point', nbits, ...
    'est_names', {est_names}, 'rx_names', {rx_names});
R.nmse_db_mean = 10*log10(squeeze(mean(nmse, 1)));                    % est x snr
R.nmse_db_ci95 = squeeze(1.96*std(10*log10(nmse), 0, 1)/sqrt(ntrial));
R.ber = squeeze(sum(berr, 1)) / nbits;                                % rx x snr
R.constellation_20dB = const;
fid = fopen(fullfile(out, 'week4_curves.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);
disp(R.nmse_db_mean); disp(R.ber);

function s = cplx(z)
s = struct('re', real(z(:)).', 'im', imag(z(:)).');
end
function z = mmse_eq(Hh, Y, s2)
M = size(Hh, 1); z = ((Hh'*Hh + s2*eye(M)) \ Hh') * Y;
end
function Hh = est_bem(Y, S, F2, lam, s2, hs)
[M, Np] = size(S); L = size(F2, 2); Nh = numel(hs); y = Y(:);
Phi = zeros(M*Np, Nh*L);
for k = 1:Nh, for p = 1:Np
    Phi((p-1)*M+(1:M), (k-1)*L+(1:L)) = shr(S(:,p) .* F2, hs(k));
end, end
w = repmat(lam(:), Nh, 1);
sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
Dw = sqrt(sc*w);
c = Dw .* (((Dw .* (Phi'*Phi)) .* Dw.' + s2*eye(numel(w))) \ (Dw .* (Phi'*y)));
Hh = zeros(M);
for k = 1:Nh, Hh = Hh + diag(ones(M-abs(hs(k)),1), -hs(k)) * diag(F2 * c((k-1)*L+(1:L))); end
end
function Z = shr(X, h)
Z = zeros(size(X)); M = size(X, 1);
if h >= 0, Z(1+h:M, :) = X(1:M-h, :); else, Z(1:M+h, :) = X(1-h:M, :); end
end
function x = qammod_(d, Mq)
m = sqrt(Mq); lv = -(m-1):2:(m-1);
x = (lv(mod(d, m) + 1) + 1j*lv(floor(d / m) + 1)) / sqrt(2*(Mq-1)/3);
end
function d = qamdemod_(z, Mq)
m = sqrt(Mq); sc = sqrt(2*(Mq-1)/3);
ir = min(max(round((real(z)*sc + m - 1)/2), 0), m-1);
iq = min(max(round((imag(z)*sc + m - 1)/2), 0), m-1);
d = ir + m*iq;
end
function x = qamdec_(z, Mq)
x = qammod_(qamdemod_(z, Mq), Mq);
end
function b = bits_(d, Mq)
% Gray bits per axis, log2(Mq) bits per symbol.
m = sqrt(Mq); nb = log2(m);
g = @(v) bitxor(v, bitshift(v, -1));
gi = g(mod(d, m)); gq = g(floor(d / m));
b = zeros([size(d), 2*nb]);
for k = 1:nb, b(:,:,k) = bitget(gi, k); b(:,:,nb+k) = bitget(gq, k); end
end
