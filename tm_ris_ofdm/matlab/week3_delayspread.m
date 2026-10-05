%WEEK3_DELAYSPREAD Robustness to delay spreads longer than assumed (reviewer check).
%   Every earlier test used 30 ns RMS per hop, well inside the estimator's
%   assumed support tau_max = 0.3 us. Here the truth uses 30, 100 and 200 ns RMS
%   (exponential tails then exceed 0.3 us), with the estimator settings UNCHANGED.
%   Timing-only (model-free) estimator, S1, Np = 2 pilots + 2 DD iterations.
%   Pre-registered bar: loss <= 1 dB at 20 dB for 100 ns RMS. 200 ns reported as is.

clear; rng(1900);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8; D = 0.75; Np = 2; Nd = 6; ntrial = 20;
base = [ones(1, round(D*K)), 2*ones(1, K - round(D*K))];
seqs = repmat(base, Q, 1);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
Bmf = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
       coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
rms_list = [30e-9 100e-9 200e-9]; snrs = [20 30]; qam = [16 64]; slope = [0.31 0.33];
R = struct('rms_s', rms_list, 'snr_db', snrs);
for ir = 1:numel(rms_list)
    popt = struct('P', 5, 'K_db', 10, 'rms', rms_list(ir));
    loss = zeros(ntrial, 2);
    for t = 1:ntrial
        [As, Ad] = phys_channel(Q, bins, df, fc, popt);
        H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
        Pc = mean(sum(abs(H(:, evalb)).^2, 1));
        for is = 1:2
            s2 = Pc / 10^(snrs(is)/10); Mq = qam(is);
            S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
            Yp = H*S + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
            Xd = qammod_(randi(Mq, M, Nd) - 1, Mq);
            Yd = H*Xd + sqrt(s2/2)*(randn(M,Nd) + 1j*randn(M,Nd));
            Hh = est_lmmse(Yp, S, Bmf, F, lamF, s2, 'linear', 15);
            for it = 1:2
                Xh = qamdec_(((Hh'*Hh + s2*eye(M)) \ Hh') * Yd, Mq);
                Hh = est_lmmse([Yp Yd], [S Xh], Bmf, F, lamF, s2, 'linear', 15);
            end
            loss(t, is) = (mmse_rate(H, H, s2, evalb) - mmse_rate(Hh, H, s2, evalb)) / slope(is);
        end
    end
    key = sprintf('rms%dns', round(rms_list(ir)*1e9));
    R.(key) = struct('loss_db_mean', mean(loss, 1), 'loss_db_ci95', 1.96*std(loss, 0, 1)/sqrt(ntrial));
    fprintf('%s: loss dB at 20/30 dB = %.2f +/- %.2f / %.2f +/- %.2f\n', key, mean(loss(:,1)), ...
        1.96*std(loss(:,1))/sqrt(ntrial), mean(loss(:,2)), 1.96*std(loss(:,2))/sqrt(ntrial));
end
fid = fopen(fullfile(out, 'week3_delayspread.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
function x = qammod_(d, Mq)
m = sqrt(Mq); lv = -(m-1):2:(m-1);
x = (lv(mod(d, m) + 1) + 1j*lv(floor(d / m) + 1)) / sqrt(2*(Mq-1)/3);
end
function x = qamdec_(z, Mq)
m = sqrt(Mq); sc = sqrt(2*(Mq-1)/3); lv = -(m-1):2:(m-1);
q = @(v) lv(min(max(round((v + m - 1)/2), 0), m-1) + 1);
x = reshape((q(real(z)*sc) + 1j*q(imag(z)*sc)) / sc, size(z));
end
