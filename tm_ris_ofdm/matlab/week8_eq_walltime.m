%WEEK8_EQ_WALLTIME (L5) Wall-clock time: matrix-free CG-MMSE vs dense MMSE.
%   Both equalise Nd = 12 data symbols of one slot with the SAME channel, given in the
%   factored form the receiver holds: H = sum_s T_s o (A_out B_s^T), inner dimension
%   Q = 16 (the exact physical model; the low-rank estimate has the same form, r <= Q).
%     dense  build H (T_s precomputed once: they depend on timing only, not timed),
%            then (H'H + s2 I) \ (H'Y)
%     cg     25 CG iterations on the normal equations, FFT operator (fast_coupling)
%   Bar: CG faster in wall-clock at M = 3276, with the same SER (16-QAM, 20 dB).
%   M = 512, 1024, 2048, 3276; median of 3 timed repetitions; S1, physical channel.

clear; rng(4100);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3; Q = 16; K = 8; Nd = 12; nit = 25; Mq = 16; snr_db = 20;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
Ms = [512 1024 2048 3276]; nrep = 3;
R = struct('M', Ms, 'Nd', Nd, 'cg_iters', nit, 'snr_db', snr_db, 'threads', maxNumCompThreads);
td = zeros(numel(Ms), nrep); tc = td; tb = td; serd = zeros(1, numel(Ms)); serc = serd;
for im = 1:numel(Ms)
    M = Ms(im); bins = -M/2:M/2-1; evalb = round(M/4)+1:round(3*M/4);
    f = fc + bins(:)*df;
    Gam = [element_response(f, 1, P), element_response(f, 2, P)];
    T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
         coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    Aout = Ad.'; B = {(As .* Gam(:,1).').', (As .* Gam(:,2).').'};          % M x Q each
    op = fast_coupling(bins, fc, df, seqs, As, Ad, P);
    Hx = T{1} .* (Aout*B{1}.') + T{2} .* (Aout*B{2}.');
    Pc = mean(sum(abs(Hx(:, evalb)).^2, 1)); s2 = Pc / 10^(snr_db/10);
    d = randi(Mq, M, Nd) - 1; X = qammod_(d, Mq);
    Y = Hx*X + sqrt(s2/2)*(randn(M,Nd) + 1j*randn(M,Nd));
    clear Hx
    for k = 1:nrep
        t0 = tic; H = T{1} .* (Aout*B{1}.') + T{2} .* (Aout*B{2}.'); tb(im, k) = toc(t0);
        Zd = (H'*H + s2*eye(M)) \ (H'*Y); td(im, k) = toc(t0);
        clear H
        t0 = tic; Zc = cg_mmse(op, Y, s2, nit); tc(im, k) = toc(t0);
    end
    serd(im) = mean(qamdemod_(Zd(evalb,:), Mq) ~= d(evalb,:), 'all');
    serc(im) = mean(qamdemod_(Zc(evalb,:), Mq) ~= d(evalb,:), 'all');
    fprintf('M = %4d: dense %.2f s (build %.2f) | CG %.2f s | speed-up %.1fx | SER dense %.4f CG %.4f\n', ...
        M, median(td(im,:)), median(tb(im,:)), median(tc(im,:)), median(td(im,:))/median(tc(im,:)), serd(im), serc(im));
    clear T
end
R.dense_s = median(td, 2).'; R.dense_build_s = median(tb, 2).'; R.cg_s = median(tc, 2).';
R.speedup = R.dense_s ./ R.cg_s; R.ser_dense = serd; R.ser_cg = serc;
fid = fopen(fullfile(out, 'week8_eq_walltime.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function z = cg_mmse(op, Y, s2, nit)
b = op.Ht(Y); z = zeros(size(b)); r = b; p = r; rs = sum(conj(r).*r, 1);
for it = 1:nit
    Ap = op.Ht(op.H(p)) + s2*p;
    a = rs ./ sum(conj(p).*Ap, 1);
    z = z + a.*p; r = r - a.*Ap;
    rn = sum(conj(r).*r, 1);
    p = r + (rn./rs).*p; rs = rn;
end
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
