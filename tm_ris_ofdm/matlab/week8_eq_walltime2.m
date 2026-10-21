%WEEK8_EQ_WALLTIME2 (L5, POST-HOC follow-up to week8_eq_walltime, which FAILED its bar)
%   The pre-registered run timed fast_coupling, which loops over elements (2 Q small FFT
%   calls per product). With identical sequences all elements share ONE kernel per state,
%   so a product is H x = sum_q a_q o IFFT( sum_s Tf_s .* FFT(b_{q,s} o x) ):
%   one batched forward FFT per state over all (element, symbol) columns, the two states
%   summed in the frequency domain, one batched inverse FFT. Same arithmetic model,
%   fewer, larger FFT calls. Also 10 CG iterations (week7: LSQR-10 = full MMSE).
%   Reported separately from the bar; dense numbers re-measured in the same run.

clear; rng(4100);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3; Q = 16; K = 8; Nd = 12; Mq = 16; snr_db = 20;
base = [ones(1, 6), 2*ones(1, 2)];
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
Ms = [1024 2048 3276 6552]; nrep = 3; its = [10 25];
R = struct('M', Ms, 'Nd', Nd, 'cg_iters', its, 'snr_db', snr_db, 'note', 'post-hoc, batched shared-kernel operator');
td = nan(numel(Ms), 1); tc = nan(numel(Ms), numel(its)); serd = nan(1, numel(Ms)); serc = nan(numel(its), numel(Ms));
for im = 1:numel(Ms)
    M = Ms(im); bins = -M/2:M/2-1; evalb = round(M/4)+1:round(3*M/4);
    f = fc + bins(:)*df;
    Gam = [element_response(f, 1, P), element_response(f, 2, P)];
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    Aout = Ad.'; B = {(As .* Gam(:,1).').', (As .* Gam(:,2).').'};
    op = batched_op(bins, base, Aout, B);
    d = randi(Mq, M, Nd) - 1; X = qammod_(d, Mq);
    Y0 = op.H(X);
    Pc = mean(abs(Y0(:)).^2) / mean(abs(X(:)).^2); s2 = Pc / 10^(snr_db/10);
    Y = Y0 + sqrt(s2/2)*(randn(M,Nd) + 1j*randn(M,Nd));
    if M <= 3276
        T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
             coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
        tt = zeros(1, nrep);
        for k = 1:nrep
            t0 = tic; H = T{1} .* (Aout*B{1}.') + T{2} .* (Aout*B{2}.');
            Zd = (H'*H + s2*eye(M)) \ (H'*Y); tt(k) = toc(t0); clear H
        end
        td(im) = median(tt); serd(im) = mean(qamdemod_(Zd(evalb,:), Mq) ~= d(evalb,:), 'all');
        clear T
    end
    for ii = 1:numel(its)
        tt = zeros(1, nrep);
        for k = 1:nrep, t0 = tic; Zc = cg_mmse(op, Y, s2, its(ii)); tt(k) = toc(t0); end
        tc(im, ii) = median(tt); serc(ii, im) = mean(qamdemod_(Zc(evalb,:), Mq) ~= d(evalb,:), 'all');
    end
    fprintf('M = %4d: dense %.2f s | CG-10 %.2f s, CG-25 %.2f s | SER dense %.4f CG-10 %.4f CG-25 %.4f\n', ...
        M, td(im), tc(im, :), serd(im), serc(:, im));
end
R.dense_s = td.'; R.cg_s = tc; R.ser_dense = serd; R.ser_cg = serc;
fid = fopen(fullfile(out, 'week8_eq_walltime2.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function op = batched_op(bins, base, Aout, B)
% H = sum_s T_s o (Aout B_s^T), shared Toeplitz kernel per state (identical sequences).
M = numel(bins); K = numel(base); Q = size(Aout, 2);
h = (-(M-1):(M-1)).'; env = exp(-1j*pi*h/K)/K .* snc(h/K);
Nf = 3*M; Tf = zeros(Nf, 2); Tfh = Tf;
for s = 1:2
    ks = find(base == s) - 1; t = env .* sum(exp(-1j*2*pi*h*ks/K), 2);
    Tf(:, s) = fft(t, Nf); Tfh(:, s) = fft(conj(flipud(t)), Nf);
end
op.H = @(x) fwd(x); op.Ht = @(y) adj(y);
    function y = fwd(x)
        n = size(x, 2); acc = 0;
        for s = 1:2
            Z = reshape(B{s}, M, Q, 1) .* reshape(x, M, 1, n);          % M x Q x n
            acc = acc + Tf(:, s) .* fft(reshape(Z, M, Q*n), Nf);
        end
        c = ifft(acc); c = reshape(c(M:2*M-1, :), M, Q, n);
        y = reshape(sum(Aout .* c, 2), M, n);
    end
    function x = adj(y)
        n = size(y, 2);
        Z = conj(Aout) .* reshape(y, M, 1, n);
        Yf = fft(reshape(Z, M, Q*n), Nf); x = 0;
        for s = 1:2
            c = ifft(Tfh(:, s) .* Yf); c = reshape(c(M:2*M-1, :), M, Q, n);
            x = x + reshape(sum(conj(B{s}) .* c, 2), M, n);
        end
    end
end
function y = snc(x), y = ones(size(x)); nz = x ~= 0; y(nz) = sin(pi*x(nz)) ./ (pi*x(nz)); end
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
