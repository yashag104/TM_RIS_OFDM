%WEEK11_EQ_AMORTISED Review H2 + H3 (pre-registered 2026-10-22). Run on an otherwise idle machine.
%   H2  Dense MMSE with amortisation: factor (Hh'Hh + s2 I) once (Cholesky) per channel estimate, then N
%       triangular solves; CG: N symbols batched as columns, stop at tolerance. Wall-clock for
%       N in {1, 12, 40, 120} at M in {1024, 2048, 3276} (dense) and M = 6552 (CG only); break-even N.
%   H3  CG tolerance in {1e-1, 3e-2, 1e-2, 3e-3, 1e-3}: SER and iterations vs the dense solve at each M, and
%       the final relative residual (the CG equaliser is "exact products, iterative solve").
%   Batched shared-kernel operator (week8_eq_walltime2), 16-QAM, 20 dB, median of 3 repetitions.
clear; rng(11800);
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; Q = 16; K = 8; Mq = 16; snr_db = 20;
base = [ones(1, 6), 2*ones(1, 2)];
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
Ms = [1024 2048 3276 6552]; Ns = [1 12 40 120]; tols = [1e-1 3e-2 1e-2 3e-3 1e-3]; nrep = 3;
R = struct('M', Ms, 'N', Ns, 'tol', tols, 'snr_db', snr_db, 'threads', maxNumCompThreads);
for im = 1:numel(Ms)
    M = Ms(im); bins = -M/2:M/2-1; evalb = round(M/4)+1:round(3*M/4); f = fc + bins(:)*df;
    Gam = [element_response(f, 1, P), element_response(f, 2, P)];
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    Aout = Ad.'; B = {(As .* Gam(:,1).').', (As .* Gam(:,2).').'};
    op = batched_op(bins, base, Aout, B);
    d = randi(Mq, M, max(Ns)) - 1; X = qammod_(d, Mq); Y0 = op.H(X);
    s2 = (mean(abs(Y0(:)).^2)/mean(abs(X(:)).^2)) / 10^(snr_db/10);
    Y = Y0 + sqrt(s2/2)*(randn(M, max(Ns)) + 1j*randn(M, max(Ns)));
    td_f = nan; td_s = nan(1, numel(Ns)); serd = nan;
    if M <= 3276
        T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
             coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
        tt = zeros(1, nrep);
        for k = 1:nrep
            t0 = tic; H = T{1} .* (Aout*B{1}.') + T{2} .* (Aout*B{2}.'); A = H'*H + s2*eye(M); Rc = chol(A); tt(k) = toc(t0);
        end
        td_f = median(tt);
        for in = 1:numel(Ns)
            tt = zeros(1, nrep);
            for k = 1:nrep, t0 = tic; Zd = Rc \ (Rc' \ (H'*Y(:, 1:Ns(in)))); tt(k) = toc(t0); end
            td_s(in) = median(tt);
        end
        serd = mean(qamdemod_(Zd(evalb, :), Mq) ~= d(evalb, 1:Ns(end)), 'all');
        clear T H A Rc
    end
    tc = nan(numel(tols), numel(Ns)); itc = nan(numel(tols), 1); serc = nan(numel(tols), 1); resc = serc;
    for it = 1:numel(tols)
        for in = 1:numel(Ns)
            tt = zeros(1, nrep);
            for k = 1:nrep, t0 = tic; [Zc, ni, rr] = cg_tol(op, Y(:, 1:Ns(in)), s2, tols(it), 200); tt(k) = toc(t0); end
            tc(it, in) = median(tt);
        end
        itc(it) = ni; resc(it) = rr; serc(it) = mean(qamdemod_(Zc(evalb, :), Mq) ~= d(evalb, 1:Ns(end)), 'all');
    end
    key = sprintf('M%d', M);
    R.(key) = struct('dense_factor_s', td_f, 'dense_solve_s', td_s, 'dense_total_s', td_f + td_s, 'ser_dense', serd, ...
                     'cg_s', tc, 'cg_iters', itc, 'cg_final_relres', resc, 'ser_cg', serc);
    fprintf('M = %d: dense factor %.2f s, solves %s s | SER dense %.4f\n', M, td_f, sprintf('%.3f ', td_s), serd);
    for it = 1:numel(tols)
        fprintf('   tol %.0e: CG %s s, %d it, relres %.1e, SER %.4f\n', tols(it), sprintf('%.2f ', tc(it, :)), itc(it), resc(it), serc(it));
    end
end
fid = fopen(fullfile(out, 'week11_eq_amortised.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function [z, it, rr] = cg_tol(op, Y, s2, tol, maxit)
b = op.Ht(Y); z = zeros(size(b)); r = b; p = r; rs = sum(conj(r).*r, 1); r0 = rs;
for it = 1:maxit
    Ap = op.Ht(op.H(p)) + s2*p; a = rs ./ sum(conj(p).*Ap, 1);
    z = z + a.*p; r = r - a.*Ap; rn = sum(conj(r).*r, 1);
    if all(sqrt(rn ./ r0) < tol), break; end
    p = r + (rn./rs).*p; rs = rn;
end
rr = max(sqrt(rn ./ r0));
end
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
