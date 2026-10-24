%WEEK12_STANDARD_CURVES Standard-format results with intermediate references (pre-registered 2026-10-24).
%   S1: NMSE, robust-MMSE rate and uncoded 16-QAM BER vs SNR 0..30 dB (shared path Np = 2, TDL-C Np = 4).
%   S2: element non-ideality continuum H_a = (1-a) H_ideal + a H_OpenRIS, a = 0..1, 20/30 dB, Np = 2.
%   S3: NMSE vs Np = 1..32 at 20 dB, shared path.
%   References between ours and the literature: genie (model-aware LMMSE with the TRUE element response),
%   mismatched element model (resonances shifted 5/10/20 MHz), CE-BEM with all harmonics at Np = 32.
%   Set parts (default 1:3) and ntrial (default 30) before calling to run a subset. Saves after every channel.

if ~exist('parts', 'var'), parts = 1:3; end
if ~exist('ntrial', 'var'), ntrial = 30; end
if ~exist('outname', 'var'), outname = 'week12_standard_curves'; end
here = fileparts(mfilename('fullpath')); out = fullfile(here, '..', 'results');
P = params_openris(); fc = 3.594e9; df = 30e3; K = 8;
M = 256; bins = -M/2:M/2-1; evalb = 65:192; Q = 16; fin = fc + bins*df;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
cfg.phys30 = struct('tw', 0.3e-6, 'L', 16, 'r', 8, 'tw2', 0.6e-6);
cfg.TDLC100 = struct('tw', 1e-6, 'L', ceil(M*df*1e-6) + 10, 'r', 16, 'tw2', 2e-6);
gam = @(Pm, a) deal((1-a)*1 + a*element_response(fin, 1, Pm), (1-a)*(-1) + a*element_response(fin, 2, Pm));
Pshift = @(d) setfield(P, 'f0', P.f0 + d); %#ok<SFLD>
Nd = 20; Mq = 16;
fp = fullfile(out, [outname '.json']);
if exist(fp, 'file'), R = jsondecode(fileread(fp)); else, R = struct(); end

%% ---------------- S1: vs SNR ----------------
if ismember(1, parts)
    snrs = 0:5:30;
    ests = {'onetap', 'ideal', 'mm20', 'mm10', 'mm5', 'genie', 'bem128', 'bemall_Np32', 'omp', 'ours'};
    chans = {'phys30', 'TDLC100'}; npc = [2 4];
    for ic = 1:numel(chans)
        cf = cfg.(chans{ic}); Np = npc(ic);
        [F, lamF] = dpss_basis(bins, df, cf.tw, cf.L);
        L2 = ceil(M*df*cf.tw2) + 12; [F2, lam2] = dpss_basis(bins, df, cf.tw2, L2);
        nm = nan(ntrial, numel(ests), numel(snrs)); rt = nm; be = nm; rp = nan(ntrial, numel(snrs)); bp = rp;
        for t = 1:ntrial
            rng(12100 + 100*ic + t);
            if ic == 1, [As, Ad] = phys_channel(Q, bins, df, fc, popt);
            else, [As, Ad] = tdl38901(Q, bins, df, fc, 'C', 100e-9); end
            H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
            Pc = mean(sum(abs(H(:, evalb)).^2, 1));
            Sp = qpsk(M, 32); W0 = randn(M, 32) + 1j*randn(M, 32);
            d = randi(Mq, M, Nd) - 1; X = qammod_(d, Mq); Wd = randn(M, Nd) + 1j*randn(M, Nd);
            for is = 1:numel(snrs)
                s2 = Pc / 10^(snrs(is)/10);
                Yall = H*Sp + sqrt(s2/2)*W0; Yd = H*X + sqrt(s2/2)*Wd;
                [rp(t, is), bp(t, is)] = metrics(H, H, s2, 0, evalb, Yd, d, Mq);
                S = Sp(:, 1:Np); Y = Yall(:, 1:Np);
                for ie = 1:numel(ests)
                    Hh = run_est(ests{ie}, Y, S, Yall, Sp, T, F, lamF, F2, lam2, s2, cf, P, Pshift, gam, bins, df, M, Q);
                    se_r = max(norm(Y - Hh*S, 'fro')^2/(M*Np) - s2, 0);
                    if strcmp(ests{ie}, 'bemall_Np32'), se_r = max(norm(Yall - Hh*Sp, 'fro')^2/(M*32) - s2, 0); end
                    nm(t, ie, is) = norm(Hh(:, evalb) - H(:, evalb), 'fro')^2 / norm(H(:, evalb), 'fro')^2;
                    [rt(t, ie, is), be(t, ie, is)] = metrics(Hh, H, s2, se_r, evalb, Yd, d, Mq);
                end
            end
            fprintf('S1 %s ch %2d | 20 dB rate %s | perfect %.2f\n', chans{ic}, t, sprintf('%.2f ', rt(t, :, 5)), rp(t, 5));
            R.S1.(chans{ic}) = struct('snr_db', snrs, 'Np', Np, 'estimators', {ests}, 'nmse', nm(1:t,:,:), ...
                'rate', rt(1:t,:,:), 'ber', be(1:t,:,:), 'rate_perfect', rp(1:t,:), 'ber_perfect', bp(1:t,:), ...
                'bits_per_channel', numel(evalb)*Nd*log2(Mq));
            save_json(fp, R);
        end
    end
end

%% ---------------- S2: element non-ideality continuum ----------------
if ismember(2, parts)
    alphas = [0 0.1 0.25 0.5 0.75 1]; snrs = [20 30]; Np = 2; cf = cfg.phys30;
    ests = {'onetap', 'ideal', 'genie', 'bem128', 'ours'};
    [F, lamF] = dpss_basis(bins, df, cf.tw, cf.L);
    L2 = ceil(M*df*cf.tw2) + 12; [F2, lam2] = dpss_basis(bins, df, cf.tw2, L2);
    rt = nan(ntrial, numel(ests), numel(alphas), numel(snrs)); nm = rt; rp = nan(ntrial, numel(alphas), numel(snrs));
    for t = 1:ntrial
        rng(12300 + t);
        [As, Ad] = phys_channel(Q, bins, df, fc, popt);
        Hr = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
        Hi = coupling_matrix(bins, fc, df, seqs, As, Ad, 'ideal', M-1);
        Sp = qpsk(M, Np); W0 = randn(M, Np) + 1j*randn(M, Np);
        for ia = 1:numel(alphas)
            a = alphas(ia); H = (1-a)*Hi + a*Hr;
            [g1, g2] = gam(P, a); Tg = {T{1} .* g1, -T{2} .* g2};
            Pc = mean(sum(abs(H(:, evalb)).^2, 1));
            for is = 1:numel(snrs)
                s2 = Pc / 10^(snrs(is)/10); Y = H*Sp + sqrt(s2/2)*W0;
                rp(t, ia, is) = metrics(H, H, s2, 0, evalb, [], [], Mq);
                for ie = 1:numel(ests)
                    switch ests{ie}
                        case 'onetap', Hh = diag(mean(Y ./ Sp, 2));
                        case 'ideal',  Hh = est_ideal(Y, Sp, T, F, lamF, s2);
                        case 'genie',  Hh = est_ideal(Y, Sp, Tg, F, lamF, s2);
                        case 'bem128', Hh = est_bem_cg(Y, Sp, F2, lam2, s2, -128:128);
                        case 'ours',   Hh = est_general(Y, Sp, T, ones(Q,1), F, lamF, s2, cf.r, 20);
                    end
                    se_r = max(norm(Y - Hh*Sp, 'fro')^2/(M*Np) - s2, 0);
                    nm(t, ie, ia, is) = norm(Hh(:, evalb) - H(:, evalb), 'fro')^2 / norm(H(:, evalb), 'fro')^2;
                    rt(t, ie, ia, is) = metrics(Hh, H, s2, se_r, evalb, [], [], Mq);
                end
            end
        end
        fprintf('S2 ch %2d | 20 dB gap (ideal, genie, ours) at a = 0: %.3f %.3f %.3f, a = 1: %.3f %.3f %.3f\n', t, ...
            rp(t,1,1) - rt(t,[2 3 5],1,1), rp(t,end,1) - rt(t,[2 3 5],end,1));
        R.S2 = struct('alpha', alphas, 'snr_db', snrs, 'Np', Np, 'estimators', {ests}, 'rate', rt(1:t,:,:,:), ...
            'nmse', nm(1:t,:,:,:), 'rate_perfect', rp(1:t,:,:));
        save_json(fp, R);
    end
end

%% ---------------- S3: NMSE vs pilots ----------------
if ismember(3, parts)
    Nps = [1 2 4 8 16 32]; snr = 20; cf = cfg.phys30;
    ests = {'ideal', 'mm10', 'genie', 'bem128', 'bemall', 'omp', 'ours'};
    [F, lamF] = dpss_basis(bins, df, cf.tw, cf.L);
    L2 = ceil(M*df*cf.tw2) + 12; [F2, lam2] = dpss_basis(bins, df, cf.tw2, L2);
    nm = nan(ntrial, numel(ests), numel(Nps)); rt = nm; rp = nan(ntrial, 1);
    for t = 1:ntrial
        rng(12500 + t);
        [As, Ad] = phys_channel(Q, bins, df, fc, popt);
        H = coupling_matrix(bins, fc, df, seqs, As, Ad, P, M-1);
        Pc = mean(sum(abs(H(:, evalb)).^2, 1)); s2 = Pc / 10^(snr/10);
        Sp = qpsk(M, 32); Yall = H*Sp + noise(M, 32, s2);
        rp(t) = metrics(H, H, s2, 0, evalb, [], [], Mq);
        for ip = 1:numel(Nps)
            Np = Nps(ip); S = Sp(:, 1:Np); Y = Yall(:, 1:Np);
            for ie = 1:numel(ests)
                if strcmp(ests{ie}, 'bemall') && Np < 8, continue; end
                try
                    Hh = run_est(ests{ie}, Y, S, Y, S, T, F, lamF, F2, lam2, s2, cf, P, Pshift, gam, bins, df, M, Q);
                catch err
                    fprintf('S3 %s Np=%d failed: %s\n', ests{ie}, Np, err.message); continue;
                end
                se_r = max(norm(Y - Hh*S, 'fro')^2/(M*Np) - s2, 0);
                nm(t, ie, ip) = norm(Hh(:, evalb) - H(:, evalb), 'fro')^2 / norm(H(:, evalb), 'fro')^2;
                rt(t, ie, ip) = metrics(Hh, H, s2, se_r, evalb, [], [], Mq);
            end
        end
        fprintf('S3 ch %2d | NMSE dB at Np=2: %s\n', t, sprintf('%.1f ', 10*log10(nm(t, :, 2))));
        R.S3 = struct('Np', Nps, 'snr_db', snr, 'estimators', {ests}, 'nmse', nm(1:t,:,:), 'rate', rt(1:t,:,:), ...
            'rate_perfect', rp(1:t));
        save_json(fp, R);
    end
end
fprintf('done\n');

%% ---------------------------------------------------------------------
function Hh = run_est(name, Y, S, Yall, Sall, T, F, lamF, F2, lam2, s2, cf, P, Pshift, gam, bins, df, M, Q)
switch name
    case 'onetap', Hh = diag(mean(Y ./ S, 2));
    case 'ideal',  Hh = est_ideal(Y, S, T, F, lamF, s2);
    case 'genie',  [g1, g2] = gam(P, 1);  Hh = est_ideal(Y, S, {T{1} .* g1, -T{2} .* g2}, F, lamF, s2);
    case 'mm5',    [g1, g2] = gam(Pshift(5e6), 1);  Hh = est_ideal(Y, S, {T{1} .* g1, -T{2} .* g2}, F, lamF, s2);
    case 'mm10',   [g1, g2] = gam(Pshift(10e6), 1); Hh = est_ideal(Y, S, {T{1} .* g1, -T{2} .* g2}, F, lamF, s2);
    case 'mm20',   [g1, g2] = gam(Pshift(20e6), 1); Hh = est_ideal(Y, S, {T{1} .* g1, -T{2} .* g2}, F, lamF, s2);
    case 'bem128', Hh = est_bem_cg(Y, S, F2, lam2, s2, -128:128);
    case 'bemall', Hh = est_bem_cg(Y, S, F2, lam2, s2, -(M-1):(M-1));
    case 'bemall_Np32', Hh = est_bem_cg(Yall, Sall, F2, lam2, s2, -(M-1):(M-1));
    case 'omp',    Hh = est_omp(Y, S, s2, -16:16, 0:1/(2*M*df):cf.tw2, bins, df);
    case 'ours',   Hh = est_general(Y, S, T, ones(Q,1), F, lamF, s2, cf.r, 20);
end
end
function [r, ber] = metrics(Hh, H, s2, se2, evalb, Yd, d, Mq)
% robust-MMSE rate (known-G detector, as week10) and uncoded 16-QAM BER after unbiased MMSE
M = size(H, 1);
W = (Hh'*Hh + (s2 + se2)*eye(M)) \ Hh'; G = W*H; g = diag(G);
it = sum(abs(G).^2, 2) - abs(g).^2 + s2*sum(abs(W).^2, 2);
r = mean(log2(1 + abs(g(evalb)).^2 ./ it(evalb)));
ber = nan;
if ~isempty(Yd)
    Z = (W*Yd) ./ diag(W*Hh);
    dh = qamdemod_(Z(evalb, :), Mq);
    ber = mean(de2bi_(dh, log2(Mq)) ~= de2bi_(d(evalb, :), log2(Mq)), 'all');
end
end
function save_json(fp, R)
fid = fopen(fp, 'w'); fprintf(fid, '%s', jsonencode(R)); fclose(fid);
end
function Hh = est_ideal(Y, S, T, F, lam, s2)
% model-aware LMMSE: H = (T{1} - T{2}) o (F C F^T); T{s} may carry the element response as column scaling
[M, Np] = size(S); L = size(F, 2); y = Y(:); Td = T{1} - T{2};
Phi = zeros(M*Np, L*L);
for p = 1:Np
    Gs = Td * (S(:,p) .* F);
    for i = 1:L, Phi((p-1)*M+(1:M), (i-1)*L+(1:L)) = F(:, i) .* Gs; end
end
w = kron(lam(:), lam(:));
sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
Dw = sqrt(sc*w);
c = Dw .* (((Dw .* (Phi'*Phi)) .* Dw.' + s2*eye(numel(w))) \ (Dw .* (Phi'*y)));
Hh = Td .* (F*reshape(c, L, L).'*F.');
end
function Hh = est_omp(Y, S, s2, hs, taus, bins, df)
[M, Np] = size(S); y = Y(:);
E = exp(-1j*2*pi*bins(:)*df*taus); nh = numel(hs); nt = numel(taus);
D = zeros(M*Np, nh*nt);
for k = 1:nh, for p = 1:Np, D((p-1)*M+(1:M), (k-1)*nt+(1:nt)) = shr(S(:,p) .* E, hs(k)); end, end
Dn = D ./ vecnorm(D); res = y; sel = []; target = s2 * numel(y);
for it = 1:300
    [~, j] = max(abs(Dn' * res)); sel(end+1) = j; %#ok<AGROW>
    c = D(:, sel) \ y; res = y - D(:, sel) * c;
    if real(res'*res) <= target, break; end
end
Hh = zeros(M);
for ii = 1:numel(sel)
    k = floor((sel(ii)-1)/nt) + 1; ti = sel(ii) - (k-1)*nt;
    Hh = Hh + diag(ones(M-abs(hs(k)),1), -hs(k)) * diag(c(ii) * E(:, ti));
end
end
function Z = shr(X, h)
Z = zeros(size(X)); M = size(X, 1);
if h >= 0, Z(1+h:M, :) = X(1:M-h, :); else, Z(1:M+h, :) = X(1-h:M, :); end
end
function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
function W = noise(M, n, s2), W = sqrt(s2/2) * (randn(M,n) + 1j*randn(M,n)); end
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
function b = de2bi_(d, nb)
m = 2^(nb/2);
ir = mod(d, m); iq = floor(d / m);
g = @(v) bitxor(v, bitshift(v, -1));
b = false([size(d), nb]);
gi = g(ir); gq = g(iq);
for k = 1:nb/2
    b(:,:,k) = bitget(gi, k); b(:,:,nb/2+k) = bitget(gq, k);
end
end
