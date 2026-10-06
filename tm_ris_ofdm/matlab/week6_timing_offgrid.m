%WEEK6_TIMING_OFFGRID Follow-up to week6_timing.
%   (c2) Blind offset recovery with offsets OFF the scan grid: true offset uniform in
%        [-0.5, 0.5] slot; coarse scan (step 0.025) then golden-section refinement of
%        the pilot residual on the bracketing interval. No truth used.
%        Same bar as before: within 0.3 dB of the known-timing case at 20 dB.
%   (b2) Fine jitter levels, to set the hardware requirement: rms 0.001, 0.002, 0.005 slot.
%   S1, D = 0.75, K = 8, Q = 16, M = 256, Np = 2, physical channel.

clear; rng(3500);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8; Np = 2;
base = [ones(1, 6), 2*ones(1, 2)]; st = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
[F, lamF] = dpss_basis(bins, df, 0.3e-6, 16);
E0 = (0:K)/K;
Tnom = {coupling_edges(bins, fc, df, E0, base, ones(1,M), ones(1,M), 'state1', M-1), ...
        coupling_edges(bins, fc, df, E0, base, ones(1,M), ones(1,M), 'state2', M-1)};
slope = [0.31 0.33]; snrs = [20 30]; R = struct();

%% (c2) off-grid blind recovery
ntrial = 12; grid = -0.5:0.025:0.5;
res = zeros(ntrial, 2, 2); offs = zeros(ntrial, 2);           % [blind known], [true est]
for t = 1:ntrial
    d0 = rand - 0.5;
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    H = coupling_edges(bins, fc, df, repmat(((0:K) + d0)/K, Q, 1), st, As, Ad, P, M-1);
    Pc = mean(sum(abs(H(:, evalb)).^2, 1));
    for is = 1:2
        s2 = Pc / 10^(snrs(is)/10);
        S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
        Y = H*S + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
        f = @(d) resid(d, Y, S, Tnom, bins, K, F, lamF, s2, Q, 8);
        r = arrayfun(f, grid); [~, i0] = min(r);
        a = grid(max(i0-1, 1)); b = grid(min(i0+1, numel(grid)));
        dh = golden(f, a, b, 14);
        Hb = est_general(Y, S, shiftk(Tnom, bins, dh, K), ones(Q,1), F, lamF, s2, 8, 20);
        Hk = est_general(Y, S, shiftk(Tnom, bins, d0, K), ones(Q,1), F, lamF, s2, 8, 20);
        rp = mmse_rate(H, H, s2, evalb);
        res(t, is, :) = (rp - [mmse_rate(Hb, H, s2, evalb), mmse_rate(Hk, H, s2, evalb)]) / slope(is);
        if is == 1, offs(t, :) = [d0, dh]; end
    end
end
R.c2_loss_db_blind_known = squeeze(mean(res, 1));            % snr x [blind known]
R.c2_offset_error_rms = sqrt(mean((offs(:,1) - offs(:,2)).^2));
R.c2_offsets = offs;
fprintf('(c2) off-grid: offset error rms %.4f slot | loss @20 blind/known %.2f/%.2f | @30 %.2f/%.2f\n', ...
    R.c2_offset_error_rms, R.c2_loss_db_blind_known(1,:), R.c2_loss_db_blind_known(2,:));

%% (b2) fine jitter
jit = [0.001 0.002 0.005]; ntrial = 10; loss = zeros(ntrial, numel(jit), 2);
for t = 1:ntrial
    [As, Ad] = phys_channel(Q, bins, df, fc, popt);
    for ij = 1:numel(jit)
        Hs = cell(1, 3);
        for n = 1:3
            E = repmat((0:K)/K, Q, 1);
            E(:, 2:K) = E(:, 2:K) + jit(ij)/K * randn(Q, K-1); E = sort(E, 2);
            Hs{n} = coupling_edges(bins, fc, df, E, st, As, Ad, P, M-1);
        end
        Pc = mean(sum(abs(Hs{3}(:, evalb)).^2, 1));
        for is = 1:2
            s2 = Pc / 10^(snrs(is)/10);
            S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
            Y = [Hs{1}*S(:,1), Hs{2}*S(:,2)] + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
            Hh = est_general(Y, S, Tnom, ones(Q,1), F, lamF, s2, 8, 20);
            loss(t, ij, is) = (mmse_rate(Hs{3}, Hs{3}, s2, evalb) - mmse_rate(Hh, Hs{3}, s2, evalb)) / slope(is);
        end
    end
end
R.b2_jitter_rms_slots = jit; R.b2_jitter_ns = jit * (1/df) / K * 1e9; R.b2_loss_db = squeeze(mean(loss, 1));
for ij = 1:numel(jit)
    fprintf('(b2) jitter %.3f slot (%.1f ns): loss @20 %.2f | @30 %.2f\n', jit(ij), R.b2_jitter_ns(ij), R.b2_loss_db(ij,:));
end
fid = fopen(fullfile(out, 'week6_timing_offgrid.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function r = resid(d, Y, S, Tnom, bins, K, F, lamF, s2, Q, it)
Hc = est_general(Y, S, shiftk(Tnom, bins, d, K), ones(Q,1), F, lamF, s2, 8, it);
r = norm(Y - Hc*S, 'fro');
end
function T = shiftk(Tnom, bins, d, K)
D = exp(-1j*2*pi*bins(:)*d/K);
T = {D .* Tnom{1} .* D', D .* Tnom{2} .* D'};
end
function x = golden(f, a, b, n)
g = (sqrt(5) - 1)/2; c = b - g*(b - a); d = a + g*(b - a); fc = f(c); fd = f(d);
for i = 1:n
    if fc < fd, b = d; d = c; fd = fc; c = b - g*(b - a); fc = f(c);
    else, a = c; c = d; fc = fd; d = a + g*(b - a); fd = f(d); end
end
x = (a + b)/2;
end
function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
