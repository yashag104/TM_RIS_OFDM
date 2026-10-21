%WEEK8_TWOSCALE_300NS_C (L1c) As week8_twoscale_300ns_b (L1b) with Nf = 100 training frames.
%   learn_input_side: initial V from the row space of all frames' stacked linear
%   estimates of the first 30 frames, then 6 ALS iterations over all 100 (and the U prior Ru).
%     single   est_general on the test frame alone, rank 16
%     ts100    learned W (and U prior Ru) from 100 frames, U-only LMMSE (est_u_frame)
%     genieW   true input side, U-only LMMSE
%   TDL-A/C/D at 300 ns, window = CP (2.34 us), Np = 8, 6 trials, 20/30 dB.
%   Bar (pre-registered, as L1): Np = 8, 20 dB: gap <= 0.3 for A, C and D (ts100),
%   and below single at 30 dB.

clear; rng(4700);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3; Tcp = 2.34e-6;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
K = 8; Q = 16; r = 16; Nfs = 100;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
f = fc + bins(:)*df; Gam = [element_response(f, 1, P), element_response(f, 2, P)];
snrs = [20 30]; ntrial = 6; Nps = 8; models = {'A', 'C', 'D'}; ds = 300e-9;
names = {'single', 'ts100', 'genieW'};
R = struct('Np', Nps, 'snr_db', snrs, 'names', {names}, 'Nf', Nfs, 'rank', r, 'ntrial', ntrial);
tw = Tcp; L = ceil(M*df*tw) + 10; [F, lamF] = dpss_basis(bins, df, tw, L);
for im = 1:3
    gap = zeros(ntrial, numel(names), numel(Nps), 2);
    for t = 1:ntrial
        [As, Adt] = tdl38901(Q, bins, df, fc, models{im}, ds);
        Htr = cell(1, max(Nfs));
        for k = 1:max(Nfs), [~, A] = tdl38901(Q, bins, df, fc, models{im}, ds); Htr{k} = coupling_matrix(bins, fc, df, seqs, As, A, P, M-1); end
        H = coupling_matrix(bins, fc, df, seqs, As, Adt, P, M-1);
        Pc = mean(sum(abs(H(:, evalb)).^2, 1));
        Wg = {(As .* Gam(:,1).').', (As .* Gam(:,2).').'};
        for is = 1:2
            s2 = Pc / 10^(snrs(is)/10); rp = mmse_rate(H, H, s2, evalb);
            for ip = 1:numel(Nps)
                Np = Nps(ip);
                Ys = cell(1, max(Nfs)); Ss = Ys;
                for k = 1:max(Nfs), Ss{k} = qpsk(M, Np); Ys{k} = Htr{k}*Ss{k} + noise(M, Np, s2); end
                S = qpsk(M, Np); Y = H*S + noise(M, Np, s2);
                Hh = cell(1, 3);
                Hh{1} = est_general(Y, S, T, ones(Q,1), F, lamF, s2, r, 20);
                for k = 1
                    st = learn_input_side(Ys, Ss, T, F, lamF, s2, r, 6, 30);
                    o = cell(1,2); [o{1}, o{2}] = est_u_frame(Y, S, T, F, lamF, s2, st.W, st.Ru); Hh{1+k} = o{1};
                end
                o = est_twoscale('frame', Y, S, T, F, lamF, s2, Wg, mean(abs(F'*Adt.').^2, 'all') / mean(lamF)); Hh{3} = o{1};
                for e = 1:3, gap(t, e, ip, is) = rp - mmse_rate(Hh{e}, H, s2, evalb); end
            end
        end
        fprintf('  TDL-%s trial %d: Np8 @20 %s| @30 %s\n', models{im}, t, sprintf('%.2f ', gap(t,:,1,1)), sprintf('%.2f ', gap(t,:,1,2)));
    end
    key = sprintf('TDL%s_ds300ns', models{im});
    g = squeeze(mean(gap, 1)); ci = squeeze(1.96*std(gap, 0, 1)/sqrt(ntrial));
    R.(key) = struct('L', L, 'gap_bits', g, 'ci95', ci);
    fprintf('%s gap (rows single/ts100/genieW, Np = 8)\n  @20:\n', key); disp(round(g(:,1), 2));
    fprintf('  @30:\n'); disp(round(g(:,2), 2));
end
fid = fopen(fullfile(out, 'week8_twoscale_300ns_c.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
function W = noise(M, n, s2), W = sqrt(s2/2) * (randn(M,n) + 1j*randn(M,n)); end
function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
