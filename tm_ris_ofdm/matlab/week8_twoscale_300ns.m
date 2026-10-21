%WEEK8_TWOSCALE_300NS (L1) Long delay spread with a learned static input side.
%   Per trial: one fixed BS -> RIS hop (As); Nf = 10 training frames with independent
%   RIS -> user hops (other users / earlier positions) and Np pilots each; then a NEW
%   user channel, estimated from its own Np pilots.
%     single    est_general on the test frame alone, rank 16 (as week6_channels_fix)
%     twoscale  input side W learned from the training frames (r = Q = 16), U-only LMMSE
%     genieW    true input side (As o Gamma_s), U-only LMMSE (upper bound for the idea)
%   TDL-A/C/D at 300 ns (and 100 ns for reference), window = CP (2.34 us) at 300 ns,
%   1 us at 100 ns, as before. Np = 2, 4, 8. S1, Q = 16, M = 256, 6 trials, 20/30 dB.
%   Bar (pre-registered): 300 ns, 20 dB, Np = 8: twoscale gap <= 0.3 bit/s/Hz for
%   A, C and D, and below single at 30 dB.

clear; rng(4300);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3; Tcp = 2.34e-6;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
K = 8; Q = 16; Nf = 10; r = 16;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1);
T = {coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state1', M-1), ...
     coupling_matrix(bins, fc, df, base, ones(1,M), ones(1,M), 'state2', M-1)};
f = fc + bins(:)*df; Gam = [element_response(f, 1, P), element_response(f, 2, P)];
snrs = [20 30]; ntrial = 6; Nps = [2 4 8];
models = {'A', 'C', 'D'}; dss = [300e-9 100e-9];
names = {'single', 'twoscale', 'genieW'};
R = struct('Np', Nps, 'snr_db', snrs, 'names', {names}, 'Nf', Nf, 'rank', r, 'ntrial', ntrial);
for id = 1:numel(dss)
    for im = 1:3
        tw = min(Tcp, 10*dss(id)); L = ceil(M*df*tw) + 10;
        [F, lamF] = dpss_basis(bins, df, tw, L);
        gap = zeros(ntrial, numel(names), numel(Nps), 2);
        for t = 1:ntrial
            [As, Adt] = tdl38901(Q, bins, df, fc, models{im}, dss(id));     % test user
            Adtr = cell(1, Nf);
            for k = 1:Nf, [~, Adtr{k}] = tdl38901(Q, bins, df, fc, models{im}, dss(id)); end
            H = coupling_matrix(bins, fc, df, seqs, As, Adt, P, M-1);
            Pc = mean(sum(abs(H(:, evalb)).^2, 1));
            Htr = cellfun(@(A) coupling_matrix(bins, fc, df, seqs, As, A, P, M-1), Adtr, 'UniformOutput', false);
            Wg = {(As .* Gam(:,1).').', (As .* Gam(:,2).').'};
            for is = 1:2
                s2 = Pc / 10^(snrs(is)/10); rp = mmse_rate(H, H, s2, evalb);
                for ip = 1:numel(Nps)
                    Np = Nps(ip);
                    Ys = cell(1, Nf); Ss = Ys;
                    for k = 1:Nf, Ss{k} = qpsk(M, Np); Ys{k} = Htr{k}*Ss{k} + noise(M, Np, s2); end
                    st = est_twoscale('learn', Ys, Ss, T, F, lamF, s2, r, 8);
                    S = qpsk(M, Np); Y = H*S + noise(M, Np, s2);
                    Hs = est_general(Y, S, T, ones(Q,1), F, lamF, s2, r, 20);
                    o = est_twoscale('frame', Y, S, T, F, lamF, s2, st.W, st.sclU); Ht = o{1};
                    og = est_twoscale('frame', Y, S, T, F, lamF, s2, Wg, ...
                                      mean(abs(F'*Adt.').^2, 'all') / mean(lamF)); Hg = og{1};
                    Hh = {Hs, Ht, Hg};
                    for e = 1:3, gap(t, e, ip, is) = rp - mmse_rate(Hh{e}, H, s2, evalb); end
                end
            end
            fprintf('  TDL-%s %d ns trial %d done\n', models{im}, round(dss(id)*1e9), t);
        end
        key = sprintf('TDL%s_ds%dns', models{im}, round(dss(id)*1e9));
        g = squeeze(mean(gap, 1)); ci = squeeze(1.96*std(gap, 0, 1)/sqrt(ntrial));   % names x Np x snr
        R.(key) = struct('window_s', tw, 'L', L, 'gap_bits', g, 'ci95', ci);
        fprintf('%s gap (rows single/twoscale/genieW, cols Np = 2/4/8)\n  @20:\n', key); disp(round(g(:,:,1), 2));
        fprintf('  @30:\n'); disp(round(g(:,:,2), 2));
    end
end
fid = fopen(fullfile(out, 'week8_twoscale_300ns.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function S = qpsk(M, n), S = (sign(randn(M,n)) + 1j*sign(randn(M,n))) / sqrt(2); end
function W = noise(M, n, s2), W = sqrt(s2/2) * (randn(M,n) + 1j*randn(M,n)); end
function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
