%WEEK8_TIMING_GUARD (L2) Absorb timing offsets with a guarded delay window, no search.
%   A switching-clock offset of d slots multiplies the timing kernel as
%   T_d = D T D^*, D = diag(exp(-j 2 pi m d / K)), so
%       T_d o ((F U)(F V).') = T o ((D F U)(D^* F V).'),
%   i.e. it is a delay of d*Tslot on the output hop and an advance on the input hop.
%   A delay basis covering [-g, tau + g] on both hops therefore represents ANY offset
%   |d| Tslot <= g with the NOMINAL timing kernels. The same holds for a fixed skew
%   per element (each element's two hops shift by its own d_q), since the factors
%   sum over elements.
%   Receivers (all with nominal timing kernels unless stated):
%     plain    window [0, 0.3 us]                     (= ignoring the offset, as before)
%     guard    window [-g, 0.3 us + g], g = guard     (new, no timing search)
%     known    true offset in the kernels, window [0, 0.3 us]   (reference)
%   Cases: (a) common offset, uniform in +-0.02 slot (83 ns), guard 100 ns;
%          (b) common offset, uniform in +-0.05 slot (208 ns), guard 250 ns;
%          (c) fixed per-element skew, uniform in +-0.01 slot (42 ns), guard 50 ns
%              (known = perfect channel here; reported as gap).
%   Bar (pre-registered): case (a) at 30 dB, guard within 1 dB of known.
%   S1, K = 8, Q = 16, M = 256, Np = 2, physical channel, 12 trials, 20/30 dB.

clear; rng(4200);
here = fileparts(mfilename('fullpath'));
out = fullfile(here, '..', 'results');
P = params_openris();
fc = 3.594e9; df = 30e3; Tslot = 1/df/8;
M = 256; bins = -M/2:M/2-1; evalb = 65:192;
Q = 16; K = 8; Np = 2;
base = [ones(1, 6), 2*ones(1, 2)]; st = repmat(base, Q, 1);
popt = struct('P', 5, 'K_db', 10, 'rms', 30e-9);
tw = 0.3e-6; [F0, l0] = dpss_basis(bins, df, tw, 16);
E0 = (0:K)/K;
Tnom = {coupling_edges(bins, fc, df, E0, base, ones(1,M), ones(1,M), 'state1', M-1), ...
        coupling_edges(bins, fc, df, E0, base, ones(1,M), ones(1,M), 'state2', M-1)};
snrs = [20 30]; slope = [0.31 0.33]; ntrial = 12;
cases = struct('name', {'a_common_83ns', 'b_common_208ns', 'c_skew_42ns'}, ...
               'dmax', {0.02, 0.05, 0.01}, 'guard', {100e-9, 250e-9, 50e-9}, 'skew', {false, false, true});
R = struct('snr_db', snrs, 'receivers', {{'plain', 'guard', 'known'}}, 'ntrial', ntrial);
for ic = 1:numel(cases)
    c = cases(ic); g = c.guard;
    Lg = 16 + ceil(M*df*2*g) + 1;
    [Fg, lg] = dpss_basis(bins, df, tw + 2*g, Lg);
    Fg = exp(1j*2*pi*bins(:)*df*g) .* Fg;                % delays [0, tw+2g] -> [-g, tw+g]
    loss = zeros(ntrial, 2, 3);
    for t = 1:ntrial
        if c.skew, dq = (2*rand(Q, 1) - 1) * c.dmax; else, dq = repmat((2*rand - 1) * c.dmax, Q, 1); end
        [As, Ad] = phys_channel(Q, bins, df, fc, popt);
        H = coupling_edges(bins, fc, df, ((0:K) + dq)/K, st, As, Ad, P, M-1);
        Pc = mean(sum(abs(H(:, evalb)).^2, 1));
        for is = 1:2
            s2 = Pc / 10^(snrs(is)/10);
            S = (sign(randn(M,Np)) + 1j*sign(randn(M,Np))) / sqrt(2);
            Y = H*S + sqrt(s2/2)*(randn(M,Np) + 1j*randn(M,Np));
            rp = mmse_rate(H, H, s2, evalb);
            Hh = {est_general(Y, S, Tnom, ones(Q,1), F0, l0, s2, 8, 20), ...
                  est_general(Y, S, Tnom, ones(Q,1), Fg, lg, s2, 8, 20), []};
            if c.skew, Hh{3} = H;
            else, Hh{3} = est_general(Y, S, shiftk(Tnom, bins, dq(1), K), ones(Q,1), F0, l0, s2, 8, 20); end
            for e = 1:3, loss(t, is, e) = (rp - mmse_rate(Hh{e}, H, s2, evalb)) / slope(is); end
        end
    end
    m = squeeze(mean(loss, 1)); ci = squeeze(1.96*std(loss, 0, 1)/sqrt(ntrial));
    R.(c.name) = struct('dmax_slot', c.dmax, 'dmax_ns', c.dmax*Tslot*1e9, 'guard_ns', g*1e9, 'L_guard', Lg, ...
                        'loss_db', m, 'ci95_db', ci);
    fprintf('%s (|d| <= %.0f ns, guard %.0f ns, L %d -> %d): loss plain/guard/known @20 %.2f/%.2f/%.2f | @30 %.2f/%.2f/%.2f dB\n', ...
        c.name, c.dmax*Tslot*1e9, g*1e9, 16, Lg, m(1,:), m(2,:));
end
fid = fopen(fullfile(out, 'week8_timing_guard.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);

function T = shiftk(Tnom, bins, d, K)
D = exp(-1j*2*pi*bins(:)*d/K);
T = {D .* Tnom{1} .* D', D .* Tnom{2} .* D'};
end
function r = mmse_rate(Hh, H, s2, evalb)
M = size(H, 1);
W = (Hh'*Hh + s2*eye(M)) \ Hh'; G = W*H; sg = abs(diag(G)).^2;
sinr = sg ./ (sum(abs(G).^2, 2) - sg + s2*sum(abs(W).^2, 2));
r = mean(log2(1 + sinr(evalb)));
end
