% Full band: error of the best flat (two-scalar) element model given the TRUE cascade, i.e. the model error that
% dispersion alone imposes on a calibrated flat baseline. Sampled columns; prints only.
P = params_openris(); fc = 3.594e9; df = 30e3; M = 3276; bins = -M/2:M/2-1; Q = 16; f = fc + bins(:)*df;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1); evalb = round(M/4)+1:round(3*M/4);
cols = evalb(round(linspace(1, numel(evalb), 48))); E = zeros(M, 48); E(sub2ind(size(E), cols, 1:48)) = 1;
g1c = element_response(fc, 1, P); g2c = element_response(fc, 2, P); G1 = element_response(f, 1, P);
for t = 1:3
    rng(13500 + t); [As, Ad] = phys_channel(Q, bins, df, fc, struct('P', 5, 'K_db', 10, 'rms', 30e-9));
    HE = fast_coupling(bins, fc, df, seqs, As, Ad, P).H(E);
    Asf = (G1.' / g1c) .* As;                                  % flat model absorbs Gamma_1(f) into the hop
    o1 = fast_coupling(bins, fc, df, seqs, g1c * Asf, Ad, 'state1'); o2 = fast_coupling(bins, fc, df, seqs, g2c * Asf, Ad, 'state2');
    fprintf('channel %d: flat-model error given the true cascade %.1f dB\n', t, 10*log10(norm(o1.H(E) + o2.H(E) - HE, 'fro')^2/norm(HE, 'fro')^2));
end
