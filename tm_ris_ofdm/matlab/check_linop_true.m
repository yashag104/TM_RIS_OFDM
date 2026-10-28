% Operator of lin_op built from the TRUE coefficients C = F^H A conj(F): if the construction is right, its error on
% sampled columns is at the delay-basis floor (~ -96 dB). Also the DPSS prior range. Prints only.
P = params_openris(); fc = 3.594e9; df = 30e3; M = 3276; bins = -M/2:M/2-1; Q = 16; f = fc + bins(:)*df;
base = [ones(1, 6), 2*ones(1, 2)]; seqs = repmat(base, Q, 1); evalb = round(M/4)+1:round(3*M/4);
L = ceil(M*df*0.2e-6) + 10; [F, lamF] = dpss_basis(bins, df, 0.2e-6, L);
rng(13501); [As, Ad] = phys_channel(Q, bins, df, fc, struct('P', 5, 'K_db', 10, 'rms', 30e-9));
op = fast_coupling(bins, fc, df, seqs, As, Ad, P);
cols = evalb(round(linspace(1, numel(evalb), 48))); E = zeros(M, 48); E(sub2ind(size(E), cols, 1:48)) = 1; HE = op.H(E);
gt = {element_response(f, 1, P), element_response(f, 2, P)};
C = F' * (Ad.' * As) * conj(F);
[Uc, Sc, Vc] = svd(C); sq = sqrt(Sc); FU = F*Uc*sq; FV = F*conj(Vc)*sq; seqv = repmat(base, L, 1);
o1 = fast_coupling(bins, fc, df, seqv, (gt{1} .* FV).', FU.', 'state1');
o2 = fast_coupling(bins, fc, df, seqv, (gt{2} .* FV).', FU.', 'state2');
fprintf('operator from true C: %.1f dB\n', 10*log10(norm(o1.H(E) + o2.H(E) - HE, 'fro')^2/norm(HE, 'fro')^2));
% same with the virtual elements given directly as the true hops (no basis): checks fast_coupling's virtual-element use
o1 = fast_coupling(bins, fc, df, seqs, gt{1}.' .* As, Ad, 'state1'); o2 = fast_coupling(bins, fc, df, seqs, gt{2}.' .* As, Ad, 'state2');
fprintf('true hops through indicator kernels: %.1f dB\n', 10*log10(norm(o1.H(E) + o2.H(E) - HE, 'fro')^2/norm(HE, 'fro')^2));
fprintf('DPSS prior lam: max %.2e, min %.2e; energy of C per DPSS index (dB): %s\n', max(lamF), min(lamF), sprintf('%.0f ', 10*log10(sum(abs(C).^2, 2)/sum(abs(C(:)).^2))));
