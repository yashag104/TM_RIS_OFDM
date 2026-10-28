% smoke test for week13_review_baselines: static-channel orientation and timing (writes results/smoke_week13.json)
P = params_openris(); fc = 3.594e9; df = 30e3; M = 256; bins = -M/2:M/2-1; Q = 16; fin = fc + bins*df;
rng(1); [As, Ad] = phys_channel(Q, bins, df, fc, struct('P', 5, 'K_db', 10, 'rms', 30e-9));
for s = 1:2
    Hs = coupling_matrix(bins, fc, df, repmat(s*ones(1, 8), Q, 1), As, Ad, P, M-1);
    h = sum(Ad .* As, 1) .* element_response(fin, s, P);
    fprintf('state %d: static diag rel err %.2e, offdiag energy %.2e\n', s, norm(diag(Hs).' - h)/norm(h), norm(Hs - diag(diag(Hs)), 'fro')/norm(Hs, 'fro'));
end
% (full one-channel run done 2026-10-28: 814 s)
