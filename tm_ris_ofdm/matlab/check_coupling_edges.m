% coupling_edges with ideal edges must equal coupling_matrix exactly.
P = params_openris(); fc = 3.594e9; df = 30e3; M = 64; bins = -M/2:M/2-1; Q = 3; K = 8;
st = randi(2, Q, K); As = randn(Q,M)+1j*randn(Q,M); Ad = randn(Q,M)+1j*randn(Q,M);
E = repmat((0:K)/K, Q, 1);
H1 = coupling_matrix(bins, fc, df, st, As, Ad, P, M-1);
H2 = coupling_edges(bins, fc, df, E, st, As, Ad, P, M-1);
fprintf('ideal edges: rel. diff %.2e\n', norm(H1-H2,'fro')/norm(H1,'fro'));
