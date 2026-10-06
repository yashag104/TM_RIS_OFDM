P = params_openris(); fc = 3.594e9; df = 30e3; M = 256; bins = -M/2:M/2-1; Q = 2; K = 8;
base = [1 1 1 1 1 1 2 2]; As = randn(Q,M)+1j*randn(Q,M); Ad = randn(Q,M)+1j*randn(Q,M);
Nt = 32768; t = ((0:Nt-1)+0.5)/Nt*K; w1 = double(base(floor(t)+1)==1);
W = zeros(Q,Nt,2); W(:,:,1) = repmat(w1,Q,1); W(:,:,2) = 1-W(:,:,1);
H1 = coupling_matrix(bins, fc, df, repmat(base,Q,1), As, Ad, P, M-1);
H2 = coupling_waveform(bins, fc, df, W, As, Ad, P, M-1);
fprintf('instant switching: rel. diff %.2e\n', norm(H1-H2,'fro')/norm(H1,'fro'));
