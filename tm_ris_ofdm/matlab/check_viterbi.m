clear; rng(1);
t = make_trellis(); n = 506;
u = randi([0 1], n, 1); c = conv_enc([u; zeros(6,1)], t);
llr = 10*(1 - 2*c);                                   % noiseless, strong
d = viterbi_soft(llr, t, n + 6); fprintf('noiseless errors: %d\n', sum(d(1:n) ~= u));
x = 1 - 2*c; sig = 0.8; y = x + sig*randn(size(x)); llr = 2*y/sig^2;
d = viterbi_soft(llr, t, n + 6); fprintf('noisy (Eb/N0 %.1f dB): bit errors %d of %d, uncoded hard errors %d\n', ...
   10*log10(1/(2*0.5*sig^2)), sum(d(1:n) ~= u), n, sum((y(1:2:end) < 0) ~= c(1:2:end)));
function t = make_trellis()
g = [bin2vec(oct2dec_('133')); bin2vec(oct2dec_('171'))];   % 2 x 7, MSB = current bit
ns = 64; t.next = zeros(ns, 2); t.out = zeros(ns, 2, 2);
for s = 0:ns-1
    for u = 0:1
        reg = [u, bitget(s, 6:-1:1)];                         % current input + 6 memory bits
        t.out(s+1, u+1, :) = mod(g * reg.', 2);
        t.next(s+1, u+1) = bitshift(s, -1) + u*32;
    end
end
end
function v = bin2vec(d), v = bitget(d, 7:-1:1); end
function d = oct2dec_(s), d = 0; for ch = s, d = d*8 + (ch - '0'); end, end

function c = conv_enc(u, t)
s = 0; c = zeros(2*numel(u), 1);
for i = 1:numel(u)
    c(2*i-1:2*i) = squeeze(t.out(s+1, u(i)+1, :));
    s = t.next(s+1, u(i)+1);
end
end

function u = viterbi_soft(llr, t, n)
% Vectorised soft Viterbi: maximise sum (1 - 2c) * llr; start and end in state 0.
% Next state ns = (s >> 1) + 32 u, so ns has predecessors 2*mod(ns,32) and +1, input floor(ns/32).
ns = 64; nsv = (0:ns-1).';
pA = 2*mod(nsv, 32) + 1; pB = pA + 1; uin = floor(nsv/32) + 1;     % 1-based
oA = [t.out(sub2ind([ns 2 2], pA, uin, ones(ns,1))), t.out(sub2ind([ns 2 2], pA, uin, 2*ones(ns,1)))];
oB = [t.out(sub2ind([ns 2 2], pB, uin, ones(ns,1))), t.out(sub2ind([ns 2 2], pB, uin, 2*ones(ns,1)))];
sA = 1 - 2*oA; sB = 1 - 2*oB;
pm = -inf(ns, 1); pm(1) = 0; choice = false(ns, n);
for i = 1:n
    l = llr(2*i-1:2*i);
    mA = pm(pA) + sA * l; mB = pm(pB) + sB * l;
    choice(:, i) = mB > mA; pm = max(mA, mB);
end
u = zeros(n, 1); st = 1;
for i = n:-1:1
    u(i) = uin(st) - 1;
    if choice(st, i), st = pB(st); else, st = pA(st); end
end
end

