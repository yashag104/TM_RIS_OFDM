function [Hh, info] = est_bem_cg(Y, S, F2, lam, s2, hs, tol, maxit)
%EST_BEM_CG CE-BEM LMMSE estimate computed matrix-free (for large harmonic sets).
%   Same estimator as the dense est_bem of week7_strong_baselines: H = sum_{h in hs} S_h diag(F2 c_h),
%   prior c_h ~ CN(0, sc diag(lam)) per harmonic, sc by moment matching. The LMMSE solution
%       c = D (D A^H A D + s2 I)^{-1} D A^H y,   D = sqrt(sc diag(w)),
%   is obtained by conjugate gradients on z with c = D z, so A (M Np x Nh L) is never formed:
%   A c = sum_h shift_h(S_p o (F2 c_h)) per pilot p.
%   Returns the dense M x M estimate (for scoring) and info.iters / info.relres.
if nargin < 7 || isempty(tol), tol = 1e-6; end
if nargin < 8 || isempty(maxit), maxit = 500; end
[M, Np] = size(S); L = size(F2, 2); Nh = numel(hs); y = Y(:);
Afun = @(c) applyA(c, S, F2, hs, M, Np, L);
Ahfun = @(r) applyAh(r, S, F2, hs, M, Np, L);
w = repmat(lam(:), Nh, 1);
% moment matching for the prior scale: column energies of A are sum_p |S_p|^2 |F2(:,l)|^2 over valid rows
colE = zeros(L, Nh);
for k = 1:Nh
    h = hs(k); rows = max(1, 1-h):min(M, M-h);                 % input indices that stay in band
    for p = 1:Np, colE(:, k) = colE(:, k) + sum(abs(S(rows, p) .* F2(rows, :)).^2, 1).'; end
end
sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(colE(:).' * w)/numel(y), eps);
D = sqrt(sc*w);
b = D .* Ahfun(y);
op = @(z) D .* Ahfun(Afun(D .* z)) + s2*z;
z = zeros(size(b)); r = b; p = r; rs = real(r'*r); nb = sqrt(real(b'*b));
it = 0;
for it = 1:maxit
    Ap = op(p); a = rs / real(p'*Ap);
    z = z + a*p; r = r - a*Ap; rn = real(r'*r);
    if sqrt(rn)/nb < tol, break; end
    p = r + (rn/rs)*p; rs = rn;
end
c = reshape(D .* z, L, Nh);
g = F2 * c;                                                 % M x Nh per-harmonic gains
Hh = zeros(M);
for k = 1:Nh
    h = hs(k); cols = max(1, 1-h):min(M, M-h);           % input index m, output m + h
    Hh(sub2ind([M M], cols + h, cols)) = Hh(sub2ind([M M], cols + h, cols)) + g(cols, k).';
end
info = struct('iters', it, 'relres', sqrt(rn)/nb);
end

function out = applyA(c, S, F2, hs, M, Np, L)
C = reshape(c, L, []); G = F2 * C; out = zeros(M, Np);
for k = 1:numel(hs), out = out + shr(G(:, k) .* S, hs(k)); end
out = out(:);
end
function c = applyAh(r, S, F2, hs, M, Np, L)
R = reshape(r, M, Np); A = zeros(M, numel(hs)); Sc = conj(S);
for k = 1:numel(hs), A(:, k) = sum(Sc .* shr(R, -hs(k)), 2); end
c = reshape(F2' * A, [], 1);
end
function Z = shr(X, h)
Z = zeros(size(X)); M = size(X, 1);
if h >= 0, Z(1+h:M, :) = X(1:M-h, :); else, Z(1:M+h, :) = X(1-h:M, :); end
end
