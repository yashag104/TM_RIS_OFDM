function [Hh, U] = est_u_frame(Y, S, T, F, lam, s2, W, Ru, U0, P0w)
%EST_U_FRAME U-only LMMSE for one frame given the learned input side W (G = 1).
%   H = sum_s T_s o ((F U) W_s^T), prior vec(U) ~ CN(U0, P0w * kron(Ru, diag(lam))).
%   Ru from LEARN_INPUT_SIDE. U0 = [] and P0w = 1 for a fresh frame; for tracking,
%   U0 = previous U and P0w = relative variance of the change.
[M, Np] = size(S); L = size(F, 2); r = size(W{1}, 2);
if nargin < 9 || isempty(U0), u0 = zeros(L*r, 1); P0w = 1; else, u0 = U0(:); end
A = zeros(M*Np, L*r);
for p = 1:Np
    g = T{1} * (S(:,p) .* W{1}) + T{2} * (S(:,p) .* W{2});
    for a = 1:r, A((p-1)*M+(1:M), (a-1)*L+(1:L)) = F .* g(:, a); end
end
Lc = sqrt(P0w) * kron(chol(Ru, 'lower'), diag(sqrt(lam(:))));      % Cov = Lc Lc'
B = A * Lc;
u = u0 + Lc * ((B'*B + s2*eye(L*r)) \ (B' * (Y(:) - A*u0)));
U = reshape(u, L, r);
FU = F*U; Hh = T{1} .* (FU*W{1}.') + T{2} .* (FU*W{2}.');
end
