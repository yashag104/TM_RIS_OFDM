function st = learn_input_side(Ys, Ss, T, F, lam, s2, r, iters)
%LEARN_INPUT_SIDE Learn the static input side W_s = F V_s from Nf frames (G = 1).
%   Improved version of est_twoscale('learn'): the initial V comes from ALL frames.
%   Each frame's linear first stage (pooled LMMSE, as in est_general) estimates
%       [C_1,f  C_2,f] = U_f [V_1^T  V_2^T]          (L x 2L),
%   so stacking the frames vertically, [C_f]_{f=1..Nf} = [U_1; ...; U_Nf] [V_1^T V_2^T],
%   has the shared input side as its row space: V from the top-r right singular
%   vectors. Then ALS (U_f per frame, V jointly over frames), as before.
%   The prior of U is also learned (empirical Bayes over frames): the learned W is an
%   arbitrary mixture of the true input side, so the U columns are correlated:
%       Cov(vec U) = kron(Ru, diag(lam)),  Ru = mean_f U_f^T diag(1/lam) conj(U_f) / L.
%   Returns W, V, sclU, U and Ru (use with EST_U_FRAME).

Nf = numel(Ys); L = size(F, 2);
pl = kron(lam(:), lam(:)); w = [pl; pl];
Cst = zeros(Nf*L, 2*L);
for f = 1:Nf
    [M, Np] = size(Ss{f}); y = Ys{f}(:); Phi = zeros(M*Np, 2*L*L);
    for s = 1:2
        for p = 1:Np
            Gs = T{s} * (Ss{f}(:,p) .* F);
            for i = 1:L, Phi((p-1)*M+(1:M), (s-1)*L*L+(i-1)*L+(1:L)) = F(:, i) .* Gs; end
        end
    end
    sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
    Dw = sqrt(sc*w);
    c = Dw .* (((Dw .* (Phi'*Phi)) .* Dw.' + s2*eye(numel(w))) \ (Dw .* (Phi'*y)));
    Cst((f-1)*L+(1:L), :) = [reshape(c(1:L*L), L, L).' reshape(c(L*L+1:end), L, L).'];
end
[~, Sg, Wc] = svd(Cst, 'econ');
r = min(r, 2*L);
sv = sqrt(Sg(1:r,1:r) / sqrt(Nf));
V = {conj(Wc(1:L, 1:r)) * sv, conj(Wc(L+1:end, 1:r)) * sv};
wl = repmat(lam(:), r, 1);
sclV = mean(abs([V{1}(:); V{2}(:)]).^2 ./ [wl; wl]); sclU = sclV;
U = cell(1, Nf); Ru = sclU*eye(r);
for it = 1:iters
    W = {F*V{1}, F*V{2}};
    for f = 1:Nf, U{f} = est_u(Ys{f}, Ss{f}, T, F, lam, W, s2, Ru); end
    Ru = learn_ru(U, lam);
    AtA = zeros(2*L*r); Aty = zeros(2*L*r, 1);
    for f = 1:Nf
        [M, Np] = size(Ss{f}); FU = F*U{f}; A = zeros(M*Np, 2*L*r);
        for p = 1:Np
            for s = 1:2
                Gs = T{s} * (Ss{f}(:,p) .* F);
                for b = 1:r, A((p-1)*M+(1:M), (s-1)*L*r+(b-1)*L+(1:L)) = FU(:, b) .* Gs; end
            end
        end
        AtA = AtA + A'*A; Aty = Aty + A'*Ys{f}(:);
    end
    Dv = sqrt(sclV*[wl; wl]);
    v = Dv .* (((Dv .* AtA) .* Dv.' + s2*eye(2*L*r)) \ (Dv .* Aty));
    V = {reshape(v(1:L*r), L, r), reshape(v(L*r+1:end), L, r)};
    sclV = mean(abs(v).^2 ./ [wl; wl]);
end
W = {F*V{1}, F*V{2}};
for f = 1:Nf, U{f} = est_u(Ys{f}, Ss{f}, T, F, lam, W, s2, Ru); end
Ru = learn_ru(U, lam); sclU = real(trace(Ru))/r;
st = struct('V', {V}, 'W', {W}, 'sclU', sclU, 'U', {U}, 'Ru', Ru);
end

function Ru = learn_ru(U, lam)
L = numel(lam); Ru = 0;
for f = 1:numel(U), Ru = Ru + U{f}.' * (conj(U{f}) ./ lam(:)) / L; end
Ru = Ru / numel(U); Ru = (Ru + Ru')/2 + 1e-6*real(trace(Ru))/size(Ru,1)*eye(size(Ru,1));
end

function U = est_u(Y, S, T, F, lam, W, s2, Ru)
[~, U] = est_u_frame(Y, S, T, F, lam, s2, W, Ru);
end
