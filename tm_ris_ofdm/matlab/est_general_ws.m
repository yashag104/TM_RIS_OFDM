function [Hh, state] = est_general_ws(Y, S, T, F, lam, s2, r, iters, init, fixV)
%EST_GENERAL_WS EST_GENERAL with warm start and an optional frozen input-side factor.
%   init  []  -> cold start (pooled linear init, as EST_GENERAL)
%         struct with fields U, V, Lam -> warm start from a previous estimate
%   fixV  true -> keep V (input side: the static BS -> RIS hop and element response)
%         and update only Lam and U, the cheap part that the moving user changes
%   state returns U, V, Lam for the next call.
%
%EST_GENERAL Timing-only estimator for ARBITRARY per-element switching patterns.
%   Elements are grouped by their (distinct) control sequence g = 1..G. For each
%   group the receiver knows the per-state Toeplitz kernels T{g,s} (timing only,
%   built with the 'state1'/'state2' indicator model; no element model).
%   Physical channels share propagation paths, so with shared factors
%       H = sum_g sum_s T{g,s} o ( (F U) Lam_g (F V_s).' )
%   U (L x r) output side, V_s (L x r) input side per state (absorbs the element
%   dispersion), Lam_g (r x r) per distinct sequence. Identical sequences (G = 1)
%   and cyclic shifts (T{g,s} = D_k T_s D_k^*) are special cases.
%   Block ALS: the Lam step is a convex linear solve; init from the pooled model.
%   Y, S: M x Np received / known pilot symbols. grp is unused except for G.
%   Returns the estimated M x M coupling matrix.

[M, Np] = size(S); L = size(F, 2); G = size(T, 1); y = Y(:); %#ok<NASGU>
r = min(r, L);                                     % rank cannot exceed the delay-basis size
Gk = cell(G, 2, Np);                               % T{g,s} (S_p .* F), M x L
for g = 1:G, for s = 1:2, for p = 1:Np
    Gk{g,s,p} = T{g,s} * (S(:,p) .* F);
end, end, end
if isempty(init)
% pooled init: H ~ sum_s (sum_g T{g,s}) o (F C_s F.')
pl = kron(lam(:), lam(:)); Phi = zeros(M*Np, 2*L*L);
for s = 1:2
    for p = 1:Np
        Gs = zeros(M, L); for g = 1:G, Gs = Gs + Gk{g,s,p}; end
        for i = 1:L, Phi((p-1)*M+(1:M), (s-1)*L*L+(i-1)*L+(1:L)) = F(:, i) .* Gs; end
    end
end
w = [pl; pl];
sc = max(real(y'*y)/numel(y) - s2, eps) / max(real(sum(abs(Phi).^2, 1) * w)/numel(y), eps);
Dw = sqrt(sc*w);
c = Dw .* (((Dw .* (Phi'*Phi)) .* Dw.' + s2*eye(numel(w))) \ (Dw .* (Phi'*y)));
[Uc, Sg, Wc] = svd([reshape(c(1:L*L), L, L).' reshape(c(L*L+1:end), L, L).']);
U = Uc(:, 1:r) * sqrt(Sg(1:r,1:r));
V = {conj(Wc(1:L, 1:r)) * sqrt(Sg(1:r,1:r)), conj(Wc(L+1:end, 1:r)) * sqrt(Sg(1:r,1:r))};
Lam = repmat({eye(r)}, 1, G);
else
    U = init.U; V = init.V; Lam = init.Lam; r = size(U, 2);
end
rr = 1e-3;
% DPSS prior for the hop factors (as in est_lmmse): scale from the init
wl = repmat(lam(:), r, 1);
scl = mean(abs([U(:); V{1}(:); V{2}(:)]).^2 ./ [wl; wl; wl]);
Du = sqrt(scl*wl); Dv = sqrt(scl*[wl; wl]);
H0 = [];
if isempty(init)                                   % pooled linear fallback, cold start only
    H0 = zeros(M);
    C1 = reshape(c(1:L*L), L, L).'; C2 = reshape(c(L*L+1:end), L, L).';
    for g = 1:G, H0 = H0 + T{g,1} .* (F*C1*F.') + T{g,2} .* (F*C2*F.'); end
end
for it = 1:iters
    FU = F*U;
    A = zeros(M*Np, G*r*r);                        % Lam step (convex, all groups)
    for g = 1:G
        for p = 1:Np
            GV = Gk{g,1,p}*V{1} + Gk{g,2,p}*V{2};
            for a = 1:r, A((p-1)*M+(1:M), (g-1)*r*r+(a-1)*r+(1:r)) = FU(:, a) .* GV; end
        end
    end
    lv = rsolve(A, y, rr);
    for g = 1:G, Lam{g} = reshape(lv((g-1)*r*r+(1:r*r)), r, r).'; end
    A = zeros(M*Np, L*r);                          % U step
    for p = 1:Np
        for g = 1:G
            GV = (Gk{g,1,p}*V{1} + Gk{g,2,p}*V{2}) * Lam{g}.';
            for a = 1:r
                A((p-1)*M+(1:M), (a-1)*L+(1:L)) = A((p-1)*M+(1:M), (a-1)*L+(1:L)) + F .* GV(:, a);
            end
        end
    end
    U = reshape(Du .* (((Du .* (A'*A)) .* Du.' + s2*eye(L*r)) \ (Du .* (A'*y))), L, r); FU = F*U;
    if fixV, continue; end
    A = zeros(M*Np, 2*L*r);                        % V step (both states)
    for p = 1:Np
        for g = 1:G
            FUL = FU * Lam{g};
            for s = 1:2
                for b = 1:r
                    idx = (s-1)*L*r + (b-1)*L + (1:L);
                    A((p-1)*M+(1:M), idx) = A((p-1)*M+(1:M), idx) + FUL(:, b) .* Gk{g,s,p};
                end
            end
        end
    end
    v = Dv .* (((Dv .* (A'*A)) .* Dv.' + s2*eye(2*L*r)) \ (Dv .* (A'*y)));
    V = {reshape(v(1:L*r), L, r), reshape(v(L*r+1:end), L, r)};
end
Hh = zeros(M); FU = F*U;
for g = 1:G
    for s = 1:2
        Hh = Hh + T{g,s} .* ((FU*Lam{g}) * (F*V{s}).');
    end
end
% choose the better fit to the pilots (no truth used)
if ~isempty(H0) && norm(Y - H0*S, 'fro') < norm(Y - Hh*S, 'fro'), Hh = H0; end
state = struct('U', U, 'V', {V}, 'Lam', {Lam});
end

function x = rsolve(A, y, rel)
Gm = A'*A;
x = (Gm + rel*real(trace(Gm))/size(Gm,1)*eye(size(Gm))) \ (A'*y);
end
