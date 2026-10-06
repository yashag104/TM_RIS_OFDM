function H = coupling_edges(m, fc, df, edges, states, As, Ad, P, hmax)
%COUPLING_EDGES Coupling matrix for switching patterns with ARBITRARY edge times.
%   Like COUPLING_MATRIX, but element q holds state states(q,k) on the interval
%   [edges(q,k), edges(q,k+1)) of the normalised period [0, 1). Edges may be
%   offset or jittered away from the ideal k/K grid. The harmonic coefficient is
%   the exact Fourier integral of the piecewise-constant pattern:
%     b_q^[h](f) = sum_s c_{q,s}[h] Gamma_s(f),
%     c_{q,s}[h] = sum_{k: state s} int_{e_k}^{e_{k+1}} exp(-j 2 pi h t) dt.
%   edges: Q x (K+1), increasing, spanning one period (wrap handled by mod 1).
%   P: element model (struct) or 'state1'/'state2' indicator, as in ELEMENT_RESPONSE.

m = m(:).'; Mn = numel(m); Q = size(states, 1);
fin = fc + m*df;
Gam = [element_response(fin, 1, P); element_response(fin, 2, P)];   % 2 x M
hs = -hmax:hmax;
H = zeros(Mn, Mn);
for q = 1:Q
    c = zeros(2, numel(hs));
    for k = 1:size(states, 2)
        a = edges(q, k); b = edges(q, k+1);
        c(states(q, k), :) = c(states(q, k), :) + seg(hs, a, b);
    end
    B = c.' * Gam;                                   % (2hmax+1) x M: b^[h](f_m)
    for ih = 1:numel(hs)
        h = hs(ih);
        src = find(m + h >= m(1) & m + h <= m(end));
        dst = src + h;
        idx = sub2ind([Mn Mn], dst, src);
        H(idx) = H(idx) + Ad(q, dst) .* B(ih, src) .* As(q, src);
    end
end
end

function v = seg(hs, a, b)
% int_a^b exp(-j 2 pi h t) dt for each h (h = 0 gives b - a).
v = zeros(1, numel(hs));
nz = hs ~= 0;
v(~nz) = b - a;
v(nz) = (exp(-1j*2*pi*hs(nz)*a) - exp(-1j*2*pi*hs(nz)*b)) ./ (1j*2*pi*hs(nz));
end
