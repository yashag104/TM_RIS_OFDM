function out = est_twoscale(mode, varargin)
%EST_TWOSCALE Two-timescale timing-only estimator (identical sequences, G = 1).
%   The coupling matrix is H = sum_s T_s o ((F U)(W_s).'), W_s = F V_s (M x r).
%   W_s is the INPUT side: the BS -> RIS hop and the element response. With a fixed
%   BS and RIS it does not change, so it is learned once (slow timescale) from Nf
%   pilot frames that share it, with the full inner dimension (r = Q). Each new frame
%   then needs only U (L x r), a LINEAR LMMSE problem (fast timescale).
%
%   st = est_twoscale('learn', Ys, Ss, T, F, lam, s2, r, iters)
%        Ys, Ss: 1 x Nf cells of M x Np pilots. Returns st.W = {W_1, W_2}, st.V,
%        st.sclU (prior power of U per DPSS weight), st.U (per training frame).
%   [Hh, U] = est_twoscale('frame', Y, S, T, F, lam, s2, W, sclU [, U0, P0w])
%        U-only LMMSE for one frame given W. Optional U0 (prior mean) and P0w
%        (prior variance scale of U - U0 relative to sclU * lam) for tracking.

switch mode
    case 'learn', out = learn(varargin{:});
    case 'frame'
        [Hh, U] = frame(varargin{:}); out = {Hh, U};
    otherwise, error('mode');
end
end

function st = learn(Ys, Ss, T, F, lam, s2, r, iters)
Nf = numel(Ys); L = size(F, 2);
% init: single-frame estimate on frame 1, Lam absorbed into U
[~, s0] = est_general_ws(Ys{1}, Ss{1}, T, F, lam, s2, r, 20, [], false);
r = size(s0.U, 2);
V = s0.V; U0 = s0.U * s0.Lam{1};
wl = repmat(lam(:), r, 1);
sclU = mean(abs(U0(:)).^2 ./ wl); sclV = mean(abs([V{1}(:); V{2}(:)]).^2 ./ [wl; wl]);
U = cell(1, Nf);
for it = 1:iters
    W = {F*V{1}, F*V{2}};
    for f = 1:Nf, U{f} = ustep(Ys{f}, Ss{f}, T, F, W, s2, sclU*wl, zeros(L*r, 1)); end
    sclU = mean(cellfun(@(u) mean(abs(u(:)).^2 ./ wl), U));
    AtA = zeros(2*L*r); Aty = zeros(2*L*r, 1);           % V step, all frames jointly
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
st = struct('V', {V}, 'W', {{F*V{1}, F*V{2}}}, 'sclU', sclU, 'U', {U});
end

function [Hh, U] = frame(Y, S, T, F, lam, s2, W, sclU, U0, P0w)
L = size(F, 2); r = size(W{1}, 2); wl = repmat(lam(:), r, 1);
if nargin < 9 || isempty(U0), u0 = zeros(L*r, 1); pv = sclU*wl;
else, u0 = U0(:); pv = P0w*sclU*wl; end
U = ustep(Y, S, T, F, W, s2, pv, u0);
FU = F*U; Hh = T{1} .* (FU*W{1}.') + T{2} .* (FU*W{2}.');
end

function U = ustep(Y, S, T, F, W, s2, pv, u0)
% LMMSE of U given W: y = A u + w, u ~ CN(u0, diag(pv))
[M, Np] = size(S); L = size(F, 2); r = size(W{1}, 2);
A = zeros(M*Np, L*r);
for p = 1:Np
    g = T{1} * (S(:,p) .* W{1}) + T{2} * (S(:,p) .* W{2});     % M x r
    for a = 1:r, A((p-1)*M+(1:M), (a-1)*L+(1:L)) = F .* g(:, a); end
end
y = Y(:) - A*u0; D = sqrt(pv);
u = u0 + D .* (((D .* (A'*A)) .* D.' + s2*eye(L*r)) \ (D .* (A'*y)));
U = reshape(u, L, r);
end
