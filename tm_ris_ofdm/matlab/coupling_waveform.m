function H = coupling_waveform(m, fc, df, W, As, Ad, Plist, hmax)
%COUPLING_WAVEFORM Coupling for arbitrary state-weight waveforms and per-element models.
%   W: Q x Nt x 2, the weight of state 1 and state 2 at Nt uniform instants over one
%   period (w1 + w2 = 1; a ramp between states models a finite switching time, with
%   the reflection taken as the weighted mix of the two state responses).
%   Plist: element model per element (cell of structs, Q x 1) or a single struct.
%   b_q^[h](f) = sum_s c_{q,s}[h] Gamma_{q,s}(f),  c_{q,s}[h] = mean_t w_s(t) e^{-j2pi h t}.

m = m(:).'; Mn = numel(m); Q = size(W, 1); Nt = size(W, 2);
fin = fc + m*df;
hs = -hmax:hmax;
H = zeros(Mn, Mn);
for q = 1:Q
    if iscell(Plist), Pq = Plist{q}; else, Pq = Plist; end
    Gam = [element_response(fin, 1, Pq); element_response(fin, 2, Pq)];
    c = zeros(2, numel(hs));
    for s = 1:2
        X = fft(squeeze(W(q, :, s))) / Nt;                 % X(h+1) = mean w e^{-j2pi h n/Nt}
        c(s, :) = X(mod(hs, Nt) + 1);
    end
    B = c.' * Gam;
    for ih = 1:numel(hs)
        h = hs(ih);
        src = find(m + h >= m(1) & m + h <= m(end));
        dst = src + h;
        idx = sub2ind([Mn Mn], dst, src);
        H(idx) = H(idx) + Ad(q, dst) .* B(ih, src) .* As(q, src);
    end
end
end
