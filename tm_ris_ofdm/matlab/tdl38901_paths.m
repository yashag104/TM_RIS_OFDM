function [As, Adp, nu] = tdl38901_paths(Q, bins, df, fc, model, ds, fD)
%TDL38901_PATHS TDL38901 with the RIS -> user hop returned tap by tap, for mobility.
%   As: static BS -> RIS hop (Q x M), exactly as TDL38901.
%   Adp: RIS -> user hop per tap (Q x M x Ltap), sum(Adp, 3) is the TDL38901 hop.
%   nu: Doppler per tap, fD cos(theta) with a random direction per tap (each TDL tap
%   is given one arrival direction, as in TDL38901), so tap l at symbol n is
%   Adp(:,:,l) exp(j 2 pi nu(l) n Tsym).

[tau_n, p_db, los] = table(model);
lam = 3e8 / fc; k = 2*pi/lam; n = round(sqrt(Q));
[gx, gz] = meshgrid(0:n-1, 0:n-1);
r = [gx(:) - (n-1)/2, zeros(Q,1), gz(:) - (n-1)/2] * lam/2;
As = sum(hop(), 3); Adp = hop();
nu = fD * cos(2*pi*rand(1, size(Adp, 3)));

    function A = hop()
        tau = tau_n * ds; pw = 10.^(p_db/10); pw = pw / sum(pw);
        L = numel(tau);
        g = sqrt(pw) .* (randn(1, L) + 1j*randn(1, L)) / sqrt(2);
        if los, g(1) = sqrt(pw(1)) * exp(1j*2*pi*rand); end
        u = randn(L, 3); u(:,2) = abs(u(:,2)); u = u ./ vecnorm(u, 2, 2);
        ph = exp(-1j*k*(r*u.'));
        A = zeros(Q, numel(bins), L);
        for l = 1:L, A(:, :, l) = (ph(:, l) * g(l)) * exp(-1j*2*pi*tau(l)*(bins(:).'*df)); end
    end
end

function [tau, p, los] = table(model)
switch upper(model)
    case 'A'
        tau = [0 0.3819 0.4025 0.5868 0.4610 0.5375 0.6708 0.5750 0.7618 1.5375 1.8978 ...
               2.2242 2.1718 2.4942 2.5119 3.0582 4.0810 4.4579 4.5695 4.7966 5.0066 5.3043 9.6586];
        p = [-13.4 0 -2.2 -4 -6 -8.2 -9.9 -10.5 -7.5 -15.9 -6.6 -16.7 -12.4 -15.2 -10.8 ...
             -11.3 -12.7 -16.2 -18.3 -18.9 -16.6 -19.9 -29.7];
        los = false;
    case 'C'
        tau = [0 0.2099 0.2219 0.2329 0.2176 0.6366 0.6448 0.6560 0.6584 0.7935 0.8213 0.9336 ...
               1.2285 1.3083 2.1704 2.7105 4.2589 4.6003 5.4902 5.6077 6.3065 6.6374 7.0427 8.6523];
        p = [-4.4 -1.2 -3.5 -5.2 -2.5 0 -2.2 -3.9 -7.4 -7.1 -10.7 -11.1 -5.1 -6.8 -8.7 -13.2 ...
             -13.9 -13.9 -15.8 -17.1 -16 -15.7 -21.6 -22.8];
        los = false;
    case 'D'
        tau = [0 0 0.035 0.612 1.363 1.405 1.804 2.596 1.775 4.042 7.937 9.424 9.708 12.525];
        p = [-0.2 -13.5 -18.8 -21 -22.8 -17.9 -20.1 -21.9 -22.9 -27.8 -23.6 -24.8 -30 -27.7];
        los = true;
    otherwise
        error('model must be A, C or D');
end
end
