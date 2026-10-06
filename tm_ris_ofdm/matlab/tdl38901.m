function [As, Ad] = tdl38901(Q, bins, df, fc, model, ds, nf_dist)
%TDL38901 3GPP TR 38.901 TDL power-delay profiles for both RIS hops.
%   model: 'A', 'C' (NLOS) or 'D' (LOS); ds: delay spread in seconds (normalised
%   delays of Table 7.7.2 are scaled by ds). Tap tables from the 38.901 TDL models
%   (as distributed with NVIDIA Sionna v0.19.0, models/TDL-*.json).
%   TDL has no angles, so each tap gets a random arrival direction that sets its
%   geometric phase across the sqrt(Q) x sqrt(Q) lambda/2 element grid; the
%   first TDL-D tap is the deterministic LoS component.
%   nf_dist (optional, metres): if > 0, the LoS path of the RIS -> user hop uses
%   exact spherical distances from each element to a user at that range
%   (near field); otherwise plane waves. Both hops have unit average power.

if nargin < 7, nf_dist = 0; end
[tau_n, p_db, los] = table(model);
lam = 3e8 / fc; k = 2*pi/lam; n = round(sqrt(Q));
[gx, gz] = meshgrid(0:n-1, 0:n-1);
r = [gx(:) - (n-1)/2, zeros(Q,1), gz(:) - (n-1)/2] * lam/2;   % centred grid
As = hop(false); Ad = hop(nf_dist > 0);

    function A = hop(nearfield)
        tau = tau_n * ds; pw = 10.^(p_db/10); pw = pw / sum(pw);
        L = numel(tau);
        g = sqrt(pw) .* (randn(1, L) + 1j*randn(1, L)) / sqrt(2);
        if los, g(1) = sqrt(pw(1)) * exp(1j*2*pi*rand); end    % LoS tap: fixed amplitude
        u = randn(L, 3); u(:,2) = abs(u(:,2)); u = u ./ vecnorm(u, 2, 2);
        ph = exp(-1j*k*(r*u.'));                             % Q x L plane-wave phases
        if nearfield && los
            pu = nf_dist * u(1, :);                          % user position for the LoS tap
            d = vecnorm(r - pu, 2, 2);
            ph(:, 1) = (nf_dist ./ d) .* exp(-1j*k*(d - nf_dist));
        end
        A = (ph .* g) * exp(-1j*2*pi*tau(:) * (bins(:).'*df));
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
