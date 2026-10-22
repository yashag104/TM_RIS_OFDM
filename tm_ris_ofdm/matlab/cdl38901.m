function [As, Ad, info] = cdl38901(Q, bins, df, fc, model, ds, opt)
%CDL38901 3GPP TR 38.901 CDL channels (A, C, D) for both RIS hops, with cluster angles.
%   Cluster delays (x ds), powers, AoD/AoA/ZoD/ZoA and intra-cluster spreads from the 38.901 CDL tables
%   (as distributed with NVIDIA Sionna v0.19.0, cdl/CDL-*.json). Each NLoS cluster has 20 rays with the
%   standard 38.901 ray offsets (Table 7.5-3) scaled by the cluster angle spreads, random pairing of the
%   azimuth/zenith offsets and random ray phases; CDL-D's first cluster is the specular LoS ray.
%   BS -> RIS hop: element phases from the ARRIVAL angles (the RIS receives); RIS -> user hop: from the
%   DEPARTURE angles (the RIS transmits). A random azimuth rotation of the whole cluster set per draw gives
%   channel diversity. Elements on a sqrt(Q) x sqrt(Q) lambda/2 grid in the x-z plane (normal along y).
%   opt.squint (false): array phase evaluated at the actual frequency fc + f (beam squint) instead of fc.
%   opt.nearfield_m (0): if > 0, the LoS ray of the RIS -> user hop uses exact element distances at that range,
%                        evaluated at each subcarrier frequency.
if nargin < 7, opt = struct(); end
if ~isfield(opt, 'squint'), opt.squint = false; end
if ~isfield(opt, 'nearfield_m'), opt.nearfield_m = 0; end
persistent tabs
if isempty(tabs), tabs = struct(); end
key = ['m' upper(model)];
if ~isfield(tabs, key)
    here = fileparts(mfilename('fullpath'));
    tabs.(key) = jsondecode(fileread(fullfile(here, 'cdl', ['CDL-' upper(model) '.json'])));
end
tb = tabs.(key);
c0 = 3e8; lam = c0/fc; n = round(sqrt(Q));
[gx, gz] = meshgrid(0:n-1, 0:n-1);
r = [gx(:) - (n-1)/2, zeros(Q,1), gz(:) - (n-1)/2] * lam/2;
fr = bins(:).' * df;                                       % 1 x M baseband frequencies
alpha = [0.0447 -0.0447 0.1413 -0.1413 0.2492 -0.2492 0.3715 -0.3715 0.5129 -0.5129 ...
         0.6797 -0.6797 0.8844 -0.8844 1.1481 -1.1481 1.5195 -1.5195 2.1551 -2.1551];
As = hop('aoa', 'zoa', tb.cASA, tb.cZSA, false);
Ad = hop('aod', 'zod', tb.cASD, tb.cZSD, opt.nearfield_m > 0);
info = struct('model', model, 'ds', ds, 'clusters', tb.num_clusters, 'los', tb.los);

    function A = hop(az, ze, cAS, cZS, nearfield)
        tau = tb.delays(:).' * ds; pw = 10.^(tb.powers(:).'/10); pw = pw / sum(pw);
        rot = 360*rand;                                    % random azimuth rotation of the cluster set
        A = zeros(Q, numel(fr));
        for c = 1:numel(tau)
            if tb.los && c == 1, nr = 1; aoff = 0; zoff = 0;
            else, nr = 20; aoff = cAS*alpha; zoff = cZS*alpha(randperm(20)); end
            for m = 1:nr
                phi = deg2rad(tb.(az)(c) + rot + aoff(min(m, numel(aoff))));
                th = deg2rad(tb.(ze)(c) + zoff(min(m, numel(zoff))));
                u = [sin(th)*cos(phi), sin(th)*sin(phi), cos(th)];
                g = sqrt(pw(c)/nr) * exp(1j*2*pi*rand);
                proj = r*u.';                              % Q x 1
                if nearfield && tb.los && c == 1
                    d0 = opt.nearfield_m; pu = d0*u; d = vecnorm(r - pu, 2, 2);
                    ph = (d0 ./ d) .* exp(-1j*2*pi*(fc + fr).*(d - d0)/c0);
                elseif opt.squint
                    ph = exp(-1j*2*pi*(fc + fr).*proj/c0);
                else
                    ph = exp(-1j*2*pi*fc*proj/c0) * ones(1, numel(fr));
                end
                A = A + g * ph .* exp(-1j*2*pi*tau(c)*fr);
            end
        end
    end
end
