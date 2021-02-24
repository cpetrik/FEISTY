%%% ENERGY AVAILABLE FOR SOMATIC GROWTH
function gam = sub_gamma_fishvec(param,fish)
    % d = predation loss
    % nmort = natural mortality rate
    % Frate = fishing mortality rate
    % selec = harvested selectivity (adults 100%, juveniles 10%)
    K = param.kappa;
    z = param.Z(param.ixFish);
    nu = fish.nu(:,param.ixFish);
    d = fish.die(:,param.ixFish);
    B = fish.bio(:,param.ixFish);
    nmrt = fish.nmort(:,param.ixFish);
    Frate = param.dfrate;
    selec = param.sel;
    
    Z = repmat(z',param.NX,1);
    kap = repmat(K,param.NX,1);
    
    % convert predation mortality to biomass specific rate
    D = (d./B) + nmrt + (Frate.*selec);
    
    gg = ((kap.*nu) - D) ./ (1 - (Z.^(1 - (D ./ (kap.*nu)))));
    gam = min(gg,nu);
    gam = max(gam,0);
    gam(isnan(gam)) = zeros;
    
end
