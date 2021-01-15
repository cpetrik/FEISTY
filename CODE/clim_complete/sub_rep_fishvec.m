%%% BIOMASS MADE FROM REPRODUCTION
function [gamma, rep] = sub_rep_fishvec(param,fish)
%nu: energy for growth or spawning
%K: proportion allocated to growth
    kap = param.kappa;
    nu = fish.nu(:,param.ixFish);
    gamma = fish.gamma(:,param.ixFish);
    
    K = repmat(kap,param.NX,1);

    % Growth of last size class included in repro
    %NOTE: Still never going to accumulate biomass as muscle tissue
    nu0 = max(nu,0.0);
    rho = (1.0-K) .* nu0;  %energy available for from eating
    rep = rho;
    rep(rho>0) = rho(rho>0) + gamma(rho>0);
    gamma(rho>0) = zeros;

end
