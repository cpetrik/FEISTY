%%% ENERGY AVAILABLE FOR SOMATIC GROWTH
function gamma = sub_gamma(K,Z,nu,d,B,nmrt,Frate,selec)
    % d = predation loss
    % nmort = natural mortality rate
    % Frate = fishing mortality rate
    % selec = harvested selectivity (adults 100%, juveniles 10%)
    
    % convert predation mortality to biomass specific rate
    if (selec > 0)
        D = (d./B) + nmrt + (Frate*selec);
    else
        D = (d./B) + nmrt;
    end
    kap=K;
    gg = ((kap.*nu) - D) ./ (1 - (Z.^(1 - (D ./ (kap.*nu)))));
    gamma = min(gg,nu);
    gneg = (gg < 0);
    gamma(gneg) = 0.0;
    gnan = (isnan(gg));
    gamma(gnan) = 0.0;
end
