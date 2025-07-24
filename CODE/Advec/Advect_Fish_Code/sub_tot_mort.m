%%% Update biomass
function tmort = sub_tot_mort(pred,nmort,fmort)
    % all inputs are in g g-1 d-1; 
    % nmort = natural mortality rate
    % fmort = fishing rate
    % pred = predation rate
    
    tmort = nmort + fmort + pred;
end