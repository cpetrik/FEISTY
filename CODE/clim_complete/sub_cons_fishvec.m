%%% Type II consumption
function con = sub_cons_fishvec(cmax,fish)
    %enc: array of all encounter rate of all food
    %cmax: maximum consumption rate
    
    ENC = sum(fish.enc,3); % total biomass encountered
    con = cmax .* fish.enc ./ (cmax + ENC); % Type II
    
end
