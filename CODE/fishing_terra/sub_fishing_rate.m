%%%% Fishing
function [bio, caught, fmort] = sub_fishing_rate(bio,F,selec)
    %bio = fish biomass
    %F = fishing rate per day
    %selec = fishery selectivity 
        
    % Linear fishing mortality
    caught = bio .* selec .* F;
    fmort = caught ./ bio;
    bio = bio - caught;
    
end
