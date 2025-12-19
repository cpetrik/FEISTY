%%%% Fishing
function [caught, fmort] = sub_fishing_rate(bio,F,selec)
    %bio = fish biomass
    %F = fishing rate per day
    %selec = fishery selectivity 
        
    % Linear fishing mortality
    caught = bio .* selec .* F;

    % Prevent divide by zero if no fish in that layer
    fmort = zeros(size(bio));
    id = bio>0;
    fmort(id) = caught(id) ./ bio(id);
    
    %bio = bio - caught; Don't update bio; Do in dB/dt calc instead
    
end
