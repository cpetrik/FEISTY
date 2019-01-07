%%%% Fishing
function [bio, caught, fmort] = sub_quad_fishing(bio,F,selec,Tp,Tb,tdif)
    %bio = fish biomass
    %F = fishing rate per day
    %selec = fishery selectivity 
    %NOTE: selec only 1 or 0 now, but could update code (here & gamma calc) so it is a fraction
     
    % Quadratic fishing mortality
    %Temp-dep
%     temp = (Tp.*tdif) + (Tb.*(1.0-tdif));
%     caught = (exp(0.063*(temp-10.0))) .* bio.^2 .* selec .* F;
    %Temp-indep
    caught = bio.^2 .* selec .* F;
    fmort = caught ./ bio;
    bio = bio - caught;
    

end
