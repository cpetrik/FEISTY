%%% Update biomass
function [bio_out, pred] = sub_update_be(Tb,bio_in,param,det,con,bio)
    %Tb = bottom temperature
    %bio_in = benthic biomass
    %det = detritus flux to bottom (g/m2/d)
    %con = biomass specific consumption rate by MD & LD
    %bio = biomass of MD & LD
    
    BE = param.bent_eff;
    eaten = con.*bio;
    pred = sum(eaten,2);
        
    %! No carrying capacity
    r = BE .* det; %Needs to be in units of per time (g/m2/d) * (g/m2)
    bio_out = bio_in + r - pred;
    
    %! Quadratic mortality from carrying capacity
    % Chemostat
%     r = BE*det;
%     bio_out = bio_in + r * (r*CC - bio_in) - pred;
    
    % Logistic
%     r = BE*det ./ bio_in; %Needs to be in units of per time (g/m2/d) * (g/m2)
%     bio_out = bio_in + r .* bio_in .* (1 - bio_in./CC) - pred;

end
