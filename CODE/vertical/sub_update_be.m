%%% Update biomass
function [bio_out, pred] = sub_update_be(bio_in,param,det,con,bio)
    %Tb = bottom temperature
    %bio_in = benthic biomass
    %det = detritus flux to bottom (g/m2/d)
    %con = biomass specific consumption rate by MD & LD
    %bio = biomass of MD & LD
    
    BE = param.bent_eff;
    CC = param.CC;
    eaten = con.*bio;
    pred = sum(eaten,2);
    pred = pred(param.NX); %bottom value only

    bio_out = bio_in;
    bent = bio_in(param.NX); %bottom value only
       
    if (CC==0)
        %! No carrying capacity
        r = BE .* det; %Needs to be in units of per time (g/m2/d) * (g/m2)
        bio_out(param.NX) = bent + r - pred;
    else
        % Logistic
        r = BE*det ./ bent; %Needs to be in units of per time (g/m2/d) * (g/m2)
        bio_out(param.NX) = bent + r .* bent .* (1 - bent./CC) - pred;
    end
    
    %! Quadratic mortality from carrying capacity
    % Chemostat
%     r = BE*det;
%     bio_out = bio_in + r * (r*CC - bio_in) - pred;

end
