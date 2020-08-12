%%% Update biomass
function [bio_out, pred] = sub_update_be(Tb,bio_in,param,det,con,bio)
    %bio_in = benthic biomass
    %con = biomass specific consumption rate by MD & LD
    %bio = biomass of MD & LD
    
    BE = param.bent_eff;
    
    eaten = con.*bio;
    pred = sum(eaten,2);
    
    %! Temp-dep
    %Incr w/T
%     r = exp(0.063*(Tb-10.0)) .* BE .* det;
    %Decr w/T
%     mx = exp(0.063*(33.1-10.0));
%     r = (mx - exp(0.063*(Tb-10.0))) .* BE .* det;
    
    %! No carrying capacity
    r = BE .* det; %Needs to be in units of per time (g/m2/d) * (g/m2)
    bio_out = bio_in + r - pred;
    
end
