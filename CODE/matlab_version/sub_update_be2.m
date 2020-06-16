%%% Update biomass
function [bio_sm,bio_md,predS,predM] = sub_update_be2(BE,det,bio_in,Dcon,Dbio)
    global fracSB
    %BE = benthic efficiency
    %det = seafloor detritus flux in g/m2/d
    %bio_in = benthic biomass in g
    %Dcon = biomass specific consumption rate by MD & LD in /m2/d
    %Dbio = biomass of MD & LD in g
    
    eaten = Dcon.*Dbio;
%     pred = sum(eaten,2);
    
    %! Temp-dep
    %Incr w/T
%     r = exp(0.063*(Tb-10.0)) .* BE .* det;
    %Decr w/T
%     mx = exp(0.063*(33.1-10.0));
%     r = (mx - exp(0.063*(Tb-10.0))) .* BE .* det;
    %Opt T
    
    %! No carrying capacity
    r = BE .* det; %Needs to be in units of per time (g/m2/d) * (g/m2)
    
    bio_sm = bio_in(:,1) + (fracSB)*r - eaten(:,1);
    bio_md = bio_in(:,2) + (1 - fracSB)*r - eaten(:,2);
    
    predS = eaten(:,1);
    predM = eaten(:,2);

end
