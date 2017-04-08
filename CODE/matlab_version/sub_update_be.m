%%% Update biomass
function [bio_out, pred] = sub_update_be(bio_in,BE,det,CC,con,bio)
    %bio_in = benthic biomass
    %con = biomass specific consumption rate by MD & LD
    %bio = biomass of MD & LD
    
    %! Peterson & Wroblewski: daily & uses dry weight
    %nmort = exp(0.063*(temp-15.0)) * 5.26e-3 * (wgt/9.0)^(-0.25);
    % Use Wei et al mean (1.3956 gWW) or median (0.7656) wgt
%     wgt = 0.7656;
%     omort = exp(0.063*(Tb-15.0)) * 5.26e-3 * (wgt/9.0)^(-0.25);
%     eaten = con.*bio;
%     pred = sum(eaten);
%     bio_out = bio_in - sum(eaten) - bio_in*omort;
    
    %! Quadratic mortality from carrying capacity
    eaten = con.*bio;
    pred = sum(eaten,2);
    
    %Half of detritus becomes small, half becomes medium benthic inverts
    %det2 = det.*0.5;
    
    % Chemostat
    %r = BE*det;
    %bio_out = bio_in + r * (r*CC - bio_in) - pred;
    
    % Logistic
    r = BE*det ./ bio_in; %Needs to be in units of per time (g/m2/d) * (g/m2)
    bio_out = bio_in + r .* bio_in .* (1 - bio_in./CC) - pred;
end
