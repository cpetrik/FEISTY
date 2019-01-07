%%% Update biomass
function [bio_sm,bio_md,predS,predM] = sub_update_be2(BE,det,CC,bio_in,con,fish)
    %bio_in = benthic biomass
    %con = biomass specific consumption rate by MD & LD
    %fish = biomass of MD & LD
    
    %! Peterson & Wroblewski: daily & uses dry weight
    %nmort = exp(0.063*(temp-15.0)) * 5.26e-3 * (wgt/9.0)^(-0.25);
    % Use Wei et al mean (1.3956 gWW) or median (0.7656) wgt
%     wgt = 0.7656;
%     omort = exp(0.063*(Tb-15.0)) * 5.26e-3 * (wgt/9.0)^(-0.25);
%     eaten = con.*bio;
%     pred = sum(eaten);
%     bio_out = bio_in - sum(eaten) - bio_in*omort;
    
    %! Quadratic mortality from carrying capacity
    eaten = con.*fish;
    pred = sum(eaten,2);
    
    %Half of detritus becomes small, half becomes medium benthic inverts
    det2 = det.*0.5;
    
    % Chemostat 1
    %r = BE*det;
    %bio_out = bio_in + r * (CC - bio_in) - pred;
    
    % Chemostat 2
%     r = BE*det2;
%     bio_out = bio_in + r.*(CC - bio_in) - eaten;
    
    % Logistic 1 
    %r = BE*det ./ bio_in; %Needs to be in units of per time (g/m2/d) * (g/m2)
    %bio_out = bio_in + r .* bio_in .* (1 - bio_in./CC) - pred;
    
    % Logistic 2
    bio_out = bio_in;
    for i=1:2
        r = BE*det2 ./ bio_in(:,i);
        bio_out(:,i) = bio_in(:,i) + r .* bio_in(:,i) .* ...
            (1 - bio_in(:,i)./CC) - eaten(:,i);
    end
    
    bio_sm = bio_out(:,1);
    bio_md = bio_out(:,2);
    predS = eaten(:,1);
    predM = eaten(:,2);
    
end
