%%%  Encounter rates
function enc = sub_enc(Tp,Tb,wgt,prey,tpel,tprey,pref)
    % Tp: pelagic temp
    % Tb: bottom temp
    % wgt: ind weight of size class
    % pred: pred biomass density,
    % prey: prey biomass density,
    % A: predator search rate,
    % tpel: time spent in pelagic,
    % tprey: time spent in area with that prey item.
    % pref: preference for prey item
    
    global efn gam NX benc ke
    
    temp = (Tp.*tpel) + (Tb.*(1.0-tpel));
    
    %Enc rates from other models
    if (efn==1)
        % Specific clearance rates from Kiorboe & Hirst (m3/g/day)
        % divide by 100m to put in m2/g/day b/c zoop is integrated over 100m depth
        A = (exp(ke.*(temp-15.0)) .* 10^(3.24) .* wgt^(-0.24)) .* (24e-3/9); %./ 100.0;
    elseif (efn==2)
        % Hartvig et al (m3/g/day) gamma = 0.8e4 ref to 10C
        A = (exp(ke*(temp-10.0)) .* 0.8e4 .* wgt^(0.8-1.0)) ./365.0;
    elseif (efn==3)
        % mizer et al (m3/g/day) gamma = 2.9e3? ref to 10C
        A = (exp(ke*(temp-10.0)) .* 2.9e3 .* wgt^(0.8-1.0)) ./365.0;
    elseif (efn==4)
        % J&C15 et al (m3/g/day) gamma = 2.9e3? ref to 10C
        A = (exp(ke*(temp-10.0)) .* 1.37e4 .* wgt^(0.9-1.0)) ./365.0;
    else
    %Enc rate
        A = (exp(ke*(temp-10.0)) .* gam .* wgt^(-benc)) ./365.0;
    end
    
    %Encounter per predator, mult by biomass later
    frac = zeros(NX,1);
    ID = (tprey>0);
    frac(ID) = 1.0;
    
    enc = prey.*A.*frac.*pref;
end
