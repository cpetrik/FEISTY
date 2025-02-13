%%%  Encounter rates
function enc = sub_enc(param,Tp,Tb,wgt,prey,tpel,tprey,pref)
    % Tp: pelagic temp
    % Tb: bottom temp
    % wgt: ind weight of size class
    % pred: pred biomass density,
    % prey: prey biomass density,
    % A: predator search rate,
    % tpel: time spent in pelagic,
    % tprey: time spent in area with that prey item.
    % pref: preference for prey item
    
    temp = (Tp.*tpel) + (Tb.*(1.0-tpel));
    
    %Enc rate
    A = (exp(param.ke * (temp-10.0)) .* param.gam .* wgt^(-param.benc)) ./365.0;
    
    %Encounter per predator, mult by biomass later
    frac = zeros(param.NX,1);
    ID = (tprey>0);
    frac(ID) = 1.0;
    
    % g/m2 * g/g/d * -- * --
    enc = prey.*A.*frac.*pref;
end
