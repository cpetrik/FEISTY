%%%  Encounter rates
function enc = sub_enc(param,Tp,wgt,prey,td,pref)
    % Tp: pelagic temp
    % Tb: bottom temp
    % wgt: ind weight of size class
    % pred: pred biomass density,
    % prey: prey biomass density,
    % A: predator search rate,
    % tpel: time spent in pelagic,
    % tprey: time spent in area with that prey item -not needed
    % pref: preference for prey item

    %Mult enc rate by 100 to go from 100m vert integral to per m3
    
    temp = (Tp);
    
    %Enc rate
    A = (exp(param.ke * (temp-10.0)) .* 100.0*param.gam .* wgt^(-param.benc)) ./365.0;
    
    %Encounter per predator, mult by biomass later
    % frac = zeros(param.NX,1);
    % ID = (tprey>0);
    % frac(ID) = 1.0;
    
    % g/m2 * g/g/d * -- * --
    %enc = td.*prey.*A.*frac.*pref;
    enc = td.*prey.*A.*pref;

end
