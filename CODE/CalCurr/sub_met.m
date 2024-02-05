%%% Metabolism
function met = sub_met(Tp,Tb,tdif,wgt,param)
    %Tp: pelagic temp
    %Tb: bottom temp
    %tdif: frac pelagic time
    %wgt: ind weight of size class
    %fcrit: feeding level to meet resting respiration rate
    %cmax: max consumption rate
    %U: swimming speed

    temp = (Tp.*tdif) + (Tb.*(1.0-tdif));

    %Own Fn ------------
    %Metabolism with its own coeff, temp-sens, mass-sens
    met = (exp(param.kt * (temp-10.0)) .* param.amet .* wgt.^(-param.bpow)) ./365.0;

end
