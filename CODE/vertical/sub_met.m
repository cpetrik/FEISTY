%%% Metabolism
function met = sub_met(Tp,td,wgt,param)
    %Tp: pelagic temp
    %Tb: bottom temp
    %tdif: frac pelagic time
    %wgt: ind weight of size class
    %fcrit: feeding level to meet resting respiration rate
    %cmax: max consumption rate
    %U: swimming speed

    temp = Tp;

    %Own Fn ------------
    %Metabolism with its own coeff, temp-sens, mass-sens
    met = td .* (exp(param.kt * (temp-10.0)) .* param.amet .* wgt.^(-param.bpow)) ./365.0;

end
