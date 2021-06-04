%%% Metabolism
function met = sub_met_midw(Tp,Tm,Tb,tdif,wgt,param)
    %Tp: epipel temp
    %Tm: mesopel temp
    %Tb: bottom temp
    %tdif: frac time (1) epi, (2) meso, (3) benthic
    %wgt: ind weight of size class

    temp = (Tp.*tdif(:,1)) + (Tm.*tdif(:,2)) + (Tb.*tdif(:,3));

    %Own Fn ------------
    %Metabolism with its own coeff, temp-sens, mass-sens
    met = (exp(param.kt * (temp-10.0)) .* param.amet .* wgt.^(-param.bpow)) ./365.0;

end
