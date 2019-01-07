%%% Metabolism
function met = sub_met(Tp,Tb,tdif,wgt)
    %Tp: pelagic temp
    %Tb: bottom temp
    %tdif: frac pelagic time
    %wgt: ind weight of size class
    %fcrit: feeding level to meet resting respiration rate
    %cmax: max consumption rate
    %U: swimming speed
    
    global mfn kt bpow amet
    
    temp = (Tp.*tdif) + (Tb.*(1.0-tdif));
       
    %Own Fn ------------
    %Metabolism from other models
    if (mfn==2)
        % Hartvig et al (g/g/day) k=10 ref to 10C
        met = (exp(0.063*(temp-10.0)) .* 10 .* wgt^(-0.25)) ./365.0;
    elseif (mfn==3)
        % mizer et al (g/g/day) k=4.8 ref to 10C
        met = (exp(0.063*(temp-10.0)) .* 4.8 .* wgt^(-0.25)) ./365.0;
    elseif (mfn==4)
        % J&C15 et al (g/g/day) k=2.5 ref to 10C
        temp2 = temp+273;
        Tref = 283;
        E=0.6;
        k=8.62e-5;
        tfact = exp((-1*E/k)*((1./temp2)-(1./Tref)));
        met = (tfact .* 2.5 .* wgt^(-0.25)) ./365.0;
    else
    %Metabolism with its own coeff, temp-sens, mass-sens
        met = (exp(kt*(temp-10.0)) .* amet .* wgt.^(-bpow)) ./365.0;
    end
    
end
