%%% Consumption/Cmax
function clev = sub_clev(con,Tp,Tb,tdif,wgt)
    % calculates consumption rate of first element of enc
    
    global cfn h bcmx kc
    
    temp = (Tp.*tdif) + (Tb.*(1.0-tdif));
    
    %Cmax rates from other models
    if (cfn==2)
        % Hartvig et al (g/g/day) ref to 10C
        cmax = (exp(kc*(temp-10.0)) .* 85.0 .* wgt^(-0.25)) ./365.0;
    elseif (cfn==3)
        % mizer (g/g/day) ref to 10C
        cmax = (exp(kc*(temp-10.0)) .* 40.0 .* wgt^(-0.25)) ./365.0;
    elseif (cfn==4)
        % J&C15 (g/g/day) ref to 10C
        temp2 = temp+273;
        Tref = 283;
        E=0.6;
        k=8.62e-5;
        tfact = exp((-1*E/k)*((1./temp2)-(1./Tref)));
        cmax = (tfact .* 25.0 .* wgt^(-0.33)) ./365.0;
    else
    %Cmax rate
        cmax = (exp(kc*(temp-10.0)) .* h .* wgt^(-bcmx)) ./365.0;
    end
    
    %Clev
    clev = con./cmax;
end
