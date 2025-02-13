%%% Consumption/Cmax
function clev = sub_clev(param,con,Tp,Tb,tdif,wgt)
    % calculates consumption rate of first element of enc
    
    temp = (Tp.*tdif) + (Tb.*(1.0-tdif));
    
    %Cmax rate
    cmax = (exp(param.kc * (temp-10.0)) .* param.h .* wgt^(-param.bcmx)) ./365.0;
    
    %Clev
    clev = con./cmax;
end
