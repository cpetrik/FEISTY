%%% Consumption/Cmax
function clev = sub_clev(param,con,Tp,tdif,wgt)
    % calculates consumption rate of first element of enc
    
    temp = (Tp);
    
    %Cmax rate
    cmax = (exp(param.kc * (temp-10.0)) .* param.h .* wgt^(-param.bcmx)) ./365.0;
    
    %Clev
    clev = tdif .* con./cmax;
end
