%%% Type II consumption
function con = sub_cons(param,Tp,Tb,tpel,wgt,enc)
    %Tp: pelagic temp
    %Tb: bottom temp
    %tpel: frac pelagic time
    %wgt: ind weight of size class
    %enc: array of all encountered food
    % calculates consumption rate of first element of enc
    
    temp = (Tp.*tpel) + (Tb.*(1.0-tpel));
    
    %Cmax rate
    cmax = (exp(param.kc * (temp-10.0)) .* param.h .* wgt^(-param.bcmx)) ./365.0;
    
    ENC = sum(enc,2); % total biomass encountered
    con = cmax .* enc(:,1) ./ (cmax + ENC); % Type II
    
end
