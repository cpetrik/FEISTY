%%% Metabolism
function met = sub_met_dem(Tp,Tb,tdif,wgt)
    %Tp: pelagic temp
    %Tb: bottom temp
    %tdif: frac pelagic time
    %wgt: ind weight of size class
    %Dact: activity reduction for demersals
    
    global kt bpow amet Dact
    
    temp = (Tp.*tdif) + (Tb.*(1.0-tdif));
    
    %Own Fn ------------
    %Metabolism
    %shares coeff "h" with cmax
    %met = fcrit .* (exp(kt*(temp-10.0)) .* h .* wgt.^(-bpow)) ./365.0;
    
    %its own coeff amet
    met = (exp(kt*(temp-10.0)) .* Dact .* amet .* wgt.^(-bpow)) ./365.0;

end
