%%% Type II consumption
function con = sub_cons(Tp,Tb,tpel,wgt,enc)
    %Tp: pelagic temp
    %Tb: bottom temp
    %tpel: frac pelagic time
    %wgt: ind weight of size class
    %enc: array of all encountered food
    % calculates consumption rate of first element of enc
    
    global cfn h bcmx
    
    temp = (Tp.*tpel) + (Tb.*(1.0-tpel));
    
    %Cmax rates from other models
    if (cfn==2)
        % Hartvig et al (g/g/day) ref to 10C
        cmax = (exp(0.063*(temp-10.0)) .* 85.0 .* wgt^(-0.25)) ./365.0;
    elseif (cfn==3)
        % mizer (g/g/day) ref to 10C
        cmax = (exp(0.063*(temp-10.0)) .* 40.0 .* wgt^(-0.25)) ./365.0;
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
        cmax = (exp(0.063*(temp-10.0)) .* h .* wgt^(-bcmx)) ./365.0;
    end
    
    ENC = sum(enc,2); % total biomass encountered
    con = cmax .* enc(:,1) ./ (cmax + ENC); % Type II
    
end
