%%% Temp-dep natural mortality
function nmort = sub_nmort(param,Tp,tpel,wgt)
    %Tp: pelagic temp
    %Tb: bottom temp
    %tpel: frac pelagic time
    
    MORT = param.MORT;
    Nat_mrt = param.Nat_mrt;
    NX = param.NX;
    
    if (MORT==0) % None
        nmort = 0.0;
    end
    if (MORT==1) % Constant
        nmort = Nat_mrt * ones(NX,1);
    end
    if (MORT==2) % Hartvig Temperature- & weight-dependent mortality
        temp = (Tp);
        % Hartvig
        nmort = exp(0.063*(temp-10.0)) .* 0.84 .* wgt^(-0.25) /365.0;
    end
    if (MORT==3) % mizer Temperature- & weight-dependent mortality
        temp = (Tp);
        % mizer
        nmort = exp(0.063*(temp-10.0)) .* 3.0 .* wgt^(-0.25) /365.0;
    end
    if (MORT==4) % Jennings & Collingridge Temperature- & weight-dependent mortality
        temp = (Tp);
        % J&C
        temp2 = temp+273;
        Tref = 283;
        E=0.6;
        k=8.62e-5;
        tfact = exp((-1*E/k)*((1./temp2)-(1./Tref)));
        nmort = tfact .* 0.5 .* wgt^(-0.33) /365.0;
    end
    if (MORT==5) % Peterson & Wrob Temperature- & weight-dependent mortality
        temp = (Tp);
        % Peterson & Wroblewski (daily & uses dry weight)
        nmort = exp(0.063*(temp-15.0)) * 5.26e-3 * (wgt/9.0)^(-0.25); 
    end
    if (MORT==6) % Temp-dep but constant by weight
        temp = (Tp);
        nmort = exp(0.063*(temp-10.0)) .* Nat_mrt;
    end
    if (MORT==7) % wgt-dep but constant by temp
        nmort = 0.5 .* wgt^(-0.25) /365.0;
    end

    nmort = nmort .* tpel;

end
