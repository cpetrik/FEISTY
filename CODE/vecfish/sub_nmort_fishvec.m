%%% Temp-dep natural mortality
function nmort = sub_nmort_fishvec(param,fish,temp)
    
    MORT = param.MORT;
    Nat_mrt = param.Nat_mrt;
    NX = param.NX;
    wgt = param.wc;
    
    if (MORT==0) % None
        nmort = zeros(size(fish.bio));
    end
    if (MORT==1) % Constant
        nmort = Nat_mrt * ones(size(fish.bio));
    end
    if (MORT==2) % Hartvig Temperature-dependent mortality
        nmort = exp(0.063*(temp-10.0)) .* 0.84 .* wgt.^(-0.25) /365.0;
    end
    if (MORT==3) % mizer Temperature-dependent mortality
        % mizer
        nmort = exp(0.063*(temp-10.0)) .* 3.0 .* wgt.^(-0.25) /365.0;
    end
    if (MORT==4) % Jennings & Collingridge Temperature-dependent mortality
        temp2 = temp+273;
        Tref = 283;
        E=0.6;
        k=8.62e-5;
        tfact = exp((-1*E/k)*((1./temp2)-(1./Tref)));
        nmort = tfact .* 0.5 .* wgt.^(-0.33) /365.0;
    end
    if (MORT==5) % Peterson & Wrob Temperature-dependent mortality
        % Peterson & Wroblewski (daily & uses dry weight)
        nmort = exp(0.063*(temp-15.0)) * 5.26e-3 * (wgt./9.0)^(-0.25); 
    end
    if (MORT==6) % Temp-dep but constant by weight\
        nmort = exp(0.063*(temp-10.0)) .* Nat_mrt;
    end
    if (MORT==7) % wgt-dep but constant by temp
        nmort = 0.5 .* wgt.^(-0.25) /365.0;
    end
end
