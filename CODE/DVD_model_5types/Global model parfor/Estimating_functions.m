clearvars

% double click RUN_TOT_global_final
load('globaldat_input.mat')
load('globaldat_temp.mat')

baseparameters
rowle = length(datcell);
colle = length(param.ixFish);
fluxin = zeros(rowle,colle);
foodfrac = zeros(rowle,2);

for gridi = 1:rowle
    bottom = depthWOA(gridi); % depth in meter
    photicm = photic(gridi);
    [theta, depthDay, depthNight] = calcpreference(param,bottom,photicm); % feeding preference matrix 

    % exclude mesopelagics and bathypelagics in shallow waters
    % if(bottom <= param.mesop)
    %    y0(param.ix1(2):param.ix2(2))=0; %mesopelagics to zero
    %    y0(param.ix1(4):param.ix2(4))=0; %bathypelagics to zero
    % end
    
    [scTemp, scTempm] = calctemperature(param,gridi,temp_grad,temp_prof,depthDay,depthNight,bottom);
    Cmax = (param.h*param.wc.^param.n)./param.wc .* scTemp;
    V = (param.gamma*param.wc.^param.q)./param.wc .* scTemp;
    Mc = (param.met*param.wc.^param.m)./param.wc .* scTempm;
        
    y = datcell{gridi,1}';
    [f, mortpr, mortpred, Eavail] = calcEncounter_functions(y, param,V,theta,Cmax,Mc);
    
    %
    % Fish:
    %
    ixFish = param.ixFish;

    % Total mortality:
    mort = mortpred(ixFish)' + param.mort0 + param.F;

    % Flux out of the size group:
    v = (Eavail(param.ixFish))';%./param.wc(param.ixFish))';
    vplus = max(0,v);
    gamma = (param.kappa'.*vplus - mort) ./ ...
        (1 - param.z(param.ixFish).^(1-mort./(param.kappa'.*vplus)) );

    % get flux into the size class
    fluxin (gridi,:) = gamma' .* y(param.ixFish)';
    
    % get pelagic/benthic food demersal fish
    mortpr_flux = mortpr .* y';
    pel = [1:2 5:25];
    bent = 3:4;
    
    foodfrac(gridi,1) = sum(sum(mortpr_flux(25:30,pel)));
    foodfrac(gridi,2) = sum(sum(mortpr_flux(25:30,bent)));
end

csvwrite('foodfrac.csv',foodfrac)
csvwrite('flux_in.csv', fluxin)