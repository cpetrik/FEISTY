clearvars

load('globaldat_input.mat')
load('globaldat_temp.mat')
idx = length(RmaxL);

datcell = cell(idx,1);

parpool(20)

baseparameters
parfor gridi = 1:idx 
    K =  [RmaxS(gridi)  RmaxL(gridi)  ben_prod(gridi)  0*ben_prod(gridi)];  % g ww/m2
    y0 = [0.1*K param.B0];
    
    bottom = depthWOA(gridi); % depth in meter
    photicm = photic(gridi);
    [theta, depthDay, depthNight] = calcpreference(param,bottom,photicm); % feeding preference matrix 

    % exclude mesopelagics and bathypelagics in shallow waters
    if(bottom <= param.mesop)
        y0(param.ix1(2):param.ix2(2))=0; %mesopelagics to zero
        y0(param.ix1(4):param.ix2(4))=0; %bathypelagics to zero
    end
    
    [scTemp, scTempm] = calctemperature(param,gridi,temp_grad,temp_prof,depthDay,depthNight,bottom);
    Cmax = (param.h*param.wc.^param.n)./param.wc .* scTemp;
    V = (param.gamma*param.wc.^param.q)./param.wc .* scTemp;
    Mc = (param.met*param.wc.^param.m)./param.wc .* scTempm;
    
    result = poem(param,y0,V,theta,Cmax,Mc,K);
    y = result.y;
    Bin = floor(0.8*length(y));
    Bio = mean(y(Bin:end,:));
    datcell{gridi,1} = Bio;
end

save RUN_TOT_global_200614