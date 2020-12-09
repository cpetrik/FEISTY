
coord = [18928 36916 26410 18271 19713]; % get region 1-5

for j = 1:5
    clearvars param bottom photicm depthDay depthNight y Bin yend
    %% load temp data and resources + depth
    load('globaldat_input.mat')
    load('globaldat_temp.mat')

    gridi = coord(j);
    bottom = depthWOA(gridi); 
    photicm = photic(gridi);
    baseparameters

    param.K =  [RmaxS(gridi)  RmaxL(gridi)  1*ben_prod(gridi)  0*ben_prod(gridi)];  % g ww/m2
    param.y0 = [0.1*param.K param.B0];
    [param.theta, depthDay, depthNight, param.avlocDay, param.avlocNight] = calcpreference(param,bottom,photicm); % feeding preference matrix 

    % exclude mesopelagics and bathypelagics in shallow waters
        if(bottom < photicm)
            param.y0(param.ix1(2):param.ix2(2))=0; %mesopelagics to zero
            param.y0(param.ix1(4):param.ix2(4))=0; %bathypelagics to zero
        end

    [scTemp, scTempm] = calctemperature(param,gridi,temp_grad,temp_prof,depthDay,depthNight,bottom);
    param.Cmax = (param.h*param.wc.^param.n)./param.wc .* scTemp;
    param.V = (param.gamma*param.wc.^param.q)./param.wc .* scTemp;
    param.Mc = (param.met*param.wc.^param.m)./param.wc .* scTempm;

    result = poem(param);

    % get results
    y = result.y;
    Bin = floor(0.8*length(y));
    yend = mean(y(Bin:end,:));

    % calculate predation interactions
    [mortpr_flux] = calcEncounter_foodwebplot(yend', param);

    clearvars temp_grad temp_prof photic depthWOA RmaxL ben_prod RmaxS
    save(['region_' num2str(j)])
end