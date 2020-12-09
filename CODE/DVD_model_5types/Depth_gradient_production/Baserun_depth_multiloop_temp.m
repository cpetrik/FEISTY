% check how changes with depth and productivity depend on temperature

%% tropical
clearvars

bRp = [5     50    100];    % change with plankton production
bRb = [130   250   380];    % change with benthic detritus flux
idx = [1      1      1];    % temperature region

bdepth = round(linspace(50,3000,119));  % 119 change with depth (needs to be a whole number)

rowle = length(bdepth);
colle = length(bRb); 

BioFf = zeros(rowle,colle);
BioMf = zeros(rowle,colle);
BioPf = zeros(rowle,colle);
BioBaP = zeros(rowle,colle);
BioDf = zeros(rowle,colle);

for j = 1:colle
    
for i = 1:length(bdepth)
    baseparameters
    param.bottom = bdepth(i);  % depth in meter
       
    % calculate resource productivities
    bent = bRb(j);
    BeR = 0.1*(bent*(param.bottom/param.photic)^-0.86);
    BeR(BeR>=(bent*0.1)) = bent*0.1;
    param.K =  [bRp(j)  bRp(j)    BeR*1    BeR*0];  % g C ww/m2/yr
    
      % starting values
    if (i==1)
        % Default initial conditions:
        param.y0 = [0.1*param.K 0.1*param.B0];
    else
        % Use "results" argument for initial conditions:
        param.y0 = result.y(end,:) + 0.01;
    end
    
    % exclude mesopelagics and bathypelagics in shallow waters
    if(param.bottom <= param.mesop)
        param.y0(param.ix1(2):param.ix2(2))=0; %mesopelagics to zero
        param.y0(param.ix1(4):param.ix2(4))=0; %bathypelagics to zero
    end
    
    % feeding preference matrix
    [param.theta, param.depthDay, param.depthNight, param.avlocDay, param.avlocNight] = calcpreference(param); % feeding preference matrix 

    % temperature
    param.region = idx(j); % 1=tropical, 2=temperate, 3=boreal, 4 = no temp scaling
    [param.scTemp, param.scTempm] = calctemperature(param);
    param.Cmax = (param.h*param.wc.^param.n)./param.wc .* param.scTemp;
    param.V = (param.gamma*param.wc.^param.q)./param.wc.* param.scTemp;
    param.Mc = (param.met*param.wc.^param.m)./param.wc .* param.scTempm;
    
    result = poem(param);
    y = result.y;
    Bin = floor(0.8*length(y));
       
    BioFf(i,j) = sum(mean(y(Bin:end,param.ix1(1):param.ix2(1))));
    BioMf(i,j) = sum(mean(y(Bin:end,param.ix1(2):param.ix2(2))));
    BioPf(i,j) = sum(mean(y(Bin:end,param.ix1(3):param.ix2(3))));
    BioBaP(i,j) = sum(mean(y(Bin:end,param.ix1(4):param.ix2(4))));
    BioDf(i,j) = sum(mean(y(Bin:end,param.ix1(5):param.ix2(5))));
end

end
save depth_tempscaling_tropical

%% temperate
clearvars

bRp = [5     50    100];    % change with plankton production
bRb = [130   250   380];    % change with benthic detritus flux
idx = [2      2      2];    % temperature region

bdepth = round(linspace(50,3000,119));  % 119 change with depth (needs to be a whole number)

rowle = length(bdepth);
colle = length(bRb); 

BioFf = zeros(rowle,colle);
BioMf = zeros(rowle,colle);
BioPf = zeros(rowle,colle);
BioBaP = zeros(rowle,colle);
BioDf = zeros(rowle,colle);

for j = 1:colle
    
for i = 1:length(bdepth)
    baseparameters
    param.bottom = bdepth(i);  % depth in meter
       
    % calculate resource productivities
    bent = bRb(j);
    BeR = 0.1*(bent*(param.bottom/param.photic)^-0.86);
    BeR(BeR>=(bent*0.1)) = bent*0.1;
    param.K =  [bRp(j)  bRp(j)    BeR*1    BeR*0];  % g C ww/m2/yr
    
      % starting values
    if (i==1)
        % Default initial conditions:
        param.y0 = [0.1*param.K 0.1*param.B0];
    else
        % Use "results" argument for initial conditions:
        param.y0 = result.y(end,:) + 0.01;
    end
    
    % exclude mesopelagics and bathypelagics in shallow waters
    if(param.bottom <= param.mesop)
        param.y0(param.ix1(2):param.ix2(2))=0; %mesopelagics to zero
        param.y0(param.ix1(4):param.ix2(4))=0; %bathypelagics to zero
    end
    
    % feeding preference matrix
    [param.theta, param.depthDay, param.depthNight, param.avlocDay, param.avlocNight] = calcpreference(param); % feeding preference matrix 

    % temperature
    param.region = idx(j); % 1=tropical, 2=temperate, 3=boreal, 4 = no temp scaling
    [param.scTemp, param.scTempm] = calctemperature(param);
    param.Cmax = (param.h*param.wc.^param.n)./param.wc .* param.scTemp;
    param.V = (param.gamma*param.wc.^param.q)./param.wc.* param.scTemp;
    param.Mc = (param.met*param.wc.^param.m)./param.wc .* param.scTempm;
    
    result = poem(param);
    y = result.y;
    Bin = floor(0.8*length(y));
       
    BioFf(i,j) = sum(mean(y(Bin:end,param.ix1(1):param.ix2(1))));
    BioMf(i,j) = sum(mean(y(Bin:end,param.ix1(2):param.ix2(2))));
    BioPf(i,j) = sum(mean(y(Bin:end,param.ix1(3):param.ix2(3))));
    BioBaP(i,j) = sum(mean(y(Bin:end,param.ix1(4):param.ix2(4))));
    BioDf(i,j) = sum(mean(y(Bin:end,param.ix1(5):param.ix2(5))));
end

end
save depth_tempscaling_temperate

%% boreal
clearvars

bRp = [5     50    100];    % change with plankton production
bRb = [130   250   380];    % change with benthic detritus flux
idx = [3      3      3];    % temperature region

bdepth = round(linspace(50,3000,119));  % 119 change with depth (needs to be a whole number)

rowle = length(bdepth);
colle = length(bRb); 

BioFf = zeros(rowle,colle);
BioMf = zeros(rowle,colle);
BioPf = zeros(rowle,colle);
BioBaP = zeros(rowle,colle);
BioDf = zeros(rowle,colle);

for j = 1:colle
    
for i = 1:length(bdepth)
    baseparameters
    param.bottom = bdepth(i);  % depth in meter
       
    % calculate resource productivities
    bent = bRb(j);
    BeR = 0.1*(bent*(param.bottom/param.photic)^-0.86);
    BeR(BeR>=(bent*0.1)) = bent*0.1;
    param.K =  [bRp(j)  bRp(j)    BeR*1    BeR*0];  % g C ww/m2/yr
    
      % starting values
    if (i==1)
        % Default initial conditions:
        param.y0 = [0.1*param.K 0.1*param.B0];
    else
        % Use "results" argument for initial conditions:
        param.y0 = result.y(end,:) + 0.01;
    end
    
    % exclude mesopelagics and bathypelagics in shallow waters
    if(param.bottom <= param.mesop)
        param.y0(param.ix1(2):param.ix2(2))=0; %mesopelagics to zero
        param.y0(param.ix1(4):param.ix2(4))=0; %bathypelagics to zero
    end
    
    % feeding preference matrix
    [param.theta, param.depthDay, param.depthNight, param.avlocDay, param.avlocNight] = calcpreference(param); % feeding preference matrix 

    % temperature
    param.region = idx(j); % 1=tropical, 2=temperate, 3=boreal, 4 = no temp scaling
    [param.scTemp, param.scTempm] = calctemperature(param);
    param.Cmax = (param.h*param.wc.^param.n)./param.wc .* param.scTemp;
    param.V = (param.gamma*param.wc.^param.q)./param.wc.* param.scTemp;
    param.Mc = (param.met*param.wc.^param.m)./param.wc .* param.scTempm;
    
    result = poem(param);
    y = result.y;
    Bin = floor(0.8*length(y));
       
    BioFf(i,j) = sum(mean(y(Bin:end,param.ix1(1):param.ix2(1))));
    BioMf(i,j) = sum(mean(y(Bin:end,param.ix1(2):param.ix2(2))));
    BioPf(i,j) = sum(mean(y(Bin:end,param.ix1(3):param.ix2(3))));
    BioBaP(i,j) = sum(mean(y(Bin:end,param.ix1(4):param.ix2(4))));
    BioDf(i,j) = sum(mean(y(Bin:end,param.ix1(5):param.ix2(5))));
end

end
save depth_tempscaling_boreal
