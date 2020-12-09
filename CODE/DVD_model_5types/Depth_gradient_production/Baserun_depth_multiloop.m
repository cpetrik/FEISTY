
%% run with depth, similar to six classes in:
% baserun_depth_multiloop_varying_number_of_classes

bRp = [5     50    100];    % change with plankton production
bRb = [130   250   380];    % change with benthic detritus flux
idx = [4      4      4];    % temperature region

bdepth = round(linspace(50,3000,119));  % 119 change with depth (needs to be a whole number)
nstage = 6; 

rowle = length(bdepth);
colle = length(bRb); 

BioFf = zeros(rowle,colle);
BioMf = zeros(rowle,colle);
BioPf = zeros(rowle,colle);
BioBaP = zeros(rowle,colle);
BioDf = zeros(rowle,colle);
peltodem = zeros(rowle,colle);
bentodem = zeros(rowle,colle);

for j = 1:colle
    
for i = 1:length(bdepth)
    baseparameters_nstage
    param.bottom = bdepth(i);  % depth in meter
    
    % starting values
    if (i==1)
        % Default initial conditions:
        param.y0 = [0.1*param.K param.B0];
    else
        % Use "results" argument for initial conditions:
        param.y0 = result.y(end,:) + 0.01;
    end
    
    % exclude mesopelagics and bathypelagics in shallow waters
    if(param.bottom < 200)
        param.y0(param.ix1(2):param.ix2(2))=0; %mesopelagics to zero
        param.y0(param.ix1(4):param.ix2(4))=0; %bathypelagics to zero
    end
    
    % calculate resource productivities
    bent = bRb(j);
    BeR = 0.1*(bent*(param.bottom/param.photic)^-0.86);
    BeR(BeR>=(bent*0.1)) = bent*0.1;
    param.K =  [bRp(j)  bRp(j)    BeR*1    BeR*0];  % g C ww/m2/yr
    
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
    
    yend = mean(y(Bin:end,:));

    % calculate predation interactions
    %[mortpr_flux] = calcEncounter_foodwebplot(yend', param);
    
    %pel = [1:2 5:25];
    %bent = 3:4;
    
    %peltodem(i,j) = sum(sum(mortpr_flux(25:30,pel)));
    %bentodem(i,j) = sum(sum(mortpr_flux(25:30,bent)));
     
end

end
save depth_multiloop_fluxes