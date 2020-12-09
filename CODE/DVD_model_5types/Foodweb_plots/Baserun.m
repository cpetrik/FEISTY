
%%% simple run with parameters as specified in baseparameters
baseparameters
% Default initial conditions:
param.y0 = [0.1*param.K param.B0];
if(param.bottom <= param.photic)
  param.y0(param.ix1(2):param.ix2(2))=0; %mesopelagics to zero
  param.y0(param.ix1(4):param.ix2(4))=0; %bathypelagics to zero
end
BeR = (13*(param.bottom/param.photic)^-0.86);
BeR(BeR>=13)=13;
param.r = [1   1   1   1];  % g ww/m2/yr
param.K = [5   5    BeR*0.9   BeR*0.1];
result = poem(param);
plotPoemf(param, result)
depthdistribution(param)
plotdiet(param,result)

y = result.y;
Bin = floor(0.6*length(y));
BiosZ = mean(y(Bin:end,1));
BiolZ = mean(y(Bin:end,2));
BioMf = sum(mean(y(Bin:end,param.ix1(1):param.ix2(5))));
 
%%% loop run in which the resource productivity of the benthic channel is increased

bRp = [0.01 linspace(0.1,10,30)]; %linspace(0.1,5,30); % change in resource
%bRp = linspace(50,3000,80); % change with depth

va = length(bRp);
BioFf = zeros(va,1);
BioMf = zeros(va,1);
BioPf = zeros(va,1);
BioBaP = zeros(va,1);
BioDf = zeros(va,1);
BioBaD = zeros(va,1);

for i = 1:va
    baseparameters
    
    param.bottom = 100;  % depth in meter
    BeR = (.1*(param.bottom/400)^-0.86)*0.1;
    BeR(BeR>=.01)=.01;
    param.r =  [bRp(i)    bRp(i)    BeR*0.9    BeR*0.1];  % g ww/m2/yr
    
    % feeding preference matrix
    [param.theta, param.depthDay, param.depthNight, param.avlocDay, param.avlocNight] = calcpreference(param); % feeding preference matrix 

    
    param.avdepth = (param.avlocDay + param.avlocNight)/2;
    param.avdepth(param.avdepth < 200) = 200;
    param.critfrac = exp(-0.25*log(param.avdepth))/exp(-0.25*log(200))*0.2; % ikeda 2016 Routine metabolic rates of pelagic marine fishes
    param.fc(5:end) = param.critfrac(5:end);
    
    % mortality stratified along the water column
    %param.light = 1*exp(-0.01*param.depth(param.ixFish)); % describes the amount of ligth at depth
    %param.mort0 = param.mort0 + 0.5*(param.extmort(param.ixFish)'.*param.wc(param.ixFish).^(-0.25)').*param.light'; 

    result = poem(param);
    y = result.y;
    Bin = floor(0.8*length(y));
       
    BioFf(i) = sum(mean(y(Bin:end,param.ix1(1):param.ix2(1))));
    BioMf(i) = sum(mean(y(Bin:end,param.ix1(2):param.ix2(2))));
    BioPf(i) = sum(mean(y(Bin:end,param.ix1(3):param.ix2(3))));
    BioBaP(i) = sum(mean(y(Bin:end,param.ix1(4):param.ix2(4))));
    BioDf(i) = sum(mean(y(Bin:end,param.ix1(5):param.ix2(5))));
    BioBaD(i) = sum(mean(y(Bin:end,param.ix1(6):param.ix2(6))));
   end

plot(bRp,BioFf,'Color',[0 0.447 0.741],'linewidth',2) 
hold on 
%plot(bRp,BioMf,'Color',[0.929 0.694 0.125],'linewidth',2)
plot(bRp,BioPf,'Color',[0.85 0.325 0.098],'linewidth',2)
plot(bRp,BioDf,'Color',[0.494 0.184 0.556],'linewidth',2)
plot(bRp,BioBaP,'Color',[0.466 0.674 0.188],'linewidth',2)
plot(bRp,BioBaD,'Color',[0  0 0],'linewidth',2)
plot(bRp,BiosDe,'Color',[0  .7 .7],'linewidth',2)
legend('Forage fish','mesopel','Pelagic predator','Demersal pred','Bathypel','Bathydem','sm dem')
ylabel('Biomass')
xlabel('Pel resource prod')
hold off
   
%depth plot 
plot(bRp,BioFf,'Color',[0 0.447 0.741],'linewidth',2) 
hold on 
%plot(bRp,BioMf,'Color',[0.929 0.694 0.125],'linewidth',2)
plot(bRp,BioPf,'Color',[0.85 0.325 0.098],'linewidth',2)
plot(bRp,BioDf,'Color',[0.494 0.184 0.556],'linewidth',2)
plot(bRp,BioBaP,'Color',[0.466 0.674 0.188],'linewidth',2)
plot(bRp,BioBaD,'Color',[0  0 0],'linewidth',2)
plot(bRp,BiosDe,'Color',[0  .7 .7],'linewidth',2)
legend('Forage fish','mesopel','Pelagic predator','Demersal pred','Bathypel','Bathydem','Small dem')
ylabel('Biomass')
xlabel('Seafloor depth (m)')
ax = gca;
hold off
