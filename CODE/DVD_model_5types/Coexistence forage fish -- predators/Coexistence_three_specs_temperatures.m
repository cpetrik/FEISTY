clearvars

% analyze how inclusion of demersal fish and benthic food web change
% coexistence dynamics of forage fish and pelagic predator

% define temperatures
temp = [0   10   20  25  30];
TQ10_1 =  1.88.^((temp-10)/10);
TQ10_2 =  2.35.^((temp-10)/10); %2.35; 1.88

% Default initial conditions:
bRb = linspace(1, 500, 50)*0.1; % flux "carrying capacity" 
bRp = 100;  % pelagic carrying capacity
rowle = length(bRb);
colle = length(temp); 

baseparameters
param.B0(param.ix1(2)-4:param.ix2(2)-4)=0; %put small mesopelagics to zero
param.B0(param.ix1(4)-4:param.ix2(4)-4)=0; %put bathypelagics to zero
param.y0 = [0.1*param.K param.B0];

BioFf = zeros(rowle,colle);
BioPf = zeros(rowle,colle);
BioDem = zeros(rowle,colle);
for jj = 1:length(temp)
   
    % calculate temperature dependency
    param.Cmax = (param.h*param.wc.^param.n)./param.wc .* TQ10_1(jj);
    param.Mc = (param.met*param.wc.^param.m)./param.wc .* TQ10_2(jj);
    param.V =(param.gamma*param.wc.^param.q)./param.wc .* TQ10_1(jj); 

    for j = 1:length(bRb)
        param.K =  [bRp   bRp   bRb(j)*1   bRb(j)*0]; 
        result = poem(param);
        y = result.y;
        Bin = floor(0.8*length(y));
        BioFf(j,jj) = sum(mean(y(Bin:end,param.ix1(1):param.ix2(1))));
        BioPf(j,jj) = sum(mean(y(Bin:end,param.ix1(3):param.ix2(3))));
        BioDem(j,jj) = sum(mean(y(Bin:end,param.ix1(5):param.ix2(5))));
    end
end
  save('scenario_3spec_unequal_70both.mat')
