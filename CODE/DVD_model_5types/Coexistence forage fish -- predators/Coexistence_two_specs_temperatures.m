%% coexistence forage fish and large pelagics
clearvars

% interactions between spel and lpel

% define temperatures
temp = [0   10   20   25   30];
TQ10_1 =  1.88.^((temp-10)/10);
TQ10_2 =  1.88.^((temp-10)/10); %2.35

% Default initial conditions:
bRp = linspace(1, 150, 150); % change with pelagic production
bRb = 0.1;  % change with benthic detritus flux
rowle = length(bRp);
colle = length(temp); 

baseparameters 
param.B0(param.ix1(2)-4:param.ix2(2)-4)=0; %put small mesopelagics to zero
param.B0(param.ix1(4)-4:param.ix2(4)-4)=0; %put bathypelagics to zero
param.B0(param.ix1(5)-4:param.ix2(5)-4)=0; %put large demersals to zero
BioFf = zeros(rowle,colle);
BioPf = zeros(rowle,colle);

for jj = 1:length(temp)
    
    for j = 1:length(bRp)
        param.K =  [bRp(j)   bRp(j)   bRb*1   bRb*0]; 
        param.y0 = [0.1*param.K param.B0.*20];
    
        % calculate temperature dependency
        param.Cmax = (param.h*param.wc.^param.n)./param.wc .* TQ10_1(jj);
        param.Mc = (param.met*param.wc.^param.m)./param.wc .* TQ10_2(jj);
        param.V =(param.gamma*param.wc.^param.q)./param.wc .* TQ10_1(jj); 
      
        result = poem(param);
        y = result.y;
        Bin = floor(0.8*length(y));
        BioFf(j,jj) = sum(mean(y(Bin:end,param.ix1(1):param.ix2(1))));
        BioPf(j,jj) = sum(mean(y(Bin:end,param.ix1(3):param.ix2(3))));
      end
end
  save('scenario_2spec.mat')
 