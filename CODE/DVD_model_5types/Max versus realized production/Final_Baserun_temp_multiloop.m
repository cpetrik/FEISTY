clearvars

bRp = [0.001 0.01 0.1 0.5 linspace(1,200,75)];  % change with plankton Rmax
rowle = length(bRp);
temp = [0 10 20 30];
TQ10 =  1.88.^((temp-10)/10);
TQ10m = 1.88.^((temp-10)/10);
colle = length(TQ10); 

BioFf = zeros(rowle,colle);
BioPf = zeros(rowle,colle);
zoopS = zeros(rowle,colle);
zoopL = zeros(rowle,colle);

for j = 1:length(TQ10)

for i = 1:length(bRp)
    baseparameters
           
    % calculate resource productivities
    bent = 1;
    BeR = (bent*(param.bottom/param.photic)^-0.86);
    BeR(BeR>=bent) = bent;
    param.K =  [bRp(i)  bRp(i)    BeR*0.9    BeR*0.1];  % g C ww/m2
    param.y0 = [0.1*param.K 0.1*param.B0]; % Default initial conditions:
    
    % calculate temperature and food dependency
    param.Cmax = (param.h*param.wc.^param.n)./param.wc .* TQ10(j);
    param.Mc = (param.met*param.wc.^param.m)./param.wc .* TQ10m(j);
    param.V =(param.gamma*param.wc.^param.q)./param.wc .* TQ10(j); 
    
    result = poem(param);
    y = result.y;
    Bin = floor(0.8*length(y));
       
    BioFf(i,j) = sum(mean(y(Bin:end,param.ix1(1):param.ix2(1))));
    BioPf(i,j) = sum(mean(y(Bin:end,param.ix1(3):param.ix2(3))));
    zoopS(i,j) = sum(mean(y(Bin:end,1)));
    zoopL(i,j) = sum(mean(y(Bin:end,2)));
end
end

M = horzcat(BioFf, BioPf, zoopL, zoopS,bRp');
csvwrite('pel_prod_temp.csv',M)
