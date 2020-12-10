
%%% simple run with parameters as specified in baseparameters
baseparameters
% Default initial conditions:
param.y0 = [0.1*param.K 0.01*param.B0];
if(param.bottom <= param.mesop)
  param.y0(param.ix1(2):param.ix2(2))=0; %mesopelagics to zero
  param.y0(param.ix1(4):param.ix2(4))=0; %mid-water pred to zero
end
bent = 150; %max detrital flux
BeR = 0.075*(bent*(param.bottom/param.photic)^-0.86); %BE * detrital flux from Martin curve
BeR(BeR>=(bent*0.075)) = bent*0.075;
param.K =  [80    80    BeR    0];  % resource carrying capacity in g C ww/m2
result = poem(param);
plotPoemf(param, result)
plotdiet(param,result)