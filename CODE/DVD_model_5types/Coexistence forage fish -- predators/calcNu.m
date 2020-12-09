function [v, nu] = calcNu(Eavail, mort, param)

v = (Eavail(param.ixFish))';
vplus = max(0,v);
nu = (vplus - mort) ./ (1 - param.z(param.ixFish).^(1-mort./vplus) );