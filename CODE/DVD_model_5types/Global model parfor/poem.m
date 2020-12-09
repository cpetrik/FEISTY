function result = poem(param,y0,V,theta,Cmax,Mc,K, result)
 
%
% Run:
%
[t,y] = ode45(@poem_deriv, 0:1:param.tEnd, y0, odeset('NonNegative',1:length(y0)),param,V,theta,Cmax,Mc,K);
%
% Construct output:
%
result.y = y;
result.R = y(:,param.ixR);
result.B = y(:,param.ixFish);
result.t = t;

result.Yield = result.B .* (ones(length(t),1)*param.F');

[result.f, result.mortpred, result.Eavail] = calcEncounter(y(end,:)', param,V,theta,Cmax,Mc);
result.mort = result.mortpred(param.ixFish)' + param.mort0 + param.F;
[result.v , result.nu] = calcNu(result.Eavail, result.mort, param);

