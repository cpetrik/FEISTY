function dydt = poem_deriv(t,y,param,V,theta,Cmax,Mc,K,tEnd)
y(y<0) = 0;
R = y(param.ixR);
B = y(param.ixFish);
%
% Feeding
%
[f, mortpred, Eavail] = calcEncounter(y, param,V,theta,Cmax,Mc);
%
% Fish:
%
ixFish = param.ixFish;

% Total mortality:
mort = mortpred(ixFish)' + param.mort0 + param.F;

% Flux out of the size group:
v = (Eavail(param.ixFish))';%./param.wc(param.ixFish))';
vplus = max(0,v);
gamma = (param.kappa'.*vplus - mort) ./ ...
    (1 - param.z(param.ixFish).^(1-mort./(param.kappa'.*vplus)) );
Fout = gamma.*B;
Repro = (1-param.kappa').*vplus.*B;

% Flux into the size group (where Fout in last stage stays in last stage):
%for i = 1:param.nSpecies
%    ix = (param.ix1(i):param.ix2(i)) - length(R);
%    %ixPrev = [ix(end) ix(1:(end-1))];
%    Fin(ix) = [0; Fout(1:(end-1))];
%    Fin(ix(end)) = Fin(end)+Fout(end);
%    % Reproduction:
%    Fin(ix(1)) = param.eRepro(i)* Repro(ix(end));
%end

% Flux into the size group(where Fout in last stage stays is included in repro):
for i = 1:param.nSpecies
    ix = (param.ix1(i):param.ix2(i)) - length(R);
    ixPrev = [ix(end) ix(1:(end-1))];
    Fin(ix) = Fout(ixPrev);
    % Reproduction:
    Fin(ix(1)) = param.eRepro(i)*(Fin(ix(1)) + Repro(ix(end)));
end

dBdt = Fin' - Fout + (v - mort).*B - Repro;
% 
% Resource
%
dRdt = param.r'.*(K'-R) - mortpred(param.ixR)'.*R;

dydt = [dRdt; dBdt];

