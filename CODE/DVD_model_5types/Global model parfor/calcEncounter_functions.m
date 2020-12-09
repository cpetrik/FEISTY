function [f, mortpr, mortpred, Eavail] = calcEncounter_functions(y,param,V,theta,Cmax,Mc)
%
% Consumption
%
for i = 1:length(y)
    Enc(i,:) = V(i) .* theta(i,:).*y';
end

Encspecies = sum(Enc');
f = Encspecies./(Cmax+Encspecies);
f(isnan(f)) = 0;
Eavail = Cmax.*param.epsAssim.*f - Mc;

%
% Mortality:
%
 for i = 1:length(y)
 mortpr(i,:) = (Cmax(i) .* V(i) .* theta(i,:) ./(Cmax(i) + Encspecies(i))) * y(i);
 end

 mortpred = sum(mortpr);

%for i = 1:length(y)
%   mortpred(i) = sum( Enc(:,i).*(y).*(1-f)' )/(eps*eps + y(i));
%end

