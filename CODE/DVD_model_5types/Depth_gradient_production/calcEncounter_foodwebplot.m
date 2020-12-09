function [mortpr_flux] = calcEncounter_foodwebplot(y, param)
%
% Consumption
%
for i = 1:length(y)
    Enc(i,:) = param.V(i) .* param.theta(i,:).*y';
end

Encspecies = sum(Enc');
f = Encspecies./(param.Cmax+Encspecies);
f(isnan(f)) = 0;
Eavail = param.Cmax.*param.epsAssim.*f - param.Mc;

%
% Mortality:
%
for i = 1:length(y)
mortpr(i,:) = (param.Cmax(i) .* param.V(i) .* param.theta(i,:) ./(param.Cmax(i) + Encspecies(i))) * y(i);
end

%mortpred = sum(mortpr);

mortpr_flux = mortpr .* y';
%for i = 1:length(y)
%   mortpred(i) = sum( Enc(:,i).*(y).*(1-f)' )/(eps*eps + y(i));
%end

