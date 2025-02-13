%POEM size classes
 
clear all
close all

%%

zmlow = 0.2;
zmhi = 2;
zllow = 2; 
zlhi = 20;
jlow = 20;
jhi = 200;
alow = 200;
ahi = 2e3;

lzm_s = linspace(log10(zmlow),log10(zmhi),3); 
lzl_s = linspace(log10(zllow),log10(zlhi),3);
lj_s = linspace(log10(jlow),log10(jhi),3);
la_s = linspace(log10(alow),log10(ahi),3);

%Length
zm_L = 10^(lzm_s(2)); % = 10^((log10(0.2)+log10(2))/2)
zl_L = 10^(lzl_s(2));
j_L = 10^(lj_s(2));
a_L = 10^(la_s(2));

%% Weight from Andersen & Beyer weight-length relationship for fish,
% from Watkins et al for zooplankton
% in g
%Med zoop
zm_w = exp(1.953 + 2.399 * log(zm_L)) * 1e-6;
%Large zoop
zl_w = exp(1.953 + 2.399 * log(zl_L)) * 1e-6;
%Larvae
l_w = 0.01 * (0.1*zl_L)^3;
%Juveniles piscivores/Adult planktivores
j_w = 0.01 * (0.1*j_L)^3;
%Adult piscivores
a_w = 0.01 * (0.1*a_L)^3;

