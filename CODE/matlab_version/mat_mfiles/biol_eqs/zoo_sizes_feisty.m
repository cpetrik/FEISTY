% Notes on zooplankton size 
% from POEMv1 developed by James Watson (Z_s)
% Modifications by Colleen Petrik (M_zm, M_zl)

%% Median Zooplankton body mass 
% 2.96 + 2.73ln(ESD) McCauly 1984, James Watkins, 
% Lars Rudstam and Kristen Holeck later
% Zm ESD = 1.1mm, Zl ESD = 11mm (from Charlie/COBALT)
%const global Z_s = [exp(2.96 + (2.73*log(1.1)))*1.0e-6; 
%					exp(2.96 + (2.73*log(11)))*1.0e-6]; 
%const global Z_s = [exp(1.953 + (2.399*log(1.1)))*1.0e-6; 
%					exp(1.953 + (2.399*log(11)))*1.0e-6]; 
Z_s = [exp(1.953 + (2.399*log(2)))*1.0e-6;...
    exp(1.953 + (2.399*log(20)))*1.0e-6];


%% ! Fish mass from length using Andersen & Beyer 2013
% Convert from mm to cm and use their const coeff = 0.01g/cm3
M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;

%% ! Median Zooplankton body mass
% James Watkins, Lars Rudstam and Kristen Holeck
% (from Charlie/COBALT)
% L-W eq is ug dry weight, *1e-6 to g, *9 to wet weight
L_zm = 10^((log10(0.2)+log10(2))/2); % lengths (ESD)
L_zl = 10^((log10(2)+log10(20))/2);
M_zm = exp(1.953 + (2.399*log(L_zm)))*1.0e-6*9; % body mass
M_zl = exp(1.953 + (2.399*log(L_zl)))*1.0e-6*9;

