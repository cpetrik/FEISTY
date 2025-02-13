% Calc growth efficiencies expected for our fish sizes

clear all
close all

L_s = 10^((log10(2)+log10(20))/2);
L_m = 10^((log10(20)+log10(200))/2);
L_l = 10^((log10(200)+log10(2000))/2);
L_max = 10^((log10(2000)+log10(20000))/2);

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;
M_max = 0.01 * (0.1*L_max)^3;
m=[M_s; M_m; M_l];

%Hartvig et al. constants
q=0.8;
n=0.75;
beta=10;
alpha=0.7;

%Avg growth effic
eibar = (beta ^ (2*n-q-1)) / (q+2-2*n);

%% If asymp size is Adult
%P&D
eipd = alpha .* (1 - (m./M_l).^(1-n));
%F
eiF = alpha .* (1 - (m(1:2)./M_m).^(1-n));

%% If asymp size is next size up
eiS = alpha .* (1 - (M_s./M_m).^(1-n));
eiM = alpha .* (1 - (M_m./M_l).^(1-n));
eiL = alpha .* (1 - (M_l./M_max).^(1-n));