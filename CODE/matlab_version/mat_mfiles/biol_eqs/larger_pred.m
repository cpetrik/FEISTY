% Compare my equations with Kiorboe & Hirst

clear all
close all

L_l = 10^((log10(200)+log10(2000))/2);
L_L = 10^((log10(2000)+log10(20000))/2);
L_2 = 2000;

M_l = 0.01 * (0.1*L_l)^3;
M_L = 0.01 * (0.1*L_L)^3;
M_2 = 0.01 * (0.1*L_2)^3;

% Kiorboe & Hirst
% T=15C
temp=15;

%Max ingestion
Imaxl = exp(0.063*(temp-15.0)) .* 2.5 .* M_l^(-0.51) .* 24e-3;
ImaxL = exp(0.063*(temp-15.0)) .* 2.5 .* M_L^(-0.51) .* 24e-3;
Imax2 = exp(0.063*(temp-15.0)) .* 2.5 .* M_2^(-0.51) .* 24e-3;

%Max clearance rate/search volume
Fmaxl = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_l^(-0.24) .* (24e-3/9);
FmaxL = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_L^(-0.24) .* (24e-3/9);
Fmax2 = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_2^(-0.24) .* (24e-3/9);

%Respiration
respl = 0.1*Imaxl;
respL = 0.1*ImaxL;
resp2 = 0.1*Imax2;

%% Higher predator consump
prey = [1e-6;5e-6;1e-5;5e-5;1e-4;5e-4;1e-3;5e-3;1e-2;5e-2;1e-1;5e-1;1;5;1e2;5e2];
A = 1.74e3 .* M_2^(-0.24) .* (24e-3/9);
cmax = 2.5 .* M_2^(-0.51) .* 24e-3;
enc = A .* 0.1 .* prey.^2;
con = exp(0.063*(temp-15.0)) .* (cmax .* A .* 0.1 .* prey.^2) ...
    ./ (cmax + A .* 0.1 .* prey.^2);

%approximate consumption
ntcon = (cmax .* A .* 0.1 .* prey.^2) ./ (cmax + A .* 0.1 .* prey.^2);
acon = ntcon ./ prey.^2;
acon2 = exp(0.063*(temp-15.0)) .* 0.03 .* prey.^2; 
k = 6.1e-3;
con3 = exp(0.063*(temp-15.0)) .* (cmax .* prey.^2) ./ (k + prey.^2);
con4 = exp(0.063*(temp-15.0)) .* (cmax .* prey) ./ (k*10 + prey);

figure(1)
plot(log10(prey),con,'k'); hold on;
plot(log10(prey),acon2,'r'); hold on;
plot(log10(prey),con3,'--b'); hold on;
plot(log10(prey),con4,'m'); hold on;
ylim([0 2e-4])

%% In code
Khp = 6.1e-2;
nmort = cmax .* prey ./ (Khp + prey);
nmortB = nmort .* prey;
Nat_mrt = (0.44 / 365) * ones(size(prey));

figure(2)
plot(log10(prey),nmortB,'k'); hold on;
plot(log10(prey),Nat_mrt,'--r'); hold on;

