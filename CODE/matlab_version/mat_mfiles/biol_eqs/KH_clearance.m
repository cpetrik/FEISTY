% Kiorboe & Hirst eqs
% Play with clearance rate for MF

clear all
close all

L_s = 10^((log10(2)+log10(20))/2);
L_m = 10^((log10(20)+log10(200))/2);
L_l = 10^((log10(200)+log10(2000))/2);

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;

%% Kiorboe & Hirst
% T=15C
temp=15;

%Max clearance rate/search volume
FmaxLarv = exp(0.063*(temp-15.0)) .* 10^(8.1) .* M_s^(0.057) .* (24e-3/9);
FmaxS = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_s^(-0.24) .* (24e-3/9);
FmaxM = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_m^(-0.24) .* (24e-3/9);
FmaxL = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_l^(-0.24) .* (24e-3/9);

m=[M_s; M_m; M_l];
fmax=[FmaxS;FmaxM;FmaxL];

figure(1)
%subplot(2,2,1)
plot(log10(m),log10(fmax),'.k','MarkerSize',25); hold on;
%plot(log10(M_s),log10(FmaxLarv),'.b','MarkerSize',25); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
%legend('K&H','K&H larv')
%legend('location','southwest')
title('Specific clearance rate')
%print('-dpng','Clearance_comp_mass.png')
print('-dpng','Clearance_mass.png')

%% vs temp
temp=-2:30;

%Max clearance rate/search volume
FmaxLarv = exp(0.063*(temp-15.0)) .* 10^(8.1) .* M_s^(0.057) .* (24e-3/9);
FmaxS = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_s^(-0.24) .* (24e-3/9);
FmaxM = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_m^(-0.24) .* (24e-3/9);
FmaxL = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_l^(-0.24) .* (24e-3/9);


%%

figure(2)
%subplot(2,2,1)
plot(temp,log10(FmaxS),'r','LineWidth',2); hold on;
plot(temp,log10(FmaxM),'k','LineWidth',2); hold on;
plot(temp,log10(FmaxL),'b','LineWidth',2); hold on;
%plot(temp,log10(FmaxLarv),'m--','LineWidth',2); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('temp (C)')
%legend('K&H S','K&H M','K&H L','K&H larv')
legend('K&H S','K&H M','K&H L')
legend('location','west')
title('Specific clearance rate')
%print('-dpng','Clearance_comp_temp.png')
print('-dpng','Clearance_temp.png')


%% K&H data
load('K&H14_clearance_rates.mat')
massC=mass;
cl_day = clearance;
cl_hr = clearance * 24;
clear mass clearance
load('K&H14_feeding_rates.mat')
Fmax_day = Fmax / 24;
Fmax_hr = Fmax;
sFmax_day = SpecFmax / 24;
sFmax_hr = SpecFmax;

W=[massC;mass*1e-3];
F=[cl_day;Fmax_day];
sF=[cl_day./massC; Fmax_day./(mass*1e-3)];

log10M=[-9:4]';
M=10.^(log10M);
temp=15;
larFmax = exp(0.063*(temp-15.0)) .* 10^(8.1) .* M.^(0.057) .* (24e-3/9);
adFmax = exp(0.063*(temp-15.0)) .* 1.74e3 .* M.^(-0.24) .* (24e-3/9);

%%
figure
%subplot(2,2,1)
plot(log10(massC),log10(cl_day),'o')
xlabel('log10 mass (gC)')
ylabel('log10 clearance')

figure
%subplot(2,2,1)
plot(log10(massC),log10(cl_day./massC),'o')
xlabel('log10 mass (gC)')
ylabel('log10 specific clearance')

figure
%subplot(2,2,1)
plot(log10(mass),log10(Fmax_day),'o')
xlabel('log10 mass (mgC)')
ylabel('log10 Fmax')

figure
%subplot(2,2,1)
plot(log10(mass),log10(sFmax_day),'o')
xlabel('log10 mass (mgC)')
ylabel('log10 specific Fmax')
%%
figure
%subplot(2,2,1)
plot(log10(W),log10(F),'o')
xlabel('log10 mass allW (gC)')
ylabel('log10 clearance allW')

%%
figure
%subplot(2,2,1)
plot(log10(W),log10(sF),'o'); hold on;
%plot(log10(M),log10(larFmax),'b'); hold on;
%plot(log10(M),log10(adFmax),'k'); hold on;
xlabel('log10 mass allW (gC)')
ylabel('log10 specific clearance allW')
legend('Data','Larval eq','Adult eq')
