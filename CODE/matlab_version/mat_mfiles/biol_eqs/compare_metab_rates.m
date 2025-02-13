% Compare diff Cmax, metab, A, ingest eqs

clear all
close all

L_s = 10^((log10(2)+log10(20))/2);
L_m = 10^((log10(20)+log10(200))/2);
L_l = 10^((log10(200)+log10(2000))/2);

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;
m=[M_s; M_m; M_l];

%Hartvig et al. constants
q=0.8;
h=85;
w=0.0025298;
gamma=1e6;

%J&C 15 constants
Q=0.9;
H=25;
w=0.0025298;
gamma=1e6;

prey=0:0.005:0.3;
temp=15;


%% Max ingestion
%Kiorboe & Hirst
ImaxS = exp(0.063*(temp-15.0)) .* 2.5 .* M_s^(-0.51) .* 24e-3;
ImaxM = exp(0.063*(temp-15.0)) .* 2.5 .* M_m^(-0.51) .* 24e-3;
ImaxL = exp(0.063*(temp-15.0)) .* 2.5 .* M_l^(-0.51) .* 24e-3;
%Hartvig
cmaxS = exp(0.063*(temp-10.0)) .* (h.*M_s^-0.25)/365;
cmaxM = exp(0.063*(temp-10.0)) .* (h.*M_m^-0.25)/365;
cmaxL = exp(0.063*(temp-10.0)) .* (h.*M_l^-0.25)/365;

%% Respiration
%Kiorboe & Hirst
% %ASSUME TABLE UNITS (mL O2/mg C/hr) ARE WRONG AND FIGURE UNITS (uL O2/mg C/hr) ARE CORRECT
% % *0.0252 -> (g/g/day)
% respS = exp(0.063*(temp-15.0)) .* 9.12 .* M_s^(-0.22) .* 0.0252; %2.8027;
% respM = exp(0.063*(temp-15.0)) .* 9.12 .* M_m^(-0.22) .* 0.0252; %2.8027;
% respL = exp(0.063*(temp-15.0)) .* 9.12 .* M_l^(-0.22) .* 0.0252; %2.8027;
%ASSUME TABLE UNITS (mL O2/mg C/hr) ARE CORRECT and conversions from Charlie
respS = exp(0.063*(temp-15.0)) .* 9.12 .* M_s^(-0.22) .* 9.12e-4; 
respM = exp(0.063*(temp-15.0)) .* 9.12 .* M_m^(-0.22) .* 9.12e-4; 
respL = exp(0.063*(temp-15.0)) .* 9.12 .* M_l^(-0.22) .* 9.12e-4;
%Hartvig
metS = 0.1*cmaxS;
metM = 0.1*cmaxM;
metL = 0.1*cmaxL;

%% Max clearance rate/search volume
%Kiorboe & Hirst
FmaxS = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_s^(-0.24) .* (24e-3/9);
FmaxM = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_m^(-0.24) .* (24e-3/9);
FmaxL = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_l^(-0.24) .* (24e-3/9);
%Hartvig
AS = exp(0.063*(temp-10.0)) .* (1e6*M_s^-0.2)/365;
AM = exp(0.063*(temp-10.0)) .* (1e6*M_m^-0.2)/365;
AL = exp(0.063*(temp-10.0)) .* (1e6*M_l^-0.2)/365;
%related to resting metab
svS = 4*respS;
svM = 4*respM;
svL = 4*respL;

%%
imax=[ImaxS;ImaxM;ImaxL];
cmax=[cmaxS;cmaxM;cmaxL];
fmax=[FmaxS;FmaxM;FmaxL];
amax=[AS;AM;AL];
resp=[respS;respM;respL];
met=[metS;metM;metL];
sv=[svS;svM;svL];

% figure(1)
% %subplot(2,2,1)
% plot(log10(m),log10(amax),'.r','MarkerSize',25); hold on;
% plot(log10(m),log10(fmax),'.b','MarkerSize',25); hold on;
% plot(log10(m),log10(sv),'.k','MarkerSize',25); hold on;
% ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
% xlabel('log10 weight (gWW)')
% legend('Hart','K&H','4met')
% legend('location','southwest')
% title('Specific clearance rate')
% print('-dpng','Clearance_comp_mass.png')
% 
% figure(2)
% %subplot(2,2,1)
% plot(log10(m),log10(cmax),'.r','MarkerSize',25); hold on;
% plot(log10(m),log10(imax),'.b','MarkerSize',25); hold on;
% ylabel('log10 ingestion (gWW gWW^-^1 d^-^1)')
% xlabel('log10 weight (gWW)')
% legend('Hart','K&H')
% legend('location','southwest')
% title('Specific ingestion rate')
% print('-dpng','Ingestion_comp_mass.png')

figure(3)
%subplot(2,2,1)
plot(log10(m),log10(met),'.r','MarkerSize',25); hold on;
plot(log10(m),log10(resp),'.b','MarkerSize',25); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('Hart','K&H')
legend('location','southwest')
title('Specific respiration rate')
print('-dpng','Respiration_comp_mass.png')

figure(4)
%subplot(2,2,1)
plot(log10(m),log10(resp),'.b','MarkerSize',25); hold on;
plot(log10(m),log10(0.75*resp),'.g','MarkerSize',25); hold on;
plot(log10(m),log10(0.5*resp),'.y','MarkerSize',25); hold on;
plot(log10(m),log10(0.25*resp),'.','Color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
plot(log10(m),log10(0.5*met),'.r','MarkerSize',25); hold on;
plot(log10(m),log10(met),'.k','MarkerSize',25); hold on;
plot(log10(m),log10(2*met),'.c','MarkerSize',25); hold on;
plot(log10(m),log10(4*met),'.m','MarkerSize',25); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('K&H','75K&H','50K&H','25K&H','fcrit05','fcrit10','fcrit20','fcrit40')
legend('location','southwest')
title('Specific respiration rate')

%% K&H metab vs. Imax
figure(5)
%subplot(2,2,1)
plot(log10(m),log10(resp),'.b','MarkerSize',25); hold on;
plot(log10(m),log10(imax),'sb','MarkerSize',10,'MarkerFaceColor','b'); hold on;
plot(log10(m),log10(met),'.r','MarkerSize',25); hold on;
plot(log10(m),log10(cmax),'sr','MarkerSize',10,'MarkerFaceColor','r'); hold on;
ylabel('log10 respiration or ingestion (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('K&H R','K&H I','Hart R','Hart I')
legend('location','southwest')
title('Specific respiration rate')
print('-dpng','KH_resp_vs_ingest_mass.png')

%% K&H metab vs. Imax using slope of -0.25

ing = exp(0.063*(temp-15.0)) .* 2.52 .* m.^(-0.25) .* 24e-3;
ingA = exp(0.063*(temp-15.0)) .* 2.52 .* m.^(-0.18) .* 24e-3;

figure(5)
%subplot(2,2,1)
plot(log10(m),log10(resp),'.b','MarkerSize',25); hold on;
plot(log10(m),log10(0.75*resp),'.g','MarkerSize',25); hold on;
plot(log10(m),log10(0.5*resp),'.y','MarkerSize',25); hold on;
plot(log10(m),log10(0.25*resp),'.','Color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
plot(log10(m),log10(ing),'.r','MarkerSize',25); hold on;
plot(log10(m),log10(ingA),'.k','MarkerSize',25); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('K&H','75K&H','50K&H','25K&H','Imax slope-.25','Imax slope-.18')
legend('location','southwest')
title('Specific respiration rate')
%print('-dpng','KH_resp_vs_ingest_mass.png')


%% vs temp
temp=-2:30;

%Max ingestion
ImaxS = exp(0.063*(temp-15.0)) .* 2.5 .* M_s^(-0.51) .* 24e-3;
ImaxM = exp(0.063*(temp-15.0)) .* 2.5 .* M_m^(-0.51) .* 24e-3;
ImaxL = exp(0.063*(temp-15.0)) .* 2.5 .* M_l^(-0.51) .* 24e-3;

cmaxS = exp(0.063*(temp-10.0)) .* (h.*M_s^-0.25)/365;
cmaxM = exp(0.063*(temp-10.0)) .* (h.*M_m^-0.25)/365;
cmaxL = exp(0.063*(temp-10.0)) .* (h.*M_l^-0.25)/365;

%Max clearance rate/search volume
FmaxLarv = exp(0.063*(temp-15.0)) .* 10^(8.1) .* M_s^(0.057) .* (24e-3/9);
FmaxS = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_s^(-0.24) .* (24e-3/9);
FmaxM = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_m^(-0.24) .* (24e-3/9);
FmaxL = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_l^(-0.24) .* (24e-3/9);

AS = exp(0.063*(temp-10.0)) .* (1e6*M_s^-0.2)/365;
AM = exp(0.063*(temp-10.0)) .* (1e6*M_m^-0.2)/365;
AL = exp(0.063*(temp-10.0)) .* (1e6*M_l^-0.2)/365;

%Respiration
respS = exp(0.063*(temp-15.0)) .* 9.12 .* M_s^(-0.22) .* 2.8027;
respM = exp(0.063*(temp-15.0)) .* 9.12 .* M_m^(-0.22) .* 2.8027;
respL = exp(0.063*(temp-15.0)) .* 9.12 .* M_l^(-0.22) .* 2.8027;

metS = 0.05*cmaxS;
metM = 0.05*cmaxM;
metL = 0.05*cmaxL;

%%

figure
subplot(2,2,1)
plot(temp,log10(AS),'r','LineWidth',2); hold on;
plot(temp,log10(AM),'k','LineWidth',2); hold on;
plot(temp,log10(AL),'b','LineWidth',2); hold on;
plot(temp,log10(FmaxS),'r--','LineWidth',2); hold on;
plot(temp,log10(FmaxM),'k--','LineWidth',2); hold on;
plot(temp,log10(FmaxL),'b--','LineWidth',2); hold on;
plot(temp,log10(FmaxLarv),'m--','LineWidth',2); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('temp (C)')
legend('Our S','Our M','Our L','K&H S','K&H M','K&H L','K&H larv')
legend('location','west')
title('Specific clearance rate')
%print('-dpng','Clearance_comp_temp.png')

%
figure
plot(temp,log10(cmaxS),'r','LineWidth',2); hold on;
plot(temp,log10(cmaxM),'k','LineWidth',2); hold on;
plot(temp,log10(cmaxL),'b','LineWidth',2); hold on;
plot(temp,log10(ImaxS),'r--','LineWidth',2); hold on;
plot(temp,log10(ImaxM),'k--','LineWidth',2); hold on;
plot(temp,log10(ImaxL),'b--','LineWidth',2); hold on;
ylabel('log10 ingestion (gWW gWW^-^1 d^-^1)')
xlabel('temp (C)')
legend('Our S','Our M','Our L','K&H S','K&H M','K&H L')
legend('location','west')
title('Specific ingestion rate')
%print('-dpng','Ingestion_comp_temp.png')

figure
plot(temp,log10(metS),'r','LineWidth',2); hold on;
plot(temp,log10(metM),'k','LineWidth',2); hold on;
plot(temp,log10(metL),'b','LineWidth',2); hold on;
plot(temp,log10(respS),'r--','LineWidth',2); hold on;
plot(temp,log10(respM),'k--','LineWidth',2); hold on;
plot(temp,log10(respL),'b--','LineWidth',2); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('temp (C)')
legend('Our S','Our M','Our L','K&H S','K&H M','K&H L')
legend('location','west')
title('Specific respiration rate')
%print('-dpng','Respiration_comp_temp.png')

%% K&H data
load('K&H14_clearance_rates.mat')
massC=mass;
load('K&H14_feeding_rates.mat')

W=[massC;mass];
F=[clearance;Fmax];
sF=[clearance./massC;SpecFmax];

figure
subplot(2,2,1)
plot(log10(massC),log10(clearance),'o')

figure
subplot(2,2,1)
plot(log10(massC),log10(clearance./massC),'o')

figure
subplot(2,2,1)
plot(log10(mass),log10(Fmax),'o')

figure
subplot(2,2,1)
plot(log10(mass),log10(SpecFmax),'o')

figure
subplot(2,2,1)
plot(log10(mass),log10(Imax),'o')

figure
subplot(2,2,1)
plot(log10(mass),log10(SpecImax),'o')

% remove outlier
nn=find(~isnan(Imax));
I2=Imax(nn);
sI2=SpecImax(nn);
m2=mass(nn);

figure
subplot(2,2,1)
plot(log10(m2(2:end)),log10(I2(2:end)),'o')

figure
subplot(2,2,1)
plot(log10(m2(2:end)),log10(sI2(2:end)),'o')


