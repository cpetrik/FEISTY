% Compare search rates (A) from
% Hartvig, J&C, mizer, Blanchard

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
q=1-0.8;
h=85;
w=0.0025298;
gamma=0.8e4;
n=1-0.75;
p=1-0.75;
k=10;
u=0.84;

%mizer constants
h2=40;
gamma2=2.9e3;
k2=4.8;
u2=3;

%J&C 15 constants
Q=1-0.9;
H=25;
W=0.0025298;
G=1.37e4;
N=1-0.67;
K=2.5;
U=0.5;

prey=0:0.005:0.3;
temp=15;

%% Mortality
%Hartvig
Hu = exp(0.063*(temp-10.0)) .* (u.*m.^(-n))/365;

%mizer
Zu = exp(0.063*(temp-10.0)) .* (u2.*m.^(-n))/365;

%J&C
Ju = exp(0.063*(temp-10.0)) .* (U.*m.^(-N))/365;


%% Max clearance rate/search volume
%Kiorboe & Hirst = me
Fmax = exp(0.063*(temp-15.0)) .* 10^(3.24) .* m.^(-0.24) .* (24e-3/9);

%Hartvig
HA = exp(0.063*(temp-10.0)) .* (gamma*m.^-q)/365;

%mizer
ZA = exp(0.063*(temp-10.0)) .* (gamma2*m.^-q)/365;

%J&C
JA = exp(0.063*(temp-10.0)) .* (G*m.^-Q)/365;

%Blanchard pelagic
BAp = (exp(0.063*(temp-15.0)) * 640 * m.^(-0.18)) ./365.0;

%Blanchard benthic
BAd = (exp(0.063*(temp-15.0)) * 64 * m.^(-0.25)) ./365.0;

%%
figure(1)
%subplot(2,2,1)
plot(log10(m),log10(Fmax),'.-r','MarkerSize',25); hold on;
plot(log10(m),log10(HA),'.-b','MarkerSize',25); hold on;
plot(log10(m),log10(JA),'.-k','MarkerSize',25); hold on;
plot(log10(m),log10(ZA),'.-','color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
plot(log10(m),log10(BAp),'.-m','MarkerSize',25); hold on;
plot(log10(m),log10(BAd),'.-c','MarkerSize',25); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('K&H-me','Hart','J&C','mizer','Bpel','Bben')
legend('location','southwest')
title('Specific clearance rate')
print('-dpng','Clearance_comp_mass_Blanchard.png')


%% vs temp
temp=-2:30;

%Max clearance rate/search volume
FmaxS = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_s^(-0.24) .* (24e-3/9);
FmaxM = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_m^(-0.24) .* (24e-3/9);
FmaxL = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_l^(-0.24) .* (24e-3/9);

HAS = exp(0.063*(temp-10.0)) .* (1e6*M_s^-q)/365;
HAM = exp(0.063*(temp-10.0)) .* (1e6*M_m^-q)/365;
HAL = exp(0.063*(temp-10.0)) .* (1e6*M_l^-q)/365;

ZAS = exp(0.063*(temp-10.0)) .* (gamma2*M_s.^-q)/365;
ZAM = exp(0.063*(temp-10.0)) .* (gamma2*M_m.^-q)/365;
ZAL = exp(0.063*(temp-10.0)) .* (gamma2*M_l.^-q)/365;

JAS = exp(0.063*(temp-10.0)) .* (1e6*M_s^-Q)/365;
JAM = exp(0.063*(temp-10.0)) .* (1e6*M_m^-Q)/365;
JAL = exp(0.063*(temp-10.0)) .* (1e6*M_l^-Q)/365;

%Blanchard pelagic
BApS = (exp(0.063*(temp-15.0)) .* 640 .* M_s.^(-0.18)) ./365.0;
BApM = (exp(0.063*(temp-15.0)) .* 640 .* M_m.^(-0.18)) ./365.0;
BApL = (exp(0.063*(temp-15.0)) .* 640 .* M_l.^(-0.18)) ./365.0;

%Blanchard benthic
BAdS = (exp(0.063*(temp-15.0)) .* 64 .* M_s.^(-0.25)) ./365.0;
BAdM = (exp(0.063*(temp-15.0)) .* 64 .* M_m.^(-0.25)) ./365.0;
BAdL = (exp(0.063*(temp-15.0)) .* 64 .* M_l.^(-0.25)) ./365.0;


%%

figure(2)
%subplot(2,2,1)
plot(temp,log10(HAS),'r','LineWidth',2); hold on;
plot(temp,log10(HAM),'k','LineWidth',2); hold on;
plot(temp,log10(HAL),'b','LineWidth',2); hold on;
plot(temp,log10(JAS),'--r','LineWidth',2); hold on;
plot(temp,log10(JAM),'--k','LineWidth',2); hold on;
plot(temp,log10(JAL),'--b','LineWidth',2); hold on;
plot(temp,log10(ZAS),'vr','LineWidth',1); hold on;
plot(temp,log10(ZAM),'vk','LineWidth',1); hold on;
plot(temp,log10(ZAL),'vb','LineWidth',1); hold on;
plot(temp,log10(FmaxS),':r','LineWidth',2); hold on;
plot(temp,log10(FmaxM),':k','LineWidth',2); hold on;
plot(temp,log10(FmaxL),':b','LineWidth',2); hold on;
plot(temp,log10(BApS),'-.r','LineWidth',2); hold on;
plot(temp,log10(BApM),'-.k','LineWidth',2); hold on;
plot(temp,log10(BApL),'-.b','LineWidth',2); hold on;
plot(temp,log10(BAdS),'or','LineWidth',1); hold on;
plot(temp,log10(BAdM),'ok','LineWidth',1); hold on;
plot(temp,log10(BAdL),'ob','LineWidth',1); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('temp (C)')
legend('HS','HM','HL','JS','JM','JL','MS','MM','ML','meS','meM','meL','BpS',...
    'BpM','BpL','BdS','BdM','BdL')
legend('location','eastoutside')
title('Specific clearance rate')
print('-dpng','Clearance_comp_temp_Blanchard.png')

