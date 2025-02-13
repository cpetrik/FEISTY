% Compare diff Cmax, metab, A, ingest eqs

clear all
close all

%!Individual Mass (g)
M_s = 10^((log10(0.001)+log10(0.5))/2);
M_m = 10^((log10(0.5)+log10(250))/2);
M_l = 10^((log10(250)+log10(125000))/2);

%! Body lengths (mm)
% Convert from mm to cm and use their const coeff = 0.01g/cm3
L_s = 10.0 * (M_s/0.01)^(1/3); % small
L_m = 10.0 * (M_m/0.01)^(1/3); % medium
L_l = 10.0 * (M_l/0.01)^(1/3); % large

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

%% Annual Mortality
%Hartvig
HuA = exp(0.063*(temp-10.0)) .* (u.*m.^(-n));

%mizer
ZuA = exp(0.063*(temp-10.0)) .* (u2.*m.^(-n));

%J&C
JuA = exp(0.063*(temp-10.0)) .* (U.*m.^(-N));


%% Max ingestion
%Kiorboe & Hirst
Imax = exp(0.063*(temp-15.0)) .* 2.5 .* m.^(-0.51) .* 24e-3;

%Hartvig
Hcmax = exp(0.063*(temp-10.0)) .* (h.*m.^-n)/365;

%mizer
Zcmax = exp(0.063*(temp-10.0)) .* (h2.*m.^-n)/365;

%J&C
Jcmax = exp(0.063*(temp-10.0)) .* (H.*m.^-N)/365;

%Me
Mcmax = exp(0.063*(temp-15.0)) .* (60.*m.^-n)/365;

%% Respiration
%Hartvig
Hmet = exp(0.063*(temp-10.0)) .* k.*m.^-p;

%mizer
Zmet = exp(0.063*(temp-10.0)) .* k2.*m.^-p;

%J&C
Jmet = exp(0.063*(temp-10.0)) .* K.*m.^-p;

%Me
Mmet = exp(0.063*(temp-15.0)) .* 0.4*Mcmax;

%% Max clearance rate/search volume
%Kiorboe & Hirst = me
Fmax = exp(0.063*(temp-15.0)) .* 10^(3.24) .* m.^(-0.24) .* (24e-3/9);

%Hartvig
HA = exp(0.063*(temp-10.0)) .* (gamma*m.^-q)/365;

%mizer
ZA = exp(0.063*(temp-10.0)) .* (gamma2*m.^-q)/365;

%J&C
JA = exp(0.063*(temp-10.0)) .* (G*m.^-Q)/365;

%%
figure(1)
%subplot(2,2,1)
plot(log10(m),log10(Fmax),'.r','MarkerSize',25); hold on;
plot(log10(m),log10(HA),'.b','MarkerSize',25); hold on;
plot(log10(m),log10(JA),'.k','MarkerSize',25); hold on;
plot(log10(m),log10(ZA),'.','color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('K&H-me','Hart','J&C','mizer')
legend('location','southwest')
title('Specific clearance rate')
print('-dpng','Clearance_comp_mass_new_sizes.png')

figure(2)
%subplot(2,2,1)
plot(log10(m),log10(Imax),'.r','MarkerSize',25); hold on;
plot(log10(m),log10(Hcmax),'.b','MarkerSize',25); hold on;
plot(log10(m),log10(Jcmax),'.k','MarkerSize',25); hold on;
plot(log10(m),log10(Zcmax),'.','color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
plot(log10(m),log10(Mcmax),'ob','MarkerSize',10); hold on;
ylabel('log10 ingestion (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('K&H','Hart','J&C','mizer','me')
legend('location','southwest')
title('Specific ingestion rate')
print('-dpng','Ingestion_comp_mass_new_sizes.png')

figure(3)
%subplot(2,2,1)
plot(log10(m),log10(Hmet),'.b','MarkerSize',25); hold on;
plot(log10(m),log10(Jmet),'.k','MarkerSize',25); hold on;
plot(log10(m),log10(Zmet),'.','color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
plot(log10(m),log10(Mmet),'ob','MarkerSize',10); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('Hart','J&C','mizer','me')
legend('location','southwest')
title('Specific respiration rate')
print('-dpng','Respiration_comp_mass_new_sizes.png')

figure(4)
%subplot(2,2,1)
plot(log10(m),log10(Hmet),'.b','MarkerSize',25); hold on;
plot(log10(m),log10(Jmet),'.k','MarkerSize',25); hold on;
plot(log10(m),log10(Zmet),'.','color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
plot(log10(m),log10(Mmet),'ob','MarkerSize',10); hold on;
plot(log10(m),log10(0.1*Mcmax),'.r','MarkerSize',25); hold on;
plot(log10(m),log10(0.2*Mcmax),'.g','MarkerSize',25); hold on;
plot(log10(m),log10(0.3*Mcmax),'.c','MarkerSize',25); hold on;
plot(log10(m),log10(0.5*Mcmax),'.m','MarkerSize',25); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('Hart','J&C','mizer','fcrit40','fcrit10','fcrit20','fcrit30','fcrit50')
legend('location','southwest')
title('Specific respiration rate')
print('-dpng','Respiration_comp_fcrit_new_sizes.png')

figure(5)
%subplot(2,2,1)
plot(log10(m),(Hu),'.b','MarkerSize',25); hold on;
plot(log10(m),(Ju),'.k','MarkerSize',25); hold on;
plot(log10(m),(Zu),'.','color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
ylabel('mortality rate (d^-^1)')
xlabel('log10 weight (gWW)')
legend('Hart','J&C','mizer')
legend('location','southwest')
title('Daily mortality rate')
print('-dpng','Mortality_comp_mass_new_sizes.png')


%% vs temp
temp=-2:30;

%Max ingestion
ImaxS = exp(0.063*(temp-15.0)) .* 2.5 .* M_s^(-0.51) .* 24e-3;
ImaxM = exp(0.063*(temp-15.0)) .* 2.5 .* M_m^(-0.51) .* 24e-3;
ImaxL = exp(0.063*(temp-15.0)) .* 2.5 .* M_l^(-0.51) .* 24e-3;

HcmaxS = exp(0.063*(temp-10.0)) .* (h.*M_s^-0.25)/365;
HcmaxM = exp(0.063*(temp-10.0)) .* (h.*M_m^-0.25)/365;
HcmaxL = exp(0.063*(temp-10.0)) .* (h.*M_l^-0.25)/365;

JcmaxS = exp(0.063*(temp-10.0)) .* (H.*M_s^-0.25)/365;
JcmaxM = exp(0.063*(temp-10.0)) .* (H.*M_m^-0.25)/365;
JcmaxL = exp(0.063*(temp-10.0)) .* (H.*M_l^-0.25)/365;

McmaxS = exp(0.063*(temp-15.0)) .* (60.*M_s^-0.25)/365;
McmaxM = exp(0.063*(temp-15.0)) .* (60.*M_m^-0.25)/365;
McmaxL = exp(0.063*(temp-15.0)) .* (60.*M_l^-0.25)/365;

%Max clearance rate/search volume
FmaxS = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_s^(-0.24) .* (24e-3/9);
FmaxM = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_m^(-0.24) .* (24e-3/9);
FmaxL = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_l^(-0.24) .* (24e-3/9);

HAS = exp(0.063*(temp-10.0)) .* (1e6*M_s^-q)/365;
HAM = exp(0.063*(temp-10.0)) .* (1e6*M_m^-q)/365;
HAL = exp(0.063*(temp-10.0)) .* (1e6*M_l^-q)/365;

JAS = exp(0.063*(temp-10.0)) .* (1e6*M_s^-Q)/365;
JAM = exp(0.063*(temp-10.0)) .* (1e6*M_m^-Q)/365;
JAL = exp(0.063*(temp-10.0)) .* (1e6*M_l^-Q)/365;

%Respiration
%Hartvig
HmetS = exp(0.063*(temp-10.0)) .* k.*M_s.^-p;
HmetM = exp(0.063*(temp-10.0)) .* k.*M_m.^-p;
HmetL = exp(0.063*(temp-10.0)) .* k.*M_l.^-p;

%J&C
JmetS = exp(0.063*(temp-10.0)) .* K.*M_s.^-p;
JmetM = exp(0.063*(temp-10.0)) .* K.*M_m.^-p;
JmetL = exp(0.063*(temp-10.0)) .* K.*M_l.^-p;

MmetS = 0.3*McmaxS;
MmetM = 0.3*McmaxM;
MmetL = 0.3*McmaxL;

%% Mortality
%Hartvig
HuS = exp(0.063*(temp-10.0)) .* (u.*M_s.^(-n))/365;
HuM = exp(0.063*(temp-10.0)) .* (u.*M_m.^(-n))/365;
HuL = exp(0.063*(temp-10.0)) .* (u.*M_l.^(-n))/365;

%mizer
ZuS = exp(0.063*(temp-10.0)) .* (u2.*M_s.^(-n))/365;
ZuM = exp(0.063*(temp-10.0)) .* (u2.*M_m.^(-n))/365;
ZuL = exp(0.063*(temp-10.0)) .* (u2.*M_l.^(-n))/365;

%J&C
JuS = exp(0.063*(temp-10.0)) .* (U.*M_s.^(-N))/365;
JuM = exp(0.063*(temp-10.0)) .* (U.*M_m.^(-N))/365;
JuL = exp(0.063*(temp-10.0)) .* (U.*M_l.^(-N))/365;

%%

figure
%subplot(2,2,1)
plot(temp,log10(HAS),'r','LineWidth',2); hold on;
plot(temp,log10(HAM),'k','LineWidth',2); hold on;
plot(temp,log10(HAL),'b','LineWidth',2); hold on;
plot(temp,log10(JAS),'--r','LineWidth',2); hold on;
plot(temp,log10(JAM),'--k','LineWidth',2); hold on;
plot(temp,log10(JAL),'--b','LineWidth',2); hold on;
plot(temp,log10(FmaxS),':r','LineWidth',2); hold on;
plot(temp,log10(FmaxM),':k','LineWidth',2); hold on;
plot(temp,log10(FmaxL),':b','LineWidth',2); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('temp (C)')
legend('HS','HM','HL','JS','JM','JL','meS','meM','meL')
legend('location','west')
title('Specific clearance rate')
print('-dpng','Clearance_comp_temp_new_sizes.png')

%
figure
plot(temp,log10(HcmaxS),'r','LineWidth',2); hold on;
plot(temp,log10(HcmaxM),'k','LineWidth',2); hold on;
plot(temp,log10(HcmaxL),'b','LineWidth',2); hold on;
plot(temp,log10(JcmaxS),'--r','LineWidth',2); hold on;
plot(temp,log10(JcmaxM),'--k','LineWidth',2); hold on;
plot(temp,log10(JcmaxL),'--b','LineWidth',2); hold on;
plot(temp,log10(McmaxS),':r','LineWidth',2); hold on;
plot(temp,log10(McmaxM),':k','LineWidth',2); hold on;
plot(temp,log10(McmaxL),':b','LineWidth',2); hold on;
ylabel('log10 ingestion (gWW gWW^-^1 d^-^1)')
xlabel('temp (C)')
legend('HS','HM','HL','JS','JM','JL','meS','meM','meL')
legend('location','west')
title('Specific ingestion rate')
print('-dpng','Ingestion_comp_temp_new_sizes.png')

figure
plot(temp,log10(HmetS),'r','LineWidth',2); hold on;
plot(temp,log10(HmetM),'k','LineWidth',2); hold on;
plot(temp,log10(HmetL),'b','LineWidth',2); hold on;
plot(temp,log10(JmetS),'--r','LineWidth',2); hold on;
plot(temp,log10(JmetM),'--k','LineWidth',2); hold on;
plot(temp,log10(JmetL),'--b','LineWidth',2); hold on;
plot(temp,log10(MmetS),':r','LineWidth',2); hold on;
plot(temp,log10(MmetM),':k','LineWidth',2); hold on;
plot(temp,log10(MmetL),':b','LineWidth',2); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('temp (C)')
legend('HS','HM','HL','JS','JM','JL','meS','meM','meL')
legend('location','west')
title('Specific respiration rate')
print('-dpng','Respiration_comp_temp_new_sizes.png')

%%
figure
subplot(2,2,1)
plot(temp,(HuS),'r','LineWidth',2); hold on;
plot(temp,(JuS),'--r','LineWidth',2); hold on;
plot(temp,(ZuS),':r','LineWidth',2); hold on;
ylabel('mortality rate (d^-^1)')
title('S')
xlim([-2 30])
subplot(2,2,2)
plot(temp,(HuM),'k','LineWidth',2); hold on;
plot(temp,(JuM),'--k','LineWidth',2); hold on;
plot(temp,(ZuM),':k','LineWidth',2); hold on;
xlabel('temp (C)')
title('M')
xlim([-2 30])
subplot(2,2,3)
plot(temp,(HuL),'b','LineWidth',2); hold on;
plot(temp,(JuL),'--b','LineWidth',2); hold on;
plot(temp,(ZuL),':b','LineWidth',2); hold on;
xlabel('temp (C)')
ylabel('mortality rate (d^-^1)')
title('L')
xlim([-2 30])
legend('Hartvig11','J&C15','mizer')
legend('location','northwest')
print('-dpng','Mortality_comp_temp_new_sizes.png')
