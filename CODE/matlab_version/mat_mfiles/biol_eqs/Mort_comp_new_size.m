% Compare diff mortality rates from
% Hartvig, mizer, J&C15, Peterson&Wrob, McGurk

clear all
close all

M_s = 10^((log10(0.001)+log10(0.5))/2);  %0.0224
M_m = 10^((log10(0.5)+log10(250))/2);    %11.1803
M_l = 10^((log10(250)+log10(125000))/2); %5.5902e3

%! Body lengths (mm)
% Convert from mm to cm and use their const coeff = 0.01g/cm3
L_s = 10.0 * (M_s/0.01)^(1/3); % small   4.6-36.8 mm, 13.1 mm
L_m = 10.0 * (M_m/0.01)^(1/3); % medium 36.8-292 mm, 10.4 cm
L_l = 10.0 * (M_l/0.01)^(1/3); % large  0.292-2.32 m, 0.82 m

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

%P&W 84
Pu = exp(0.063*(temp-15.0)) .* 5.26e-3 .* (m/9).^(-0.25);

%McGurk 86
Mu = exp(0.063*(temp-15.0)) .* 2.2e-4 .* (m/9).^(-0.85);

%%
figure(3)
%subplot(2,2,1)
plot(log10(m),(Hu)*365,'.-b','MarkerSize',25); hold on;
plot(log10(m),(Ju)*365,'.-k','MarkerSize',25); hold on;
plot(log10(m),(Zu)*365,'.-','color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
plot(log10(m),(Pu)*365,'.-m','MarkerSize',25); hold on;
plot(log10(m),(Mu)*365,'.-r','MarkerSize',25); hold on;
ylabel('mortality rate (y^-^1)')
xlabel('log10 weight (gWW)')
legend('Hart','J&C','mizer','P&W','McG')
legend('location','northeast')
title('Annual mortality rate')
print('-dpng','Mortality_comp_mass_eggs_year.png')

figure(4)
%subplot(2,2,1)
plot(log10(m),(Hu),'.-b','MarkerSize',25); hold on;
plot(log10(m),(Ju),'.-k','MarkerSize',25); hold on;
plot(log10(m),(Zu),'.-','color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
plot(log10(m),(Pu),'.-m','MarkerSize',25); hold on;
plot(log10(m),(Mu),'.-r','MarkerSize',25); hold on;
ylabel('mortality rate (d^-^1)')
xlabel('log10 weight (gWW)')
legend('Hart','J&C','mizer','P&W','McG')
legend('location','northeast')
title('Daily mortality rate')
print('-dpng','Mortality_comp_mass_eggs.png')

figure(5)
%subplot(2,2,1)
plot(log10(m),(Hu),'.-b','MarkerSize',25); hold on;
plot(log10(m),(Ju),'.-k','MarkerSize',25); hold on;
plot(log10(m),(Zu),'.-','color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
plot(log10(m),(Pu),'.-m','MarkerSize',25); hold on;
plot(log10(m),(Mu),'.-r','MarkerSize',25); hold on;
ylim([0 0.06])
ylabel('mortality rate (d^-^1)')
xlabel('log10 weight (gWW)')
legend('Hart','J&C','mizer','P&W','McG')
legend('location','northeast')
title('Daily mortality rate')
print('-dpng','Mortality_comp_mass_ML.png')
%%
figure(6)
%subplot(2,2,1)
plot(log10(m),(Hu)*365,'.-b','MarkerSize',25); hold on;
plot(log10(m),(Ju)*365,'.-k','MarkerSize',25); hold on;
plot(log10(m),(Zu)*365,'.-','color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
plot(log10(m),(Pu)*365,'.-m','MarkerSize',25); hold on;
plot(log10(m),(Mu)*365,'.-r','MarkerSize',25); hold on;
ylim([0 4])
ylabel('mortality rate (d^-^1)')
xlabel('log10 weight (gWW)')
legend('Hart','J&C','mizer','P&W','McG')
legend('location','northeast')
title('Annual mortality rate')
print('-dpng','Mortality_comp_mass_ML_year.png')


%% vs temp
temp=-2:30;

% Mortality
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

%P&W 84
PuS = exp(0.063*(temp-15.0)) .* 5.26e-3 .* (M_s/9).^(-0.25);
PuM = exp(0.063*(temp-15.0)) .* 5.26e-3 .* (M_m/9).^(-0.25);
PuL = exp(0.063*(temp-15.0)) .* 5.26e-3 .* (M_l/9).^(-0.25);

%McGurk 86
MuS = exp(0.063*(temp-15.0)) .* 2.2e-4 .* (M_s/9).^(-0.85);
MuM = exp(0.063*(temp-15.0)) .* 2.2e-4 .* (M_m/9).^(-0.85);
MuL = exp(0.063*(temp-15.0)) .* 2.2e-4 .* (M_l/9).^(-0.85);


%%
figure(6)
subplot(2,2,1)
plot(temp,(HuS),'r','LineWidth',2); hold on;
plot(temp,(JuS),'--r','LineWidth',2); hold on;
plot(temp,(ZuS),':r','LineWidth',2); hold on;
plot(temp,(PuS),'-.r','LineWidth',2); hold on;
plot(temp,(MuS),'or','LineWidth',1); hold on;
ylabel('mortality rate (d^-^1)')
title('S')
xlim([-2 30])
subplot(2,2,2)
plot(temp,(HuM),'k','LineWidth',2); hold on;
plot(temp,(JuM),'--k','LineWidth',2); hold on;
plot(temp,(ZuM),':k','LineWidth',2); hold on;
plot(temp,(PuM),'-.k','LineWidth',2); hold on;
plot(temp,(MuM),'ok','LineWidth',1); hold on;
xlabel('temp (C)')
title('M')
xlim([-2 30])
subplot(2,2,3)
plot(temp,(HuL),'b','LineWidth',2); hold on;
plot(temp,(JuL),'--b','LineWidth',2); hold on;
plot(temp,(ZuL),':b','LineWidth',2); hold on;
plot(temp,(PuL),'-.b','LineWidth',2); hold on;
plot(temp,(MuL),'ob','LineWidth',1); hold on;
xlabel('temp (C)')
ylabel('mortality rate (d^-^1)')
title('L')
xlim([-2 30])
legend('Hartvig11','J&C15','mizer','P&W','McG')
legend('location','northwest')
print('-dpng','Mortality_comp_temp_eggs.png')
