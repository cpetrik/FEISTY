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
plot(log10(m),log10(Fmax),'r','MarkerSize',25); hold on;
plot(log10(m),log10(HA),'b','MarkerSize',25); hold on;
plot(log10(m),log10(JA),'k','MarkerSize',25); hold on;
plot(log10(m),log10(ZA),'color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('K&H-me','Hart','J&C','mizer')
legend('location','southwest')
title('Specific clearance rate')
print('-dpng','Clearance_comp_mass_new_sizes.png')

figure(2)
%subplot(2,2,1)
plot(log10(m),log10(Imax),'r','MarkerSize',25); hold on;
plot(log10(m),log10(Hcmax),'b','MarkerSize',25); hold on;
plot(log10(m),log10(Jcmax),'k','MarkerSize',25); hold on;
plot(log10(m),log10(Zcmax),'color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
plot(log10(m),log10(Mcmax),'ob','MarkerSize',10); hold on;
ylabel('log10 ingestion (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('K&H','Hart','J&C','mizer','me')
legend('location','southwest')
title('Specific ingestion rate')
print('-dpng','Ingestion_comp_mass_new_sizes.png')

figure(3)
%subplot(2,2,1)
plot(log10(m),log10(Hmet),'b','MarkerSize',25); hold on;
plot(log10(m),log10(Jmet),'k','MarkerSize',25); hold on;
plot(log10(m),log10(Zmet),'color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
plot(log10(m),log10(Mmet),'ob','MarkerSize',10); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('Hart','J&C','mizer','me')
legend('location','southwest')
title('Specific respiration rate')
print('-dpng','Respiration_comp_mass_new_sizes.png')

figure(4)
%subplot(2,2,1)
plot(log10(m),log10(Hmet),'b','MarkerSize',25); hold on;
plot(log10(m),log10(Jmet),'k','MarkerSize',25); hold on;
plot(log10(m),log10(Zmet),'color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
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
plot(log10(m),(Hu),'b','MarkerSize',25); hold on;
plot(log10(m),(Ju),'k','MarkerSize',25); hold on;
plot(log10(m),(Zu),'color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
ylabel('mortality rate (d^-^1)')
xlabel('log10 weight (gWW)')
legend('Hart','J&C','mizer')
legend('location','southwest')
title('Daily mortality rate')
print('-dpng','Mortality_comp_mass_new_sizes.png')


%% theoretical
% Respiration
met1 = 0.2 * exp(0.063*(temp-10.0)) .* 10.*m.^-p;
met2 = 0.2 * exp(0.063*(temp-10.0)) .* 20.*m.^-p;
met3 = 0.2 * exp(0.063*(temp-10.0)) .* 30.*m.^-p;
met4 = 0.2 * exp(0.063*(temp-10.0)) .* 40.*m.^-p;
met5 = 0.2 * exp(0.063*(temp-10.0)) .* 25.*m.^-0.33;

% Max clearance rate/search volume
A1 = exp(0.063*(temp-10.0)) .* (gamma*m.^-0.25)/365;
A2 = exp(0.063*(temp-10.0)) .* (gamma*m.^-0.2)/365;
A3 = exp(0.063*(temp-10.0)) .* (gamma*m.^-0.1)/365;

figure(6)
%subplot(2,2,1)
plot(log10(m),log10(A1),'r','MarkerSize',25); hold on;
plot(log10(m),log10(A2),'b','MarkerSize',25); hold on;
plot(log10(m),log10(A3),'k','MarkerSize',25); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('q=0.25','q=0.2','q=0.1')
legend('location','southwest')
title('Specific clearance rate')
print('-dpng','Clearance_theor_mass_new_sizes.png')

figure(7)
%subplot(2,2,1)
plot(log10(m),log10(met1),'r','MarkerSize',25); hold on;
plot(log10(m),log10(met2),'b','MarkerSize',25); hold on;
plot(log10(m),log10(met3),'color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
plot(log10(m),log10(met4),'k','MarkerSize',10); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('h=10','h=20','h=30','h=40')
legend('location','southwest')
title('Specific respiration rate')
print('-dpng','Respiration_theor1_mass_new_sizes.png')

figure(8)
%subplot(2,2,1)
plot(log10(m),log10(met2),'b','MarkerSize',25); hold on;
plot(log10(m),log10(met5),'m','MarkerSize',10); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('h=20,p=0.25','h=25,p=0.33')
legend('location','southwest')
title('Specific respiration rate')
print('-dpng','Respiration_theor2_mass_new_sizes.png')

