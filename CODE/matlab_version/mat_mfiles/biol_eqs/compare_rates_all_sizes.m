% Compare diff Cmax, metab, A, ingest eqs
% Old sizes, big sizes, new sizes

clear all
close all

cm21=[1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    0 1 1;...     %c
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.75 0.75 0.75;... %lt grey
    0.5 0.5 0.5;...    %med grey
    49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    1 1 0;...      %yellow
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255]; %peach

set(groot,'defaultAxesColorOrder',cm21);

%% OLD
L_s = 10^((log10(2)+log10(20))/2);
L_m = 10^((log10(20)+log10(200))/2);
L_l = 10^((log10(200)+log10(2000))/2);

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;
Om=[M_s; M_m; M_l];

%% BIG
L_s = 10.0; % small
L_m = 200.0; % medium
L_l = 1.0e3;% large

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;
Bm=[M_s; M_m; M_l];

%% NEW
%!Individual Mass (g)
M_s = 10^((log10(0.001)+log10(0.5))/2);
M_m = 10^((log10(0.5)+log10(250))/2);
M_l = 10^((log10(250)+log10(125000))/2);

%! Body lengths (mm)
% Convert from mm to cm and use their const coeff = 0.01g/cm3
L_s = 10.0 * (M_s/0.01)^(1/3); % small
L_m = 10.0 * (M_m/0.01)^(1/3); % medium
L_l = 10.0 * (M_l/0.01)^(1/3); % large

Nm=[M_s; M_m; M_l];

M = [Om, Bm, Nm];

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
Hu = exp(0.063*(temp-10.0)) .* (u.*M.^(-n))/365;

%mizer
Zu = exp(0.063*(temp-10.0)) .* (u2.*M.^(-n))/365;

%J&C
Ju = exp(0.063*(temp-10.0)) .* (U.*M.^(-N))/365;

%% Annual Mortality
%Hartvig
HuA = exp(0.063*(temp-10.0)) .* (u.*M.^(-n));

%mizer
ZuA = exp(0.063*(temp-10.0)) .* (u2.*M.^(-n));

%J&C
JuA = exp(0.063*(temp-10.0)) .* (U.*M.^(-N));


%% Max ingestion
%Kiorboe & Hirst
Imax = exp(0.063*(temp-15.0)) .* 2.5 .* M.^(-0.51) .* 24e-3;

%Hartvig
Hcmax = exp(0.063*(temp-10.0)) .* (h.*M.^-n)/365;

%mizer
Zcmax = exp(0.063*(temp-10.0)) .* (h2.*M.^-n)/365;

%J&C
Jcmax = exp(0.063*(temp-10.0)) .* (H.*M.^-N)/365;

%Me
Mcmax = exp(0.063*(temp-15.0)) .* (60.*M.^-n)/365;

%% Respiration
%Hartvig
Hmet = exp(0.063*(temp-10.0)) .* k.*M.^-p;

%mizer
Zmet = exp(0.063*(temp-10.0)) .* k2.*M.^-p;

%J&C
Jmet = exp(0.063*(temp-10.0)) .* K.*M.^-p;

%Me
Mmet = exp(0.063*(temp-15.0)) .* 0.4*Mcmax;

%% Max clearance rate/search volume
%Kiorboe & Hirst = me
Fmax = exp(0.063*(temp-15.0)) .* 10^(3.24) .* M.^(-0.24) .* (24e-3/9);

%Hartvig
HA = exp(0.063*(temp-10.0)) .* (gamma*M.^-q)/365;

%mizer
ZA = exp(0.063*(temp-10.0)) .* (gamma2*M.^-q)/365;

%J&C
JA = exp(0.063*(temp-10.0)) .* (G*M.^-Q)/365;

%%
figure(50)
bar(log10(M))
xlabel('Size class')
ylabel('log10 weight (gWW)')
legend('Old','Big','New')
legend('location','northwest')
title('Mass')
print('-dpng','Comp_mass_all_sizes.png')

%%
figure(1)
%subplot(2,2,1)
plot(log10(M),log10(Fmax),'.','MarkerSize',25); hold on;
plot(log10(M),log10(HA),'o','MarkerSize',15); hold on;
plot(log10(M),log10(JA),'s','MarkerSize',15); hold on;
plot(log10(M),log10(ZA),'^','MarkerSize',15); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('oK&H','bK&H','nK&H','oHart','bHart','nHart','oJ&C','bJ&C','nJ&C','omizer','bmizer','nmizer')
legend('location','southwest')
title('Specific clearance rate')
print('-dpng','Clearance_comp_mass_all_sizes.png')

%%
figure(2)
%subplot(2,2,1)
plot(log10(M),log10(Imax),'.','MarkerSize',25); hold on;
plot(log10(M),log10(Hcmax),'o','MarkerSize',15); hold on;
plot(log10(M),log10(Jcmax),'s','MarkerSize',15); hold on;
plot(log10(M),log10(Zcmax),'^','MarkerSize',15); hold on;
plot(log10(M),log10(Mcmax),'*','MarkerSize',15); hold on;
ylabel('log10 ingestion (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('oK&H','bK&H','nK&H','oHart','bHart','nHart','oJ&C','bJ&C','nJ&C','omizer','bmizer','nmizer','oMe','bMe','nMe')
legend('location','southwest')
title('Specific ingestion rate')
print('-dpng','Ingestion_comp_mass_all_sizes.png')

%%
figure(3)
%subplot(2,2,1)
plot(log10(M),log10(Hmet),'o','MarkerSize',15); hold on;
plot(log10(M),log10(Jmet),'s','MarkerSize',15); hold on;
plot(log10(M),log10(Zmet),'^','MarkerSize',15); hold on;
plot(log10(M),log10(Mmet),'*','MarkerSize',15); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('oHart','bHart','nHart','oJ&C','bJ&C','nJ&C','omizer','bmizer','nmizer','oMe','bMe','nMe')
legend('location','eastoutside')
title('Specific respiration rate')
print('-dpng','Respiration_comp_mass_all_sizes.png')

%%
figure(4)
%subplot(2,2,1)
plot(log10(M),log10(Hmet),'o','MarkerSize',15); hold on;
plot(log10(M),log10(Jmet),'s','MarkerSize',15); hold on;
plot(log10(M),log10(Zmet),'^','MarkerSize',15); hold on;
plot(log10(M),log10(Mmet),'*','MarkerSize',15); hold on;
plot(log10(M),log10(0.2*Mcmax),'.','MarkerSize',25); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('oHart','bHart','nHart','oJ&C','bJ&C','nJ&C','omizer','bmizer','nmizer','oMe4','bMe4','nMe4','oMe2','bMe2','nMe2')
legend('location','eastoutside')
title('Specific respiration rate')
print('-dpng','Respiration_comp_fcrit_all_sizes.png')

%%
figure(5)
%subplot(2,2,1)
plot(log10(M),(Hu),'o','MarkerSize',15); hold on;
plot(log10(M),(Ju),'s','MarkerSize',15); hold on;
plot(log10(M),(Zu),'^','MarkerSize',15); hold on;
ylabel('mortality rate (d^-^1)')
xlabel('log10 weight (gWW)')
legend('oHart','bHart','nHart','oJ&C','bJ&C','nJ&C','omizer','bmizer','nmizer')
legend('location','northeast')
title('Daily mortality rate')
print('-dpng','Mortality_comp_mass_all_sizes.png')


