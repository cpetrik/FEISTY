% Compare diff Cmax, metab, A, ingest eqs

clear 
close all

prey=0:0.005:0.3;
temp=15;

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
temp2 = temp+273;
Tref = 283;
E=0.6;
k=8.62e-5;


%% Mortality
%Hartvig
Hu = exp(0.063*(temp-10.0)) .* (u.*m.^(-n))/365;

%mizer
Zu = exp(0.063*(temp-10.0)) .* (u2.*m.^(-n))/365;

%J&C
Ju = exp((-1*E/k)*((1./temp2)-(1./Tref))) .* (U.*m.^(-N))/365;

%me
Cu = 0.1/365 * ones(size(m));

%% Annual Mortality
%Hartvig
HuA = exp(0.063*(temp-10.0)) .* (u.*m.^(-n));

%mizer
ZuA = exp(0.063*(temp-10.0)) .* (u2.*m.^(-n));

%J&C
JuA = exp((-1*E/k)*((1./temp2)-(1./Tref))) .* (U.*m.^(-N));

%me
CuA = 0.1 * ones(size(m));


%% Max ingestion
%Kiorboe & Hirst
Imax = exp(0.063*(temp-15.0)) .* 2.5 .* m.^(-0.51) .* 24e-3;

%Hartvig
Hcmax = exp(0.063*(temp-10.0)) .* (h.*m.^-n)/365;

%mizer
Zcmax = exp(0.063*(temp-10.0)) .* (h2.*m.^-n)/365;

%J&C
Jcmax = exp((-1*E/k)*((1./temp2)-(1./Tref))) .* (H.*m.^-N)/365;

%Me
Ccmax = (exp(0.063*(temp-10.0)) .* 20 .* m.^(-0.25)) ./365.0;

%% Respiration
%Hartvig
Hmet = exp(0.063*(temp-10.0)) .* k.*m.^-p ./365.0;

%mizer
Zmet = exp(0.063*(temp-10.0)) .* k2.*m.^-p ./365.0;

%J&C
Jmet = exp((-1*E/k)*((1./temp2)-(1./Tref))) .* K.*m.^-p ./365.0;

%Me
Cmet = 0.2 * (exp(0.08555*(temp-10.0)) .* 20 .* m.^(-0.175)) ./365.0;

%% Max clearance rate/search volume
%Kiorboe & Hirst
Fmax = exp(0.063*(temp-15.0)) .* 10^(3.24) .* m.^(-0.24) .* (24e-3/9);

%Hartvig
HA = exp(0.063*(temp-10.0)) .* (gamma*m.^-q)/365;

%mizer
ZA = exp(0.063*(temp-10.0)) .* (gamma2*m.^-q)/365;

%J&C
JA = exp(0.063*(temp-10.0)) .* (G*m.^-Q)/365;

%me
CA = (exp(0.063*(temp-10.0)) .* 70 .* m.^(-0.2)) ./365.0;

%%
figure(1)
subplot(2,2,1)
plot(log10(m),log10(Fmax),'.-r','MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(HA),'.-b','MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(JA),'.-k','MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(ZA),'.-','color',[0.5 0.5 0.5],'MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(CA),'.-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
% legend('K&H','Hart','J&C','mizer','me')
% legend('location','southwest')
title('Specific clearance rate')

subplot(2,2,2)
plot(log10(m),log10(Imax),'.-r','MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(Hcmax),'.-b','MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(Jcmax),'.-k','MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(Zcmax),'.-','color',[0.5 0.5 0.5],'MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(Ccmax),'.-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
ylabel('log10 ingestion (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('K&H','Hart','J&C','mizer','me')
legend('location','southwest')
title('Specific ingestion rate')

subplot(2,2,3)
plot(log10(m),log10(Hmet),'.-b','MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(Jmet),'.-k','MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(Zmet),'.-','color',[0.5 0.5 0.5],'MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(Cmet),'.-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
% legend('Hart','J&C','mizer','me')
% legend('location','southwest')
title('Specific respiration rate')

subplot(2,2,4)
plot(log10(m),log10(Hu),'.-b','MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(Ju),'.-k','MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(Zu),'.-','color',[0.5 0.5 0.5],'MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(Cu),'.-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
ylabel('log10 mortality rate (d^-^1)')
xlabel('log10 weight (gWW)')
% legend('Hart','J&C','mizer','me')
% legend('location','northeast')
title('Daily mortality rate')
print('-dpng','Rates_comp_mass.png')

%% TEMP-DEP
mass=[M_s; M_m; M_l];
temp = -2:30;
temp2 = temp+273;

for s=1:3
    m = mass(s);
    
    %% Mortality
    %Hartvig
    Hu = exp(0.063*(temp-10.0)) .* (u.*m.^(-n))/365;
    
    %mizer
    Zu = exp(0.063*(temp-10.0)) .* (u2.*m.^(-n))/365;
    
    %J&C
    Ju = exp((-1*E/k)*((1./temp2)-(1./Tref))) .* (U.*m.^(-N))/365;
    
    %me
    Cu = 0.1/365 * ones(size(temp));
    
    %% Annual Mortality
    %Hartvig
    HuA = exp(0.063*(temp-10.0)) .* (u.*m.^(-n));
    
    %mizer
    ZuA = exp(0.063*(temp-10.0)) .* (u2.*m.^(-n));
    
    %J&C
    JuA = exp((-1*E/k)*((1./temp2)-(1./Tref))) .* (U.*m.^(-N));
    
    %me
    CuA = 0.1 * ones(size(temp));
    
    
    %% Max ingestion
    %Kiorboe & Hirst
    Imax = exp(0.063*(temp-15.0)) .* 2.5 .* m.^(-0.51) .* 24e-3;
    
    %Hartvig
    Hcmax = exp(0.063*(temp-10.0)) .* (h.*m.^-n)/365;
    
    %mizer
    Zcmax = exp(0.063*(temp-10.0)) .* (h2.*m.^-n)/365;
    
    %J&C
    Jcmax = exp((-1*E/k)*((1./temp2)-(1./Tref))) .* (H.*m.^-N)/365;
    
    %Me
    Ccmax = (exp(0.063*(temp-10.0)) .* 20 .* m.^(-0.25)) ./365.0;
    
    %% Respiration
    %Hartvig
    Hmet = exp(0.063*(temp-10.0)) .* k.*m.^-p ./365.0;
    
    %mizer
    Zmet = exp(0.063*(temp-10.0)) .* k2.*m.^-p ./365.0;
    
    %J&C
    Jmet = exp((-1*E/k)*((1./temp2)-(1./Tref))) .* K.*m.^-p ./365.0;
    
    %Me
    Cmet = 0.2 * (exp(0.0905*(temp-10.0)) .* 20 .* m.^(-0.175)) ./365.0;
    
    %% Max clearance rate/search volume
    %Kiorboe & Hirst
    Fmax = exp(0.063*(temp-15.0)) .* 10^(3.24) .* m.^(-0.24) .* (24e-3/9);
    
    %Hartvig
    HA = exp(0.063*(temp-10.0)) .* (gamma*m.^-q)/365;
    
    %mizer
    ZA = exp(0.063*(temp-10.0)) .* (gamma2*m.^-q)/365;
    
    %J&C
    JA = exp(0.063*(temp-10.0)) .* (G*m.^-Q)/365;
    
    %me
    CA = (exp(0.063*(temp-10.0)) .* 70 .* m.^(-0.2)) ./365.0;
    
    %%
    f5=figure(5);
    subplot(3,4,4*s-3)
    plot(temp,log10(Fmax),'-r','MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(HA),'-b','MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(JA),'-k','MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(ZA),'-','color',[0.5 0.5 0.5],'MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(CA),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
    xlim([-2 30])
    if (s==1)
%         legend('K&H','Hart','J&C','mizer','me')
%         legend('location','northwest')
        title('Clearance')
    end
    
    subplot(3,4,4*s-2)
    plot(temp,log10(Imax),'-r','MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Hcmax),'-b','MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Jcmax),'-k','MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Zcmax),'-','color',[0.5 0.5 0.5],'MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Ccmax),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
    xlim([-2 30])
    if (s==1)
%         legend('K&H','Hart','J&C','mizer','me')
%         legend('location','northwest')
        title('Ingestion rate')
    end
    
    subplot(3,4,4*s-1)
    plot(temp,log10(Hmet),'-b','MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Jmet),'-k','MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Zmet),'-','color',[0.5 0.5 0.5],'MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Cmet),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
    xlim([-2 30])
    if (s==1)
%         legend('Hart','J&C','mizer','me')
%         legend('location','southeast')
        title('Respiration rate')
    end
    
    subplot(3,4,4*s)
    plot(temp,log10(Hu),'-b','MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Ju),'-k','MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Zu),'-','color',[0.5 0.5 0.5],'MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Cu),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
    xlim([-2 30])
    if (s==1)
%         legend('Hart','J&C','mizer','me')
%         legend('location','east')
        title('Mortality rate')
    end
    
    
end
subplot(3,4,1); ylabel('Small')
subplot(3,4,5); ylabel('Medium')
subplot(3,4,9); ylabel('Large')
subplot(3,4,10); xlabel('temp (^oC)')
print('-dpng','Rates_comp_temp.png')   

%%

