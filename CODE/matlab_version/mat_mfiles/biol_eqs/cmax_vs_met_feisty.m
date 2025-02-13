% Compare diff Cmax, metab, A, ingest eqs

clear all
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


%% Max ingestion %Me
cmax = (exp(0.063*(temp-10.0)) .* 20 .* m.^(-0.25)) ./365.0;

%% Respiration %Me
met = (exp(0.0855*(temp-10.0)) .* 4 .* m.^(-0.175)) ./365.0;

%%
figure(1)
plot(log10(m),log10(cmax),'.-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(met),'-.','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
legend('Cmax','Met')
xlabel('log10 weight (gWW)')
ylabel('log10 ingestion/respiration (gWW gWW^-^1 d^-^1)')
title('FEISTY')
ylim([-3 0])
print('-dpng','Cmax_vs_met_mass_feisty.png')

%% TEMP-DEP
mass=[M_s; M_m; M_l];
temp = -2:40;
temp2 = temp+273;

for s=1:3
    m = mass(s);
    
    %% Max ingestion
    %Me
    Ccmax = (exp(0.063*(temp-10.0)) .* 20 .* m.^(-0.25)) ./365.0;
    
    %% Respiration
    %Me
    Cmet = (exp(0.0855*(temp-10.0)) .* 4 .* m.^(-0.175)) ./365.0;
    
    %%
    f9=figure(9);
    subplot(3,3,s)
    plot(temp,log10(Ccmax),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Cmet),'-.','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
    xlabel('temp (^oC)')
    xlim([-2 40])
    if (s==1)
        ylabel('log10 gWW gWW^-^1 d^-^1')
        title('S')
    elseif (s==2)
        title('M')
    else
        title('L')
    end
    
    f2=figure(2);
    subplot(2,2,s)
    plot(temp,log10(Ccmax),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Cmet),'-.','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
    xlabel('temp (^oC)')
    xlim([-2 40])
    if (s==1)
        ylabel('log10 gWW gWW^-^1 d^-^1')
        title('S')
    elseif (s==2)
        title('M')
    else
        ylabel('log10 gWW gWW^-^1 d^-^1')
        title('L')
    end
    
end
print(f9,'-dpng','Cmax_met_temp_FEISTY.png')

%%
figure(2)
subplot(2,2,4)
plot(log10(mass),log10(cmax),'.-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
plot(log10(mass),log10(met),'-.','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
legend('Cmax','Met')
xlabel('log10 weight (gWW)')
% ylabel('log10 ingestion/respiration (gWW gWW^-^1 d^-^1)')
% title('FEISTY')
ylim([-3 0])
print(f2,'-dpng','Cmax_met_mass_temp_FEISTY.png')
