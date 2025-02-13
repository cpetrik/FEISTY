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
Cmet = (exp(0.0855*(temp-10.0)) .* 4 .* m.^(-0.175)) ./365.0;

%%
figure(1)
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
print('-dpng','Ingestion_comp_mass_new_sizes_v4.png')

figure(2)
plot(log10(m),log10(Hmet),'.-b','MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(Jmet),'.-k','MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(Zmet),'.-','color',[0.5 0.5 0.5],'MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(Cmet),'.-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('Hart','J&C','mizer','me')
legend('location','southwest')
title('Specific respiration rate')
print('-dpng','Respiration_comp_mass_new_sizes_v4.png')


figure(3)
subplot(2,2,1)
plot(log10(m),log10(Hcmax),'.-b','MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(Hmet),'.-b','MarkerSize',15,'LineWidth',2); hold on;
ylabel('log10 ingestion/respiration (gWW gWW^-^1 d^-^1)')
title('Hartvig')
ylim([-3 0])

subplot(2,2,2)
plot(log10(m),log10(Jcmax),'.-k','MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(Jmet),'.-k','MarkerSize',15,'LineWidth',2); hold on;
title('J&C')
ylim([-3 0])

subplot(2,2,3)
plot(log10(m),log10(Zcmax),'.-','color',[0.5 0.5 0.5],'MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(Zmet),'.-','color',[0.5 0.5 0.5],'MarkerSize',15,'LineWidth',2); hold on;
title('mizer')
ylabel('log10 ingestion/respiration (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
ylim([-3 0])

subplot(2,2,4)
plot(log10(m),log10(Ccmax),'.-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
plot(log10(m),log10(Cmet),'.-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
xlabel('log10 weight (gWW)')
title('POEM')
ylim([-3 0])
print('-dpng','Cmax_vs_met_temp_new_sizes_v4.png')

%% TEMP-DEP
mass=[M_s; M_m; M_l];
temp = -2:30;
temp2 = temp+273;

for s=1:3
    m = mass(s);
    
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
    Cmet = (exp(0.0855*(temp-10.0)) .* 4 .* m.^(-0.175)) ./365.0;
    
    %%
    f4=figure(4);
    subplot(2,2,s)
    plot(temp,log10(Imax),'-r','MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Hcmax),'-b','MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Jcmax),'-k','MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Zcmax),'-','color',[0.5 0.5 0.5],'MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Ccmax),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
    ylabel('log10 ingestion (gWW gWW^-^1 d^-^1)')
    xlabel('temp (^oC)')
    xlim([-2 30])
    if (s==1)
    legend('K&H','Hart','J&C','mizer','me')
    legend('location','northwest')
    end
    if (s==1)
        title('S Specific ingestion rate')
    elseif (s==2)
        title('M Specific ingestion rate')
    else
        title('L Specific ingestion rate')
    end
    
    f5=figure(5);
    subplot(2,2,s)
    plot(temp,log10(Hmet),'-b','MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Jmet),'-k','MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Zmet),'-','color',[0.5 0.5 0.5],'MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Cmet),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
    ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
    xlabel('temp (^oC)')
    xlim([-2 30])
    if (s==1)
    legend('Hart','J&C','mizer','me')
    legend('location','southeast')
    end
    if (s==1)
        title('S Specific respiration rate')
    elseif (s==2)
        title('M Specific respiration rate')
    else
        title('L Specific respiration rate')
    end
    
    f6=figure(6);
    subplot(2,2,s)
    plot(temp,log10(Hcmax),'-b','MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Hmet),'-b','MarkerSize',15,'LineWidth',2); hold on;
    ylabel('log10 gWW gWW^-^1 d^-^1')
    xlabel('temp (^oC)')
    xlim([-2 30])
    if (s==1)
        title('S Hartvig')
    elseif (s==2)
        title('M Hartvig')
    else
        title('L Hartvig')
    end
    
    f7=figure(7);
    subplot(2,2,s)
    plot(temp,log10(Jcmax),'-k','MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Jmet),'-k','MarkerSize',15,'LineWidth',2); hold on;
    ylabel('log10 gWW gWW^-^1 d^-^1')
    xlabel('temp (^oC)')
    xlim([-2 30])
    if (s==1)
        title('S J&C')
    elseif (s==2)
        title('M J&C')
    else
        title('L J&C')
    end
    
    f8=figure(8);
    subplot(2,2,s)
    plot(temp,log10(Zcmax),'-','color',[0.5 0.5 0.5],'MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Zmet),'-','color',[0.5 0.5 0.5],'MarkerSize',15,'LineWidth',2); hold on;
    ylabel('log10 gWW gWW^-^1 d^-^1')
    xlabel('temp (^oC)')
    xlim([-2 30])
    if (s==1)
        title('S mizer')
    elseif (s==2)
        title('M mizer')
    else
        title('L mizer')
    end
    
    f9=figure(9);
    subplot(2,2,s)
    plot(temp,log10(Ccmax),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
    plot(temp,log10(Cmet),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
    ylabel('log10 gWW gWW^-^1 d^-^1')
    xlabel('temp (^oC)')
    xlim([-2 30])
    if (s==1)
        title('S POEM')
    elseif (s==2)
        title('M POEM')
    else
        title('L POEM')
    end
    
end
print(f4,'-dpng','Ingestion_comp_temp_new_sizes_v4.png')
print(f5,'-dpng','Respiration_comp_temp_new_sizes_v4.png')
print(f6,'-dpng','Cmax_met_temp_new_sizes_v4_Hartvig.png')
print(f7,'-dpng','Cmax_met_temp_new_sizes_v4_JC15.png')
print(f8,'-dpng','Cmax_met_temp_new_sizes_v4_mizer.png')
print(f9,'-dpng','Cmax_met_temp_new_sizes_v4_POEM.png')
