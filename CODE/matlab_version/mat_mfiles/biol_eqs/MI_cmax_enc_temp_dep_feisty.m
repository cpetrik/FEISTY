% Compare diff Cmax, metab, A, ingest eqs

clear all
close all

%!Individual Mass (g)
M_s = 10^((log10(0.001)+log10(0.5))/2);
M_m = 10^((log10(0.5)+log10(250))/2);
M_l = 10^((log10(250)+log10(125000))/2);
m=[M_s; M_m; M_l];

kt=0.0855;
amet = 4;
bpow = 0.175;

% 0.5 0.5 0;... %tan/army
%     0 0.7 0;...   %g
%     1 0 1;...     %m
%     1 0 0;...     %r
%     0.5 0 0;...   %maroon
cm10=[0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0.5 0.5 0.5; ...    %med grey
    0 0 0];...      %black

set(groot,'defaultAxesColorOrder',cm10);

%% MI params
load('MI_fntypes3_redux.mat')
%1:all, 2:D, 3:F, 4:P

bmet = Emet/ (8.6173e-5 * 273^2);
bsup = Esup/ (8.6173e-5 * 273^2);

%%  Met as fn of temp
Temps = [0,10,20,30];
for i=1:length(Temps)
    T=Temps(i);
    % Metabolism
    metT = (exp(kt*(T-10.0)) .* amet .* m.^(-bpow)) ./365.0;
    metTmi = (exp(bmet*(T-10.0)) .* alphaD .* m'.^(eps1)) ./365.0;
    
    figure(1)
    subplot(2,2,i)
    plot(log10(m),log10(metTmi),'-','LineWidth',2); hold on;
    plot(log10(m),log10(metT),'-k','LineWidth',2); hold on;
    if i>=3
        xlabel('log10 weight (gWW)')
    end
%     if i==1
%         legend('MImmet','MImetD','MImetF','MImetP','FEISTY')
%     end
    if i==3
        ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
    end
    title(['T = ' num2str(T) '^oC'])
    ylim([-3 0.5])
    
end
print('-dpng','MI_feisty_met_mass_4temps.png')

figure(2)
plot(log10(m),log10(metTmi),'-','LineWidth',2); hold on;
plot(log10(m),log10(metT),'-k','LineWidth',2); hold on;
legend('MImmet','MImetD','MImetF','MImetP','FEISTY')
print('-dpng','MI_feisty_met_mass_4temps_legend.png')
    
%% Fn of temp
temp = 0:5:35;
tmat = repmat(temp,3,1);
mmat = repmat(m,1,length(temp));
pO2 = 0.21;

% Max supply %Me=cmax
cmaxO = (exp(0.063*(tmat-10.0)) .* 20 .* mmat.^(-0.25)) ./365.0;
supA = (exp(bsup(1)*(tmat-10.0)) .* alphaS(1) .* m.^(del(1)) *pO2) ./365.0;
supD = (exp(bsup(2)*(tmat-10.0)) .* alphaS(2) .* m.^(del(2)) *pO2) ./365.0;
supF = (exp(bsup(3)*(tmat-10.0)) .* alphaS(3) .* m.^(del(3)) *pO2) ./365.0;
supP = (exp(bsup(4)*(tmat-10.0)) .* alphaS(4) .* m.^(del(4)) *pO2) ./365.0;

% Met
metO = (exp(kt*(tmat-10.0)) .* amet .* m.^(-bpow)) ./365.0;
metA = (exp(bmet(1)*(tmat-10.0)) .* alphaD(1) .* m.^(eps1(1))) ./365.0;
metD = (exp(bmet(2)*(tmat-10.0)) .* alphaD(2) .* m.^(eps1(2))) ./365.0;
metF = (exp(bmet(3)*(tmat-10.0)) .* alphaD(3) .* m.^(eps1(3))) ./365.0;
metP = (exp(bmet(4)*(tmat-10.0)) .* alphaD(4) .* m.^(eps1(4))) ./365.0;

%% Supp
figure(3)
subplot(3,1,1)
plot(temp,log10(supA(1,:)),'LineWidth',2); hold on;
plot(temp,log10(supD(1,:)),'LineWidth',2); hold on;
plot(temp,log10(supF(1,:)),'LineWidth',2); hold on;
plot(temp,log10(supP(1,:)),'LineWidth',2); hold on;
plot(temp,log10(cmaxO(1,:)),'LineWidth',2); hold on;
legend('MIall','MID','MIF','MIP','FEISTY')
legend('location','eastoutside')
title('Small')
text(37,0.2,['pO2 = ' num2str(pO2)])
%xlim([0 4])

subplot(3,1,2)
plot(temp,log10(supA(2,:)),'LineWidth',2); hold on;
plot(temp,log10(supD(2,:)),'LineWidth',2); hold on;
plot(temp,log10(supF(2,:)),'LineWidth',2); hold on;
plot(temp,log10(supP(2,:)),'LineWidth',2); hold on;
plot(temp,log10(cmaxO(2,:)),'LineWidth',2); hold on;
legend('MIall','MID','MIF','MIP','FEISTY')
legend('location','eastoutside')
title('Medium')
ylabel('Supply')
%xlim([0 4])

subplot(3,1,3)
plot(temp,log10(supA(3,:)),'LineWidth',2); hold on;
plot(temp,log10(supD(3,:)),'LineWidth',2); hold on;
plot(temp,log10(supF(3,:)),'LineWidth',2); hold on;
plot(temp,log10(supP(3,:)),'LineWidth',2); hold on;
plot(temp,log10(cmaxO(3,:)),'LineWidth',2); hold on;
legend('MIall','MID','MIF','MIP','FEISTY')
legend('location','eastoutside')
title('Large')
%xlim([0 4])
xlabel('Temperature')
print('-dpng','Sup_vs_temp_feisty_MI.png')

%% Met vs temp
figure(4)
subplot(3,1,1)
plot(temp,log10(metA(1,:)),'LineWidth',2); hold on;
plot(temp,log10(metD(1,:)),'LineWidth',2); hold on;
plot(temp,log10(metF(1,:)),'LineWidth',2); hold on;
plot(temp,log10(metP(1,:)),'LineWidth',2); hold on;
plot(temp,log10(metO(1,:)),'LineWidth',2); hold on;
legend('MIall','MID','MIF','MIP','FEISTY')
legend('location','eastoutside')
title('Small')
%xlim([0 4])

subplot(3,1,2)
plot(temp,log10(metA(2,:)),'LineWidth',2); hold on;
plot(temp,log10(metD(2,:)),'LineWidth',2); hold on;
plot(temp,log10(metF(2,:)),'LineWidth',2); hold on;
plot(temp,log10(metP(2,:)),'LineWidth',2); hold on;
plot(temp,log10(metO(2,:)),'LineWidth',2); hold on;
legend('MIall','MID','MIF','MIP','FEISTY')
legend('location','eastoutside')
title('Medium')
ylabel('Demand')
%xlim([0 4])

subplot(3,1,3)
plot(temp,log10(metA(3,:)),'LineWidth',2); hold on;
plot(temp,log10(metD(3,:)),'LineWidth',2); hold on;
plot(temp,log10(metF(3,:)),'LineWidth',2); hold on;
plot(temp,log10(metP(3,:)),'LineWidth',2); hold on;
plot(temp,log10(metO(3,:)),'LineWidth',2); hold on;
legend('MIall','MID','MIF','MIP','FEISTY')
legend('location','eastoutside')
title('Large')
%xlim([0 4])
xlabel('Temperature')
print('-dpng','Met_vs_temp_feisty_MI.png')

%% Temp-dep metab vs cmax
T = temp;

figure(5)
subplot(3,4,1)
plot(T,log10(cmaxO(1,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metO(1,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(supA(1,:)),'-','color',cm10(1,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metA(1,:)),'-.','color',cm10(1,:),'MarkerSize',15,'LineWidth',2); hold on;
title('S')
xlim([T(1) T(end)])
subplot(3,4,5)
plot(T,log10(cmaxO(2,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metO(2,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(supA(2,:)),'-','color',cm10(1,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metA(2,:)),'-.','color',cm10(1,:),'MarkerSize',15,'LineWidth',2); hold on;
ylabel('ingestion/respiration (gWW gWW^-^1 d^-^1)')
title('M')
xlim([T(1) T(end)])
subplot(3,4,9)
plot(T,log10(cmaxO(3,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metO(3,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(supA(3,:)),'-','color',cm10(1,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metA(3,:)),'-.','color',cm10(1,:),'MarkerSize',15,'LineWidth',2); hold on;
xlabel('temperature')
title('L')
xlim([T(1) T(end)])

subplot(3,4,2)
plot(T,log10(cmaxO(1,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metO(1,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(supD(1,:)),'-','color',cm10(2,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metD(1,:)),'-.','color',cm10(2,:),'MarkerSize',15,'LineWidth',2); hold on;
title('SD')
xlim([T(1) T(end)])
subplot(3,4,6)
plot(T,log10(cmaxO(2,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metO(2,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(supD(2,:)),'-','color',cm10(2,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metD(2,:)),'-.','color',cm10(2,:),'MarkerSize',15,'LineWidth',2); hold on;
title('MD')
xlim([T(1) T(end)])
subplot(3,4,10)
plot(T,log10(cmaxO(3,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metO(3,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(supD(3,:)),'-','color',cm10(2,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metD(3,:)),'-.','color',cm10(2,:),'MarkerSize',15,'LineWidth',2); hold on;
xlabel('temperature')
title('LD')
xlim([T(1) T(end)])
text(45,0.2,['pO2 = ' num2str(pO2)])

subplot(3,4,3)
plot(T,log10(cmaxO(1,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metO(1,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(supF(1,:)),'-','color',cm10(3,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metF(1,:)),'-.','color',cm10(3,:),'MarkerSize',15,'LineWidth',2); hold on;
title('SF')
xlim([T(1) T(end)])
subplot(3,4,7)
plot(T,log10(cmaxO(2,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metO(2,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(supF(2,:)),'-','color',cm10(3,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metF(2,:)),'-.','color',cm10(3,:),'MarkerSize',15,'LineWidth',2); hold on;
title('MF')
xlim([T(1) T(end)])
% subplot(3,4,11)
% plot(T,log10(cmaxO(3,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
% plot(T,log10(metO(3,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
% plot(T,log10(supF(3,:)),'-','color',cm10(3,:),'MarkerSize',15,'LineWidth',2); hold on;
% plot(T,log10(metF(3,:)),'-.','color',cm10(3,:),'MarkerSize',15,'LineWidth',2); hold on;
% xlabel('temperature')
% title('LF')
% xlim([T(1) T(end)])

subplot(3,4,4)
plot(T,log10(cmaxO(1,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metO(1,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(supP(1,:)),'-','color',cm10(4,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metP(1,:)),'-.','color',cm10(4,:),'MarkerSize',15,'LineWidth',2); hold on;
title('SP')
xlim([T(1) T(end)])
subplot(3,4,8)
plot(T,log10(cmaxO(2,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metO(2,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(supP(2,:)),'-','color',cm10(4,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metP(2,:)),'-.','color',cm10(4,:),'MarkerSize',15,'LineWidth',2); hold on;
title('MP')
xlim([T(1) T(end)])
subplot(3,4,12)
plot(T,log10(cmaxO(3,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metO(3,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(supP(3,:)),'-','color',cm10(4,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metP(3,:)),'-.','color',cm10(4,:),'MarkerSize',15,'LineWidth',2); hold on;
xlabel('temperature')
title('LP')
xlim([T(1) T(end)])
print('-dpng','Sup_vs_met_temp_feisty_MI_size.png')

%% Calc MI both ways
T0 = 273;
T = tmat+T0;
Tref = 15+T0; 
R=8.6173e-5; 
Patm=0.209;
atm2kPa=101.325;

FEmet = 0.0855 * R * T0^2;
FEsup = 0.063 * R * T0^2;
FE0 = FEmet - FEsup;
FalphaD = 4;
FalphaC = 20;
FalphaE = 70;
Falpha0 = FalphaC / FalphaD;
Fdel = -0.25;
Feps = -0.175;
Fn = Fdel - Feps;
MIn = del - eps1;

% from Ao and Eo
%Phi = A0 * w^n *pO2 .* (exp((E0/R)*((1/T) - (1/Tref))));
PhiO = Falpha0 .* mmat.^(Fn) ./ (exp((FE0./R)*((1./T) - (1./Tref))));
PhiA = Ao(1) .* mmat.^(MIn(1)) *pO2 ./ (exp((Eo(1)/R)*((1./T) - (1./Tref))));
PhiD = Ao(2) .* mmat.^(MIn(2)) *pO2 ./ (exp((Eo(2)/R)*((1./T) - (1./Tref))));
PhiF = Ao(3) .* mmat.^(MIn(3)) *pO2 ./ (exp((Eo(3)/R)*((1./T) - (1./Tref))));
PhiP = Ao(4) .* mmat.^(MIn(4)) *pO2 ./ (exp((Eo(4)/R)*((1./T) - (1./Tref))));

% from ind supp vs demand
PhiO2 = cmaxO ./ metO;
PhiA2 = supA ./ metA;
PhiD2 = supD ./ metD;
PhiF2 = supF ./ metF;
PhiP2 = supP ./ metP;

%% Temp-dep Phi

figure(6)
subplot(3,4,1)
plot(temp,log10(PhiO(1,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiO2(1,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiA(1,:)),'-','color',cm10(1,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiA2(1,:)),'-.','color',cm10(1,:),'MarkerSize',15,'LineWidth',2); hold on;
title('S')
xlim([temp(1) temp(end)])
subplot(3,4,5)
plot(temp,log10(PhiO(2,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiO2(2,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiA(2,:)),'-','color',cm10(1,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiA2(2,:)),'-.','color',cm10(1,:),'MarkerSize',15,'LineWidth',2); hold on;
ylabel('Metabolic Index')
title('M')
xlim([temp(1) temp(end)])
subplot(3,4,9)
plot(temp,log10(PhiO(3,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiO2(3,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiA(3,:)),'-','color',cm10(1,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiA2(3,:)),'-.','color',cm10(1,:),'MarkerSize',15,'LineWidth',2); hold on;
xlabel('temperature')
title('L')
xlim([temp(1) temp(end)])

subplot(3,4,2)
plot(temp,log10(PhiO(1,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiO2(1,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiD(1,:)),'-','color',cm10(2,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiD2(1,:)),'-.','color',cm10(2,:),'MarkerSize',15,'LineWidth',2); hold on;
title('SD')
xlim([temp(1) temp(end)])
subplot(3,4,6)
plot(temp,log10(PhiO(2,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiO2(2,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiD(2,:)),'-','color',cm10(2,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiD2(2,:)),'-.','color',cm10(2,:),'MarkerSize',15,'LineWidth',2); hold on;
title('MD')
xlim([temp(1) temp(end)])
subplot(3,4,10)
plot(temp,log10(PhiO(3,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiO2(3,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiD(3,:)),'-','color',cm10(2,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiD2(3,:)),'-.','color',cm10(2,:),'MarkerSize',15,'LineWidth',2); hold on;
xlabel('temperature')
title('LD')
xlim([temp(1) temp(end)])
text(45,0.5,['pO2 = ' num2str(pO2)])
text(45,1,['Tref = ' num2str(Tref-T0) 'C'])

subplot(3,4,3)
plot(temp,log10(PhiO(1,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiO2(1,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiF(1,:)),'-','color',cm10(3,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiF2(1,:)),'-.','color',cm10(3,:),'MarkerSize',15,'LineWidth',2); hold on;
title('SF')
xlim([temp(1) temp(end)])
subplot(3,4,7)
plot(temp,log10(PhiO(2,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiO2(2,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiF(2,:)),'-','color',cm10(3,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiF2(2,:)),'-.','color',cm10(3,:),'MarkerSize',15,'LineWidth',2); hold on;
title('MF')
xlim([temp(1) temp(end)])
% subplot(3,4,11)
% plot(temp,log10(PhiO(3,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
% plot(temp,log10(PhiO2(3,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
% plot(temp,log10(PhiF(3,:)),'-','color',cm10(3,:),'MarkerSize',15,'LineWidth',2); hold on;
% plot(temp,log10(PhiF2(3,:)),'-.','color',cm10(3,:),'MarkerSize',15,'LineWidth',2); hold on;
% xlabel('temperature')
% title('LF')
% xlim([temp(1) temp(end)])

subplot(3,4,4)
plot(temp,log10(PhiO(1,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiO2(1,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiP(1,:)),'-','color',cm10(4,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiP2(1,:)),'-.','color',cm10(4,:),'MarkerSize',15,'LineWidth',2); hold on;
title('SP')
xlim([temp(1) temp(end)])
subplot(3,4,8)
plot(temp,log10(PhiO(2,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiO2(2,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiP(2,:)),'-','color',cm10(4,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiP2(2,:)),'-.','color',cm10(4,:),'MarkerSize',15,'LineWidth',2); hold on;
title('MP')
xlim([temp(1) temp(end)])
subplot(3,4,12)
plot(temp,log10(PhiO(3,:)),'-','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiO2(3,:)),'-.','color',cm10(5,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiP(3,:)),'-','color',cm10(4,:),'MarkerSize',15,'LineWidth',2); hold on;
plot(temp,log10(PhiP2(3,:)),'-.','color',cm10(4,:),'MarkerSize',15,'LineWidth',2); hold on;
xlabel('temperature')
title('LP')
xlim([temp(1) temp(end)])
print('-dpng','MIeq_vs_MIcalc_temp_feisty_MI_size.png')
