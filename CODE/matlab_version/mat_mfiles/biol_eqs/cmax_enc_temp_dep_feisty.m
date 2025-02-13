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

%% Cmax & Met as fn of temp
Temps = [0,10,20,30];
for i=1:length(Temps)
    T=Temps(i);
    % Metabolism
    metT = (exp(kt*(T-10.0)) .* amet .* m.^(-bpow)) ./365.0;
    % Max ingestion
    cmaxT = (exp(0.063*(T-10.0)) .* 20 .* m.^(-0.25)) ./365.0;
    
    figure(1)
    subplot(2,2,i)
    plot(log10(m),log10(cmaxT),'.-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
    plot(log10(m),log10(metT),'-.','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
    if i>=3
        xlabel('log10 weight (gWW)')
    end
    if i==1
        legend('Cmax','Met')
    end
    if i==3
        ylabel('log10 ingestion/respiration (gWW gWW^-^1 d^-^1)')
    end
    title(['T = ' num2str(T) '^oC'])
    ylim([-3 0])
    
end

print('-dpng','Cmax_vs_met_mass_feisty_4temps.png')

%% Fn of prey density
prey = 0:0.05:5;
%temp = 15;
temp = [0,10,20,30];
frac = 1;
pref = 1;
pmat = repmat(prey,12,1);
pmat = reshape(pmat,3,4,101);
tmat = repmat(temp,3,1,101);
mmat = repmat(m,1,4,101);

% Max ingestion %Me
%cmaxP = (exp(0.063*(temp-10.0)) .* 20 .* m.^(-0.25)) ./365.0;
%cmax = (exp(kc*(temp-10.0)) .* h .* wgt^(-bcmx)) ./365.0;
cmaxP = (exp(0.063*(tmat-10.0)) .* 20 .* mmat.^(-0.25)) ./365.0;

% Encounter rate %Me
%A = (exp(0.063*(temp-10.0)) .* 70 .* m.^(-0.20)) ./365.0;
%A = (exp(ke*(temp-10.0)) .* gam .* wgt^(-benc)) ./365.0;
%encP = prey.*A.*frac.*pref;
A = (exp(0.063*(tmat-10.0)) .* 70 .* mmat.^(-0.20)) ./365.0;
encP = pmat.*A.*frac.*pref;

% Ingestion
%type II fn
ENC = encP;
%conP = repmat(cmaxP,1,length(encP)) .* encP ./ ((repmat(cmaxP,1,length(encP))) + ENC);
conP = cmaxP .* encP ./ (cmaxP + ENC);
%ENC = sum(enc,2); % total biomass encountered
%con = cmax .* enc(:,1) ./ (cmax + ENC); % Type II

%
figure(2)
subplot(3,1,1)
plot(prey,squeeze(conP(1,:,:)),'LineWidth',2); hold on;
%plot(prey,squeeze(cmaxP(1,:,:)),'--','color',[0 0.5 0.75],'LineWidth',2); hold on;
%legend('Con','Cmax')
legend('0','10','20','30')
legend('location','eastoutside')
title('Small')
xlim([0 4])

subplot(3,1,2)
plot(prey,squeeze(conP(2,:,:)),'LineWidth',2); hold on;
legend('0','10','20','30')
legend('location','eastoutside')
title('Medium')
xlim([0 4])

subplot(3,1,3)
plot(prey,squeeze(conP(3,:,:)),'LineWidth',2); hold on;
legend('0','10','20','30')
legend('location','eastoutside')
title('Large')
xlim([0 4])
xlabel('prey concentration')
print('-dpng','Con_vs_prey_mass_feisty_4temps.png')

%% TEMP-DEP
cmax = (exp(0.063*(temp-10.0)) .* 20 .* m.^(-0.25)) ./365.0;
temps = -2:40;
P = cmax/2;
P = repmat(P,1,length(temps));
P2 = 0.75*cmax;
P2 = repmat(P2,1,length(temps));

% Encounter rate %Me
A = (exp(0.063*(temps-10.0)) .* 70 .* m.^(-0.20)) ./365.0;
%A = (exp(ke*(temp-10.0)) .* gam .* wgt^(-benc)) ./365.0;
encT = P.*A.*frac.*pref;
encT2 = P2.*A.*frac.*pref;

% Max ingestion %Me
cmaxT = (exp(0.063*(temps-10.0)) .* 20 .* m.^(-0.25)) ./365.0;
%cmax = (exp(kc*(temp-10.0)) .* h .* wgt^(-bcmx)) ./365.0;

% Ingestion
%type II fn
ENC = encT;
conT = cmaxT .* encT ./ (cmaxT + ENC);
%ENC = sum(enc,2); % total biomass encountered
%con = cmax .* enc(:,1) ./ (cmax + ENC); % Type II
ENC2 = encT2;
conT2 = cmaxT .* encT2 ./ (cmaxT + ENC2);

%
figure(3)
subplot(3,2,1)
plot(temps,encT(1,:),'b','MarkerSize',15,'LineWidth',2); hold on;
plot(temps,cmaxT(1,:),'k','MarkerSize',15,'LineWidth',2); hold on;
plot(temps,conT(1,:),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
legend('Enc','Cmax','Con')
legend('location','northwest')
title('S, prey=1/2 Cmax')
xlim([-2 40])

subplot(3,2,3)
colororder({'b','k'})
yyaxis left
plot(temps,encT(2,:),'b','MarkerSize',15,'LineWidth',2); hold on;
plot(temps,conT(2,:),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
xlim([-2 40])
ylim([0 0.03])
ylabel('encounter/consumption (gWW gWW^-^1 d^-^1)')
yyaxis right
plot(temps,cmaxT(2,:),'k','MarkerSize',15,'LineWidth',2); hold on;
xlim([-2 40])
ylim([0 0.2])
title('M, prey=1/2 Cmax')

subplot(3,2,5)
yyaxis left
plot(temps,encT(3,:),'b','MarkerSize',15,'LineWidth',2); hold on;
plot(temps,conT(3,:),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
xlim([-2 40])
ylim([0 1.5e-3])
yyaxis right
colororder('k')
plot(temps,cmaxT(3,:),'k','MarkerSize',15,'LineWidth',2); hold on;
xlim([-2 40])
ylim([0 0.05])
title('L, prey=1/2 Cmax')

subplot(3,2,2)
plot(temps,encT2(1,:),'b','MarkerSize',15,'LineWidth',2); hold on;
plot(temps,cmaxT(1,:),'k','MarkerSize',15,'LineWidth',2); hold on;
plot(temps,conT2(1,:),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
title('S, prey=3/4 Cmax')
xlim([-2 40])

subplot(3,2,4)
colororder({'b','k'})
yyaxis left
plot(temps,encT2(2,:),'b','MarkerSize',15,'LineWidth',2); hold on;
plot(temps,conT2(2,:),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
xlim([-2 40])
ylim([0 0.03])
yyaxis right
plot(temps,cmaxT(2,:),'k','MarkerSize',15,'LineWidth',2); hold on;
xlim([-2 40])
ylim([0 0.2])
ylabel('max consumption (gWW gWW^-^1 d^-^1)')
title('M, prey=3/4 Cmax')

subplot(3,2,6)
yyaxis left
plot(temps,encT2(3,:),'b','MarkerSize',15,'LineWidth',2); hold on;
plot(temps,conT2(3,:),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
xlim([-2 40])
ylim([0 1.5e-3])
yyaxis right
colororder('k')
plot(temps,cmaxT(3,:),'k','MarkerSize',15,'LineWidth',2); hold on;
xlim([-2 40])
ylim([0 0.05])
title('L, prey=3/4 Cmax')
print('-dpng','Con_vs_temp_mass_feisty.png')


%% Temp-dep metab vs cmax
T = [-2 52];
% Metabolism
metT = (exp(kt*(T-10.0)) .* amet .* m.^(-bpow)) ./365.0;
% Max ingestion
cmaxT = (exp(0.063*(T-10.0)) .* 20 .* m.^(-0.25)) ./365.0;

figure(4)
subplot(3,2,1)
plot(T,log10(cmaxT(1,:)),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metT(1,:)),'-.','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
legend('Cmax','Met')
legend('location','northwest')
%ylabel('log_1_0 (ingestion/respiration (gWW gWW^-^1 d^-^1)')
title('S')
xlim([T(1) T(end)])

subplot(3,2,3)
plot(T,log10(cmaxT(2,:)),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metT(2,:)),'-.','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
ylabel('ingestion/respiration (gWW gWW^-^1 d^-^1)')
title('M')
xlim([T(1) T(end)])

subplot(3,2,5)
plot(T,log10(cmaxT(3,:)),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
plot(T,log10(metT(3,:)),'-.','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
xlabel('temperature')
title('L')
xlim([T(1) T(end)])
%ylim([-4 0])
print('-dpng','Cmax_vs_met_mass_feisty_temp.png')
