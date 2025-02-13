% Test diff half-sat coeffs for consumption eq
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

%% Fn of prey density
prey = 0:0.05:5;
temp = 10;
frac = 1;
pref = 1;

% Encounter rate %Me
A = (exp(0.063*(temp-10.0)) .* 70 .* m.^(-0.20)) ./365.0;
%A = (exp(ke*(temp-10.0)) .* gam .* wgt^(-benc)) ./365.0;
encP = prey.*A.*frac.*pref;

% Max ingestion %Me
cmaxP = (exp(0.063*(temp-10.0)) .* 20 .* m.^(-0.25)) ./365.0;
%cmax = (exp(kc*(temp-10.0)) .* h .* wgt^(-bcmx)) ./365.0;

% Ingestion
%type II fn
ENC = encP;
conP = repmat(cmaxP,1,length(encP)) .* encP ./ ((repmat(cmaxP,1,length(encP))) + ENC);
%ENC = sum(enc,2); % total biomass encountered
%con = cmax .* enc(:,1) ./ (cmax + ENC); % Type II

%
figure(2)
subplot(3,1,1)
plot(prey,(conP(1,:)),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
plot(prey,repmat(cmaxP(1),1,length(encP)),'--','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
%legend('Con','Cmax')
title('S, T=10^oC')
xlim([0 4])

subplot(3,1,2)
plot(prey,(conP(2,:)),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
plot(prey,repmat(cmaxP(2),1,length(encP)),'--','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
%legend('Con','Cmax')
ylabel('ingestion (gWW gWW^-^1 d^-^1)')
title('M, T=10^oC')
xlim([0 4])

subplot(3,1,3)
plot(prey,(conP(3,:)),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
plot(prey,repmat(cmaxP(3),1,length(encP)),'--','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
%legend('Con','Cmax')
xlabel('prey concentration')
title('L, T=10^oC')
xlim([0 4])
print('-dpng','Con_vs_prey_mass_feisty_10C.png')

%% diff half sats
enc = 0:0.01:2;
Ccmax = cmaxP;
con0 = enc ./ (Ccmax(1)+enc);
con1 = enc ./ (1.5*Ccmax(1)+enc);
con2 = enc ./ (2.0*Ccmax(1)+enc);
con3 = enc ./ (1.75*Ccmax(1)+enc);
con4 = enc ./ (1.90*Ccmax(1)+enc);
%conP = repmat(cmaxP,1,length(encP)) .* encP ./ ((repmat(cmaxP,1,length(encP))) + ENC);

figure(4)
%subplot(2,2,1)
plot(enc,con0,'k'); hold on;
plot(enc,con1,'b'); hold on;
plot(enc,con2,'r'); hold on;
% plot(enc,con3,'c'); hold on;
% plot(enc,con4,'m'); hold on;
%xlim([0 1])

