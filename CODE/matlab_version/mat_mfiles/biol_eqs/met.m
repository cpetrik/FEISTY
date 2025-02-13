%Metabolism

clear all
close all

L_s = 10^((log10(2)+log10(20))/2);
L_m = 10^((log10(20)+log10(200))/2);
L_l = 10^((log10(200)+log10(2000))/2);

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;

fcrit=0.10;

temp=-2:30;

%! Specific ingestion rate from Kiorboe & Hirst (g/g/day)
cmaxS = (exp(0.063*(temp-15.0)) * 10^(0.4) * M_s^(-0.51)) .* 24e-3;
cmaxM = (exp(0.063*(temp-15.0)) * 10^(0.4) * M_m^(-0.51)) .* 24e-3;
cmaxL = (exp(0.063*(temp-15.0)) * 10^(0.4) * M_l^(-0.51)) .* 24e-3;
%Metabolism
basS = fcrit * cmaxS;
basM = fcrit * cmaxM;
basL = fcrit * cmaxL;

%%
cmax=[cmaxS;cmaxM;cmaxL];
met=[basS;basM;basL];
mass=[M_s;M_m;M_l];

figure
%subplot(2,2,1)
plot(temp,log10(basS),'r','LineWidth',2); hold on;
plot(temp,log10(basM),'k','LineWidth',2); hold on;
plot(temp,log10(basL),'b','LineWidth',2); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('temperature (^oC)')
legend('S','M','L')
legend('location','northwest')
title('Specific metabolism')
print('-dpng','Metab_fcrit10_temp.png')

figure
plot(log10(M_s),log10(basS(18)),'.k','MarkerSize',25); hold on;
plot(log10(M_m),log10(basM(18)),'.k','MarkerSize',25); hold on;
plot(log10(M_l),log10(basL(18)),'.k','MarkerSize',25); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('log10 mass (gWW)')
title('Specific metabolism')
print('-dpng','Metab_fcrit10_mass.png')

figure
%subplot(2,2,1)
plot(temp,log10(cmaxS),'r','LineWidth',2); hold on;
plot(temp,log10(cmaxM),'k','LineWidth',2); hold on;
plot(temp,log10(cmaxL),'b','LineWidth',2); hold on;
ylabel('log10 max ingestion (gWW gWW^-^1 d^-^1)')
xlabel('temperature (^oC)')
legend('S','M','L')
legend('location','northwest')
title('Specific ingestion')
print('-dpng','Ingest_temp.png')

figure
plot(log10(M_s),log10(cmaxS(18)),'.k','MarkerSize',25); hold on;
plot(log10(M_m),log10(cmaxM(18)),'.k','MarkerSize',25); hold on;
plot(log10(M_l),log10(cmaxL(18)),'.k','MarkerSize',25); hold on;
ylabel('log10 max ingestion (gWW gWW^-^1 d^-^1)')
xlabel('log10 mass (gWW)')
title('Specific ingestion')
print('-dpng','Ingest_mass.png')

