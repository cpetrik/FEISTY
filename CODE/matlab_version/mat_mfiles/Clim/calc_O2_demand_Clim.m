% Calc O2 demand in FEISTY
% Plot demand vs. temp (mean over Clim)
% Contour plot colors with volume of ocean with that demand/temp comb

clear all
close all

%% Params
%!Individual Mass (g) = geometric mean
M_s = 10^((log10(0.001)+log10(0.5))/2);  %0.0224
M_m = 10^((log10(0.5)+log10(250))/2);    %11.1803
M_l = 10^((log10(250)+log10(125000))/2); %5.5902e3

%%%! Metabolism constants (activity and basal)
amet = 4;       % coeff on met
kt = 0.0855;    % coeff on met T-dep fn (orig 0.063) %0.0855
bpow = 0.175;   % power on metab fn (orig 0.25)

%% Climatol
cpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/grid_cobalt/';
load([cpath 'esm26_area_1deg.mat']);
AREA_OCN = max(area,1);

load('/Volumes/MIP/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_gridspec.mat','ID');
load('/Volumes/MIP/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_grid.mat','GRD');
NX = length(GRD.Z);
%ID = 1:NX;
vAREA = AREA_OCN(ID);

load('/Volumes/MIP/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_daily.mat','COBALT');

%temp: COBALT.Tp & COBALT.Tb(ID,DY)
%area: area(ID)

%%
DAYS = 365;
Sml_p = nan*ones(NX,DAYS);
Med_p = nan*ones(NX,DAYS);
Lrg_p = nan*ones(NX,DAYS);
Med_d = nan*ones(NX,DAYS);
Lrg_d = nan*ones(NX,DAYS);

%% Calc metab at each grid cell over climatol
Sml_p = (exp(kt*(COBALT.Tp-10.0)) .* amet .* M_s.^(-bpow)) ./365.0;
Med_p = (exp(kt*(COBALT.Tp-10.0)) .* amet .* M_m.^(-bpow)) ./365.0;
Lrg_p = (exp(kt*(COBALT.Tp-10.0)) .* amet .* M_l.^(-bpow)) ./365.0;

Med_d = (exp(kt*(COBALT.Tb-10.0)) .* amet .* M_m.^(-bpow)) ./365.0;
Lrg_d = (exp(kt*(COBALT.Tb-10.0)) .* amet .* M_l.^(-bpow)) ./365.0;

%% Convert from gWWC d-1 to mol O2 d-1
%32.00 g O2 in 1 mol O2
%106/138 mol C in 1 mol O2 (Redfield ratio) OR SIMPLIFY TO 1:1
%12.01 g C in 1 mol C
%1 g dry W in 9 g wet

%1 g dry in 9 g wet; 12.01 g C in 1 mol C; 106/138 mol C in 1 mol O2
Omet_sp = Sml_p * (1/9.0) * (1/12.01) * (106/138);
Omet_mp = Med_p * (1/9.0) * (1/12.01) * (106/138);
Omet_lp = Lrg_p * (1/9.0) * (1/12.01) * (106/138);
Omet_md = Med_d * (1/9.0) * (1/12.01) * (106/138);
Omet_ld = Lrg_d * (1/9.0) * (1/12.01) * (106/138);

%% Take mean
Dmet_sp = nanmean(Omet_sp,2);
Dmet_mp = nanmean(Omet_mp,2);
Dmet_lp = nanmean(Omet_lp,2);
Dmet_md = nanmean(Omet_md,2);
Dmet_ld = nanmean(Omet_ld,2);

Dmet_sp(:,2) = nanmean(COBALT.Tp,2);
Dmet_mp(:,2) = nanmean(COBALT.Tp,2);
Dmet_lp(:,2) = nanmean(COBALT.Tp,2);
Dmet_md(:,2) = nanmean(COBALT.Tb,2);
Dmet_ld(:,2) = nanmean(COBALT.Tb,2);

%% Assign each grid cell its area
Dmet_sp(:,3) = vAREA;
Dmet_mp(:,3) = vAREA;
Dmet_lp(:,3) = vAREA;
Dmet_md(:,3) = vAREA;
Dmet_ld(:,3) = vAREA;

%% Fig
ppath = '/Users/cpetrik/Dropbox/Princeton/MAPP-METF/NCAR3/metab_index/O2_methods/';

figure(1)
subplot(3,2,1)
scatter(Dmet_sp(:,2),Dmet_sp(:,1),'filled'); hold on;
xlabel('Pelagic temperature')
ylabel('O_2 demand (mol O_2 d^-^1)')
title('Small fish in pelagic')

subplot(3,2,2)
scatter(Dmet_sp(:,2),Dmet_sp(:,1),'b','filled'); hold on;
scatter(Dmet_mp(:,2),Dmet_mp(:,1),'k','filled'); hold on;
scatter(Dmet_lp(:,2),Dmet_lp(:,1),'r','filled'); hold on;
legend('S','M','L')
legend('location','northwest')
xlabel('Pelagic temperature')
ylabel('O_2 demand (mol O_2 d^-^1)')
title('All fish')

subplot(3,2,3)
scatter(Dmet_mp(:,2),Dmet_mp(:,1),'filled'); hold on;
xlabel('Pelagic temperature')
ylabel('O_2 demand (mol O_2 d^-^1)')
title('Medium fish in pelagic')

subplot(3,2,4)
scatter(Dmet_md(:,2),Dmet_md(:,1),'filled'); hold on;
xlabel('Bottom temperature')
ylabel('O_2 demand (mol O_2 d^-^1)')
title('Medium demersal fish')

subplot(3,2,5)
scatter(Dmet_lp(:,2),Dmet_lp(:,1),'filled'); hold on;
xlabel('Pelagic temperature')
ylabel('O_2 demand (mol O_2 d^-^1)')
title('Large pelagic fish')

subplot(3,2,6)
scatter(Dmet_ld(:,2),Dmet_ld(:,1),'filled'); hold on;
xlabel('Bottom temperature')
ylabel('O_2 demand (mol O_2 d^-^1)')
title('Large demersal fish')

print('-dpng',[ppath 'Clim_O2_demand_vs_temp_size_habitat_molO2day.png'])




