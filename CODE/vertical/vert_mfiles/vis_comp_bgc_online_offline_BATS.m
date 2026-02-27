%Compare BGC offline and online

clear
close all

%%
cfile ='NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/',...
    cfile,'/Cobalt_Feisty/'];

%% Offline
offp = '/Volumes/petrik-lab/Feisty/NC/MOM6-1D/BATS_vert/cobalt_only/';

%update with results from my run next week
%load([offp '20040101.ocean_cobalt_10yr_BATS_remy.mat'])
load([offp '20040101.ocean_cobalt_tracers_daily_means_10yr_remy.mat'])
load([offp '20040101.ocean_grid_12mo_BATS.mat'])

% %What is DP? and MP->MZprod
afDE  = mean(tB2(1099:2928));
afO2  = mean(tO2(1099:2928));
afNO3 = mean(tNO2(1099:2928));
afNH4 = mean(tNH2(1099:2928));
afDI = mean(tDP2(1099:2928));
afSP = mean(tSP2(1099:2928));
afLP = mean(tLP2(1099:2928));
afSZ = mean(tSZ2(1099:2928));

clear tB2 tO2 tNO2 tNH2 tDP2 tSP2 tLP2 

%%
load([offp '20040101.ocean_feisty_forcing_means_10yr_remy.mat'])

afMZ = mean(tMZ2(1099:2928));
afLZ = mean(tLZ2(1099:2928));
afMH = mean(tMH2(1099:2928));
afLH = mean(tLH2(1099:2928));

clear tSZ2 tMZ2 tLZ2 tMH2 tLH2 tMP2 tLP2

%% Online
onp = '/Volumes/petrik-lab/Feisty/NC/MOM6-1D/BATS_vert/cobalt_feisty/test_BATS_10yr/';

load([onp '20040101.ocean_cobalt_tracers_month_z_means.mat'])

anDE  = mean(tB2(37:96));
anO2  = mean(tO2(37:96));
anNO3 = mean(tNO2(37:96));
anNH4 = mean(tNH2(37:96));
anDI = mean(tDP2(37:96));
anSP = mean(tSP2(37:96));
anLP = mean(tLP2(37:96));
anSZ = mean(tSZ2(37:96));

clear tB2 tO2 tNO2 tNH2 tDP2 tSP2 tLP2 

%%
load([onp '20040101.ocean_feisty_forcing_means_10yr.mat'])

anMZ = mean(tMZ2(1099:2928));
anLZ = mean(tLZ2(1099:2928));
anMH = mean(tMH2(1099:2928));
anLH = mean(tLH2(1099:2928));

clear tSZ2 tMZ2 tLZ2 tMH2 tLH2 tMP2 tLP2

%% Percent Diffs Online - Offline / Offline
ltex = {'O2','POC','NO3','NH4','Diaz','SP','LP','SZ','MZ','LZ'};

mglob(1) = 100*(anO2 - afO2) ./ afO2;
mglob(2) = 100* (anDE - afDE) ./ afDE;
mglob(3) = 100*(anNO3 - afNO3) ./ afNO3;
mglob(4) = 100*(anNH4 - afNH4) ./ afNH4;
mglob(5) = 100* (anDI - afDI) ./ afDI;
mglob(6) = 100* (anSP - afSP) ./ afSP;
mglob(7) = 100* (anLP - afLP) ./ afLP;
mglob(8) = 100* (anSZ - afSZ) ./ afSZ;
mglob(9) = 100*(anMZ - afMZ) ./ afMZ;
mglob(10) = 100*(anLZ - afLZ) ./ afLZ;

dmMH = 100*(anMH - afMH) ./ afMH;
dmLH = 100*(anLH - afLH) ./ afLH;

save([offp 'ocean_cobalt_feisty_BATS_10y_Pdiff.mat'],...
    'mglob','dmMH','dmLH','ltex')
save([onp 'ocean_cobalt_feisty_BATS_10y_Pdiff.mat'],...
    'mglob','dmMH','dmLH','ltex')

%% Bar plot
f2 = figure('Units','inches','Position',[1 3 7 8]);
barh(mglob,'FaceColor',[0 0.5 0.75])
set(gca,'YTickLabel',ltex)
title('% Change at BATS')
xlim([-210 210])

print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_mean_Pdiff_all_BATS.png'])

%%
f3 = figure('Units','inches','Position',[1 3 7 8]);
barh(mglob(2:10),'FaceColor',[0 0.5 0.75])
set(gca,'YTickLabel',ltex(2:10))
title('% Change at BATS')
xlim([-210 210])

print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_mean_Pdiff_all_BATS_noO2.png'])

