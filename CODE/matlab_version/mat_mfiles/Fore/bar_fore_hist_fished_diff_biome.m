% Visualize difference between
% ESM2M Hindcast of 1951-2000 and Forecast of 2051-2100
% Biome size and biomass in each biome

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

load([bpath 'COBALT_biomes_50yr_means_5100.mat'],'biome_hist','biome_fore');

load([cpath 'hindcast_gridspec.mat'],'AREA_OCN','geolon_t','geolat_t');
AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);

% colors
load('cmap_ppt_angles.mat')
cmap1(1,:)=cmap_ppt(1,:);
cmap1(2,:)=cmap_ppt(3,:);
cmap1(3,:)=cmap_ppt(5,:);
cmap1(4,:)=[1 1 1];

cmap2(1,:)=cmap_ppt(1,:);
cmap2(2,:)=cmap_ppt(3,:);
cmap2(3,:)=cmap_ppt(5,:);

cmap4(1,:)=cmap_ppt(1,:);
cmap4(2,:)=cmap_ppt(3,:);
cmap4(3,:)=cmap_ppt(5,:);
cmap4(4,:)=[0 0 0];

cmap3(1,:)=cmap_ppt(3,:);
cmap3(2,:)=cmap_ppt(1,:);
cmap3(3,:)=cmap_ppt(5,:);

%% Hindcast grid
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
%fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
fpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';

% Hindcast & Forecast together
load([fpath 'Means_hist_fore_',harv,'_cobalt_' cfile '.mat']);

cAll = cF+cP+cD;
hAll = hF+hP+hD;

zprod_hist = mzprod_hist + lzprod_hist;
zprod_fore = mzprod_fore + lzprod_fore;

ZpDet_hist = (zprod_hist./det_hist);
ZpDet_fore = (zprod_fore./det_fore);

%% area-integrated biomass
H_AF = hF .* AREA_OCN;
H_AP = hP .* AREA_OCN;
H_AD = hD .* AREA_OCN;
H_AA = hAll .* AREA_OCN;

F_AF = cF .* AREA_OCN;
F_AP = cP .* AREA_OCN;
F_AD = cD .* AREA_OCN;
F_AA = cAll .* AREA_OCN;

%% Calc area and change
% biome area
hLC = (biome_hist==1);
hCS = (biome_hist==2);
hSS = (biome_hist==3);

fLC = (biome_fore==1);
fCS = (biome_fore==2);
fSS = (biome_fore==3);

hLCA = hLC .* AREA_OCN;
hCSA = hCS .* AREA_OCN;
hSSA = hSS .* AREA_OCN;

fLCA = fLC .* AREA_OCN;
fCSA = fCS .* AREA_OCN;
fSSA = fSS .* AREA_OCN;

delLCA = sum(fLCA(:)) - sum(hLCA(:));
delCSA = sum(fCSA(:)) - sum(hCSA(:));
delSSA = sum(fSSA(:)) - sum(hSSA(:));

tbar(1,1) = delLCA;
tbar(1,2) = delCSA;
tbar(1,3) = delSSA;
tbar(2,1) = delLCA/sum(hLCA(:));
tbar(2,2) = delCSA/sum(hCSA(:));
tbar(2,3) = delSSA/sum(hSSA(:));

%% Global totals
gbar(1,1) = nansum(F_AF(:)) - nansum(H_AF(:));
gbar(2,1) = nansum(F_AP(:)) - nansum(H_AP(:));
gbar(3,1) = nansum(F_AD(:)) - nansum(H_AD(:));
gbar(4,1) = nansum(F_AA(:)) - nansum(H_AA(:));
gbar(1,2) = (nansum(F_AF(:)) - nansum(H_AF(:))) / nansum(H_AF(:));
gbar(2,2) = (nansum(F_AP(:)) - nansum(H_AP(:))) / nansum(H_AP(:));
gbar(3,2) = (nansum(F_AD(:)) - nansum(H_AD(:))) / nansum(H_AD(:));
gbar(4,2) = (nansum(F_AA(:)) - nansum(H_AA(:))) / nansum(H_AA(:));

% biomass per area
gbar(1,3) = (nanmean(cF(:))-nanmean(hF(:)));
gbar(2,3) = (nanmean(cP(:))-nanmean(hP(:)));
gbar(3,3) = (nanmean(cD(:))-nanmean(hD(:)));
gbar(4,3) = (nanmean(cAll(:))-nanmean(hAll(:)));
gbar(1,4) = gbar(1,3) ./ nanmean(hF(:));
gbar(2,4) = gbar(2,3) ./ nanmean(hP(:));
gbar(3,4) = gbar(3,3) ./ nanmean(hD(:));
gbar(4,4) = gbar(4,3) ./ nanmean(hAll(:));

%% Each biome separate
% area-integrated biomass
tbar(3,1) = nansum(F_AF(fLC)) - nansum(H_AF(hLC));
tbar(4,1) = nansum(F_AP(fLC)) - nansum(H_AP(hLC));
tbar(5,1) = nansum(F_AD(fLC)) - nansum(H_AD(hLC));
tbar(6,1) = nansum(F_AA(fLC)) - nansum(H_AA(hLC));
tbar(7,1) = (nansum(F_AF(fLC)) - nansum(H_AF(hLC))) / nansum(H_AF(hLC));
tbar(8,1) = (nansum(F_AP(fLC)) - nansum(H_AP(hLC))) / nansum(H_AP(hLC));
tbar(9,1) = (nansum(F_AD(fLC)) - nansum(H_AD(hLC))) / nansum(H_AD(hLC));
tbar(10,1) = (nansum(F_AA(fLC)) - nansum(H_AA(hLC))) / nansum(H_AA(hLC));

tbar(3,2) = nansum(F_AF(fCS)) - nansum(H_AF(hCS));
tbar(4,2) = nansum(F_AP(fCS)) - nansum(H_AP(hCS));
tbar(5,2) = nansum(F_AD(fCS)) - nansum(H_AD(hCS));
tbar(6,2) = nansum(F_AA(fCS)) - nansum(H_AA(hCS));
tbar(7,2) = (nansum(F_AF(fCS)) - nansum(H_AF(hCS))) / nansum(H_AF(hCS));
tbar(8,2) = (nansum(F_AP(fCS)) - nansum(H_AP(hCS))) / nansum(H_AP(hCS));
tbar(9,2) = (nansum(F_AD(fCS)) - nansum(H_AD(hCS))) / nansum(H_AD(hCS));
tbar(10,2) = (nansum(F_AA(fCS)) - nansum(H_AA(hCS))) / nansum(H_AA(hCS));

tbar(3,3) = nansum(F_AF(fSS)) - nansum(H_AF(hSS));
tbar(4,3) = nansum(F_AP(fSS)) - nansum(H_AP(hSS));
tbar(5,3) = nansum(F_AD(fSS)) - nansum(H_AD(hSS));
tbar(6,3) = nansum(F_AA(fSS)) - nansum(H_AA(hSS));
tbar(7,3) = (nansum(F_AF(fSS)) - nansum(H_AF(hSS))) / nansum(H_AF(hSS));
tbar(8,3) = (nansum(F_AP(fSS)) - nansum(H_AP(hSS))) / nansum(H_AP(hSS));
tbar(9,3) = (nansum(F_AD(fSS)) - nansum(H_AD(hSS))) / nansum(H_AD(hSS));
tbar(10,3) = (nansum(F_AA(fSS)) - nansum(H_AA(hSS))) / nansum(H_AA(hSS));

% biomass per area
tbar(11,1) = (nanmean(cF(fLC))-nanmean(hF(hLC)));
tbar(12,1) = (nanmean(cP(fLC))-nanmean(hP(hLC)));
tbar(13,1) = (nanmean(cD(fLC))-nanmean(hD(hLC)));
tbar(14,1) = (nanmean(cAll(fLC)) -nanmean(hAll(hLC)));
tbar(15,1) = tbar(11,1) ./ nanmean(hF(hLC));
tbar(16,1) = tbar(12,1) ./ nanmean(hP(hLC));
tbar(17,1) = tbar(13,1) ./ nanmean(hD(hLC));
tbar(18,1) = tbar(14,1) ./ nanmean(hAll(hLC));

tbar(11,2) = (nanmean(cF(fCS))-nanmean(hF(hCS)));
tbar(12,2) = (nanmean(cP(fCS))-nanmean(hP(hCS)));
tbar(13,2) = (nanmean(cD(fCS))-nanmean(hD(hCS)));
tbar(14,2) = (nanmean(cAll(fCS)) -nanmean(hAll(hCS)));
tbar(15,2) = tbar(11,2) ./ nanmean(hF(hCS));
tbar(16,2) = tbar(12,2) ./ nanmean(hP(hCS));
tbar(17,2) = tbar(13,2) ./ nanmean(hD(hCS));
tbar(18,2) = tbar(14,2) ./ nanmean(hAll(hCS));

tbar(11,3) = (nanmean(cF(fSS))-nanmean(hF(hSS)));
tbar(12,3) = (nanmean(cP(fSS))-nanmean(hP(hSS)));
tbar(13,3) = (nanmean(cD(fSS))-nanmean(hD(hSS)));
tbar(14,3) = (nanmean(cAll(fSS)) -nanmean(hAll(hSS)));
tbar(15,3) = tbar(11,3) ./ nanmean(hF(hSS));
tbar(16,3) = tbar(12,3) ./ nanmean(hP(hSS));
tbar(17,3) = tbar(13,3) ./ nanmean(hD(hSS));
tbar(18,3) = tbar(14,3) ./ nanmean(hAll(hSS));

%% tables
gstats = gbar;
grows = {'F','P','D','All'};
gcols = {'DiffAmass','PdiffAmass','DiffMass','PdiffMass'};

tstats = tbar;
tcols = {'LC','ECCS','ECSS'};
trows = {'DiffArea','PdiffArea',...
    'DiffAmassF','DiffAmassP','DiffAmassD','DiffAmassAll',...
    'PdiffAmassF','PdiffAmassP','PdiffAmassD','PdiffAmassAll',...
    'DiffMassF','DiffMassP','DiffMassD','DiffMassAll',...
    'PdiffMassF','PdiffMassP','PdiffMassD','PdiffMassAll'};

Stab = array2table(tstats,'VariableNames',tcols,'RowNames',trows);
writetable(Stab,[fpath 'Hist_Fore_All_fish03_change_area_biome_biomass.csv'],...
    'Delimiter',',','WriteRowNames',true);

Gtab = array2table(gstats,'VariableNames',gcols,'RowNames',grows);
writetable(Gtab,[fpath 'Hist_Fore_All_fish03_change_biomass_global.csv'],...
    'Delimiter',',','WriteRowNames',true);

%% pure change area, biomass, and area-integrated biomass
figure(1)
subplot(2,2,1)
bar(tbar(1,:)*1e-6*1e-7); hold on;
colormap(cmap4)
ylim([-2 2])
set(gca,'XTickLabel',{'LC','ECCS','ECSS'})
ylabel('Area (km^2) x e7')
title('Global change in area')

subplot(2,2,3)
bar(gbar(:,3)); hold on;
colormap(cmap4)
%ylim([-12e13 2e13])
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Biomass (g m^-^2)')
title('Change in Biomass')

subplot('position',[0.6 0.33 0.35 0.35])
bar(gbar(:,1)*1e-13); hold on;
colormap(cmap4)
%ylim([-12e13 2e13])
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Area-integrated Biomass (g) x e13')
title('Change in Area-integrated Biomass')
print('-dpng',[pp 'Hist_Fore_All_fish03_area_biome_mass_AImass_global_change.png'])

%% percent change area, biomass, and area-integrated biomass
figure(2)
subplot(2,2,1)
bar(tbar(2,:)); hold on;
colormap(cmap4)
ylim([-0.22 0.22])
set(gca,'XTickLabel',{'LC','ECCS','ECSS'})
ylabel('Area')
title('Global %change in area')

subplot(2,2,3)
bar(gbar(:,4)); hold on;
colormap(cmap4)
ylim([-0.25 0])
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Biomass')
title('%Change in Biomass')

subplot('position',[0.6 0.33 0.35 0.35])
bar(gbar(:,2)); hold on;
colormap(cmap4)
ylim([-0.25 0])
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Area-integrated Biomass')
title('%Change in Area-integrated Biomass')
print('-dpng',[pp 'Hist_Fore_All_fish03_area_biome_mass_AImass_global_pchange.png'])

%% pure change area and area-integrated biomass
figure(3)
subplot(2,2,1)
bar(tbar(1,:)*1e-6,'k'); hold on;
%colormap('gray')
ylim([-2e7 2e7])
set(gca,'XTickLabel',{'LC','ECCS','ECSS'})
ylabel('Area (km^2)')
title('Global change in area')

subplot(2,2,2)
bar(tbar(3:6,1),'FaceColor',cmap_ppt(3,:)); hold on;
colormap(cmap4)
ylim([-15e13 2e13])
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Area-integrated Biomass (g)')
title('Change in LC areas')

subplot(2,2,3)
bar(tbar(3:6,2),'FaceColor',cmap_ppt(1,:)); hold on;
colormap(cmap4)
ylim([-15e13 2e13])
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Area-integrated Biomass (g)')
title('Change in ECCS areas')

subplot(2,2,4)
bar(tbar(3:6,3),'FaceColor',cmap_ppt(5,:)); hold on;
colormap(cmap4)
ylim([-15e13 2e13])
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Area-integrated Biomass (g)')
title('Change in ECSS areas')
print('-dpng',[pp 'Hist_Fore_All_fish03_area_biome_change_AIbiomass.png'])

%% percent change area and area-integrated biomass
figure(4)
subplot(2,2,1)
bar(tbar(2,:),'k'); hold on;
ylim([-0.25 0.25])
set(gca,'XTickLabel',{'LC','ECCS','ECSS'})
ylabel('Area')
title('Global %change in area')

subplot(2,2,2)
bar(tbar(7:10,1),'FaceColor',cmap_ppt(3,:)); hold on;
ylim([-0.25 0])
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Area-integrated Biomass')
title('%Change in LC areas')

subplot(2,2,3)
bar(tbar(7:10,2),'FaceColor',cmap_ppt(1,:)); hold on;
ylim([-0.25 0])
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Area-integrated Biomass')
title('%Change in ECCS areas')

subplot(2,2,4)
bar(tbar(7:10,3),'FaceColor',cmap_ppt(5,:)); hold on;
ylim([-0.25 0])
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Area-integrated Biomass')
title('%Change in ECSS areas')
print('-dpng',[pp 'Hist_Fore_All_fish03_area_biome_pchange_AIbiomass.png'])


%% pure change area and biomass
figure(5)
subplot(2,2,1)
bar(tbar(1,:)*1e-6,'k'); hold on;
%colormap('gray')
ylim([-2e7 2e7])
set(gca,'XTickLabel',{'LC','ECCS','ECSS'})
ylabel('Area (km^2)')
title('Global change in area')

subplot(2,2,2)
bar(tbar(11:14,1),'FaceColor',cmap_ppt(3,:)); hold on;
colormap(cmap4)
ylim([-0.8 0.2])
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Biomass (g m-^2)')
title('Change in LC areas')

subplot(2,2,3)
bar(tbar(11:14,2),'FaceColor',cmap_ppt(1,:)); hold on;
colormap(cmap4)
ylim([-0.8 0.2])
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Biomass (g m-^2)')
title('Change in ECCS areas')

subplot(2,2,4)
bar(tbar(11:14,3),'FaceColor',cmap_ppt(5,:)); hold on;
colormap(cmap4)
ylim([-0.8 0.2])
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Biomass (g m-^2)')
title('Change in ECSS areas')
print('-dpng',[pp 'Hist_Fore_All_fish03_area_biome_change_biomass.png'])

%% percent change area and biomass
figure(6)
subplot(2,2,1)
bar(tbar(2,:),'k'); hold on;
ylim([-0.25 0.25])
set(gca,'XTickLabel',{'LC','ECCS','ECSS'})
ylabel('Area')
title('Global %change in area')

subplot(2,2,2)
bar(tbar(15:18,1),'FaceColor',cmap_ppt(3,:)); hold on;
ylim([-0.4 0.1])
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Biomass (g m-^2)')
title('%Change in LC areas')

subplot(2,2,3)
bar(tbar(15:18,2),'FaceColor',cmap_ppt(1,:)); hold on;
ylim([-0.4 0.1])
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Biomass (g m-^2)')
title('%Change in ECCS areas')

subplot(2,2,4)
bar(tbar(15:18,3),'FaceColor',cmap_ppt(5,:)); hold on;
ylim([-0.4 0.1])
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Biomass (g m-^2)')
title('%Change in ECSS areas')
print('-dpng',[pp 'Hist_Fore_All_fish03_area_biome_pchange_biomass.png'])


%% percent change in biomass and area-integrated biomass
figure(7)
b=bar(gbar(:,[2,4])'); hold on;
b(1).FaceColor = cmap4(2,:);
b(2).FaceColor = cmap4(1,:);
b(3).FaceColor = cmap4(3,:);
b(4).FaceColor = cmap4(4,:);
%colormap(cmap4)
%ylim([-12e13 2e13])
set(gca,'XTickLabel',{'Biomass','Area-integrated Biomass'})
ylabel('%Change in Biomass')
title('%Change in Biomass')
print('-dpng',[pp 'Hist_Fore_All_fish03_mass_AImass_global_pchange.png'])


%% pure change area and area-integrated biomass
figure(8)
% subplot(2,2,1)
b=bar(tbar(3:6,:),'FaceColor',cmap_ppt(3,:)); hold on;
b(1).FaceColor = cmap4(2,:);
b(2).FaceColor = cmap4(1,:);
b(3).FaceColor = cmap4(3,:);
legend(tcols)
legend('location','southwest')
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Area-integrated Biomass (g)')
title('Change in Area-integrated Biomass')
print('-dpng',[pp 'Hist_Fore_All_fish03_change_AIbiomass_byType.png'])

%%
figure(9)
% subplot(2,2,1)
b=bar(tbar(3:6,:)','FaceColor',cmap_ppt(3,:)); hold on;
b(1).FaceColor = cmap4(2,:);
b(2).FaceColor = cmap4(1,:);
b(3).FaceColor = cmap4(3,:);
b(4).FaceColor = cmap4(4,:);
legend(grows)
legend('location','southwest')
set(gca,'XTickLabel',tcols)
ylabel('Area-integrated Biomass (g)')
title('Change in Area-integrated Biomass')
print('-dpng',[pp 'Hist_Fore_All_fish03_change_AIbiomass_byBiome.png'])

%% percent change area and area-integrated biomass
figure(10)
% subplot(2,2,1)
b=bar(tbar(7:10,:),'FaceColor',cmap_ppt(3,:)); hold on;
b(1).FaceColor = cmap4(2,:);
b(2).FaceColor = cmap4(1,:);
b(3).FaceColor = cmap4(3,:);
legend(tcols)
legend('location','southwest')
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Area-integrated Biomass')
title('% Change in Area-integrated Biomass')
print('-dpng',[pp 'Hist_Fore_All_fish03_pchange_AIbiomass_byType.png'])

%%
figure(11)
% subplot(2,2,1)
b=bar(tbar(7:10,:)','FaceColor',cmap_ppt(3,:)); hold on;
b(1).FaceColor = cmap4(2,:);
b(2).FaceColor = cmap4(1,:);
b(3).FaceColor = cmap4(3,:);
b(4).FaceColor = cmap4(4,:);
legend(grows)
legend('location','southwest')
set(gca,'XTickLabel',tcols)
ylabel('Area-integrated Biomass')
title('% Change in Area-integrated Biomass')
print('-dpng',[pp 'Hist_Fore_All_fish03_pchange_AIbiomass_byBiome.png'])


%% pure change biomass
figure(12)
% subplot(2,2,1)
b=bar(tbar(11:14,:),'FaceColor',cmap_ppt(3,:)); hold on;
b(1).FaceColor = cmap4(2,:);
b(2).FaceColor = cmap4(1,:);
b(3).FaceColor = cmap4(3,:);
legend(tcols)
legend('location','southwest')
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Biomass (g m^-^2)')
title('Change in Biomass')
print('-dpng',[pp 'Hist_Fore_All_fish03_change_biomass_byType.png'])

%%
figure(13)
% subplot(2,2,1)
b=bar(tbar(11:14,:)','FaceColor',cmap_ppt(3,:)); hold on;
b(1).FaceColor = cmap4(2,:);
b(2).FaceColor = cmap4(1,:);
b(3).FaceColor = cmap4(3,:);
b(4).FaceColor = cmap4(4,:);
legend(grows)
legend('location','southwest')
set(gca,'XTickLabel',tcols)
ylabel('Biomass (g m^-^2)')
title('Change in Biomass')
print('-dpng',[pp 'Hist_Fore_All_fish03_change_biomass_byBiome.png'])

%% percent change biomass
figure(14)
% subplot(2,2,1)
b=bar(tbar(7:10,:),'FaceColor',cmap_ppt(3,:)); hold on;
b(1).FaceColor = cmap4(2,:);
b(2).FaceColor = cmap4(1,:);
b(3).FaceColor = cmap4(3,:);
legend(tcols)
legend('location','southwest')
set(gca,'XTickLabel',{'F','P','D','All'})
ylabel('Biomass')
title('% Change in Biomass')
print('-dpng',[pp 'Hist_Fore_All_fish03_pchange_biomass_byType.png'])

%%
figure(15)
% subplot(2,2,1)
b=bar(tbar(7:10,:)','FaceColor',cmap_ppt(3,:)); hold on;
b(1).FaceColor = cmap4(2,:);
b(2).FaceColor = cmap4(1,:);
b(3).FaceColor = cmap4(3,:);
b(4).FaceColor = cmap4(4,:);
legend(grows)
legend('location','southwest')
set(gca,'XTickLabel',tcols)
ylabel('Biomass')
title('% Change in Biomass')
print('-dpng',[pp 'Hist_Fore_All_fish03_pchange_biomass_byBiome.png'])

