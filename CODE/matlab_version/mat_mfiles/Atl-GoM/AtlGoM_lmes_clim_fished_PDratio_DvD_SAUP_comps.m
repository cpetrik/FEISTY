% P:D ratio by LME 
% Climatology
% 150 years
% Saved as mat files
% Compare to Daniel's model results &SAUP

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cdir 'temp_100_1deg_ESM26_5yr_clim_191_195.mat'])
load([cdir 'btm_temp_1deg_ESM26_5yr_clim_191_195.mat'])
load([cpath 'LME_clim_temp_zoop_det.mat']);

ptemp_mean_clim=squeeze(nanmean(temp_100,1));
btemp_mean_clim=squeeze(nanmean(btm_temp,1));
tlme = lme_mask_onedeg;
AREA_OCN = max(area,1);

%% FEISTY 
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';
ppath = [pp cfile '/AtlGoM/'];
dpath = [dp cfile '/'];
load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat'],...
    'lme_mcatch','lme_area','rPD_biom','rPD_catch','rPD_catch_mtkm2');

lme_area_km2 = lme_area * 1e-6;

% FEISTY  LME biomass in MT
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
% MT/km2
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;

plme_rPDcatch = plme_Pmcatch ./ (plme_Pmcatch+plme_Dmcatch);

%% DvD on grid
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/DanielVD_PelDem/Colleen_modeledfish_LME.mat')
dlme_Pfrac = NaN*ones(180,360);
for L=5:7
    lid = find(tlme==L);
    dlme_Pfrac(lid) = FracLP(L);
end

%% SAUP
load([spath 'SAUP_LME_Catch_top10_Stock.mat']);
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')

sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

%on grid
tlme = lme_mask_onedeg;
sFracPD_grid = NaN*ones(180,360);
for L=5:7
    lid = find(tlme==L);
    sFracPD_grid(lid) = sFracPD(L);
end

%% Comparison stats
did=5:7;
notLELC = 5:7;

diffD = rPD_catch - dlme_Pfrac;
diffS = rPD_catch - sFracPD_grid;

%r
[rall,pall]=corr(FracLP(did),plme_rPDcatch(did));
[rPD,pPD]=corr(sFracPD(notLELC),plme_rPDcatch(notLELC));

%root mean square error
o=FracLP(did);
p=plme_rPDcatch(did);
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);

o=sFracPD(notLELC);
p=plme_rPDcatch(notLELC);
n = length(o);
num=nansum((p-o).^2);
rmsePD = sqrt(num/n);

%Fmed
Fall=10^(median(FracLP(did)-plme_rPDcatch(did)));
FPD=10^(median(sFracPD(notLELC)-plme_rPDcatch(notLELC)));

% Table
fish_stat(1,1) = rall;
fish_stat(2,1) = rmse;
fish_stat(3,1) = Fall;
fish_stat(1,2) = rPD;
fish_stat(2,2) = rmsePD;
fish_stat(3,2) = FPD;

Fstat = array2table(fish_stat,'RowNames',{'r','RMSE','Fmed'},...
    'VariableNames',{'DvD','SAU'});
writetable(Fstat,[dpath 'AtlGoM_LME_DvD_SAU_stats_' cfile '.csv'],'Delimiter',',','WriteRowNames',true)
save([dpath 'AtlGoM_LME_DvD_SAU_stats_' cfile '.mat'],'fish_stat')

%% Plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=15; %Set these bounds for your data
plotmaxlat=50;
plotminlon=-100;
plotmaxlon=-50;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac
% ENTER -100 TO MAP ORIGIN LONG

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

cmYOR=cbrewer('seq','YlOrRd',28);
cmRP=cbrewer('seq','RdPu',28);
cmPR=cbrewer('seq','PuRd',28);

x=0:0.1:1;

%% Subplot with maps and corr 
figure(1)
%SAU
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,diffS)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.475 0.6 0.025 0.34])
set(gcf,'renderer','painters')
text(-0.4,0.95,'A','FontSize',14)
%title('FEISTY  - SAU difference')

%DvD
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,diffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-0.4,0.95,'B','FontSize',14)
%title('FEISTY  - vanD difference')

%SAU corr
subplot('Position',[0.075 0.1 0.38 0.4])
plot(x,x,'--k');hold on;
scatter(sFracPD(notLELC),plme_rPDcatch(notLELC),50,'k','filled'); hold on;
text(0.725,0.55,['r = ' sprintf('%2.2f',rPD)])
text(0.725,0.49,['RMSE = ' sprintf('%2.2f',rmsePD)])
axis([0 1.0 0 1.0])
xlabel('SAU')
ylabel('FEISTY ')
text(0.05,0.925,'C','FontSize',14)
%title('Fraction Large Pelagics')

%DvD Corr
subplot('Position',[0.56 0.1 0.38 0.4])
plot(x,x,'--k');hold on;
scatter(FracLP(notLELC),plme_rPDcatch(notLELC),50,'k','filled'); hold on;
text(0.725,0.55,['r = ' sprintf('%2.2f',rall)])
text(0.725,0.49,['RMSE = ' sprintf('%2.2f',rmse)])
axis([0 1.0 0 1.0])
xlabel('vanD')
%ylabel('FEISTY ')
%title('Fraction Large Pelagics')
%stamp(cfile)
text(0.05,0.925,'D','FontSize',14)
print('-dpng',[ppath 'Clim_' harv '_AtlGoM_LME_fracPD_catch_SAUP_DvD_comp_subplot.png'])

%% Subplot with maps and corr just SAU
figure(2)
%SAU
subplot('Position',[0 0.5 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,diffS)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.39 0.55 0.025 0.3])
set(gcf,'renderer','painters')
text(-0.4,0.95,'A','FontSize',14)
%title('FEISTY  - SAU difference')

%Corr
subplot('Position',[0.55 0.55 0.4 0.4])
%SAU corr
plot(x,x,'--k');hold on;
scatter(sFracPD(notLELC),plme_rPDcatch(notLELC),50,'k','filled'); hold on;
text(0.725,0.55,['r = ' sprintf('%2.2f',rPD)])
text(0.725,0.49,['RMSE = ' sprintf('%2.2f',rmsePD)])
axis([0 1.0 0 1.0])
xlabel('SAU')
ylabel('FEISTY ')
text(0.05,0.925,'B','FontSize',14)
print('-dpng',[ppath 'Clim_' harv '_AtlGoM_LME_fracPD_catch_SAUP_comp_subplot.png'])

