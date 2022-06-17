% P:D ratio by LME 
% CORE-forced
% Observed effort
% Saved as mat files
% Compare to Daniel's model results &SAUP

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/SAUP/';

%% CORE-forced
load('/Volumes/MIP/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');
load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

ID = GRD.ID;

% ESM2M = same grid as CORE
gpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/cobalt_data/';
load([gpath 'hindcast_gridspec.mat'],'AREA_OCN');
load([gpath 'lme_mask_esm2m.mat']);
load([cpath 'LME_hist9095_temp_zoop_det.mat'],'lme_ptemp','lme_area');

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);
tlme = lme_mask_esm2m';

%% FEISTY LME biomass in MT/km2
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/CORE/'];
ppath = [pp cfile '/CORE/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

harv = 'fish_Fobs050_Pobs100_Dobs050'; %'fished_obs';

load([fpath 'LME_core_',harv,'_Catch_top10.mat'])

plme_mcatch = alme_mcatch10;
plme_Fmcatch = Flme_mcatch10;
plme_Pmcatch = Plme_mcatch10;
plme_Dmcatch = Dlme_mcatch10;

pFracPD = sFracPD;

l10p=log10(plme_mcatch);
l10pF=log10(plme_Fmcatch);
l10pP=log10(plme_Pmcatch);
l10pD=log10(plme_Dmcatch);

clear Flme_mcatch10 Plme_mcatch10 Dlme_mcatch10 sFracPD

%% DvD on grid
load('/Users/cpetrik/Dropbox/Princeton/FEISTY_other/DanielVD_PelDem/Colleen_modeledfish_LME.mat')
dlme_Pfrac = NaN*ones(size(tlme));
for L=1:63
    lid = find(tlme==L);
    dlme_Pfrac(lid) = FracLP(L);
end

%% SAUP
load([spath 'SAUP_LME_Catch_top10_Stock.mat']);
load(['/Users/cpetrik/Dropbox/Princeton/FEISTY_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')

sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

%on grid
sFracPD_grid = NaN*ones(size(tlme));
rPD_catch = NaN*ones(size(tlme));
for L=1:66
    lid = find(tlme==L);
    sFracPD_grid(lid) = sFracPD(L);
    rPD_catch(lid) = pFracPD(L);
end

%% Comparison stats
did=[1:61,63];
did2 = notLELC(notLELC<=63);

diffD = rPD_catch - dlme_Pfrac;
diffS = rPD_catch - sFracPD_grid;

%r
[rall,pall]=corr(FracLP(did),pFracPD(did));
[rall2,pall2]=corr(FracLP(did2),pFracPD(did2));
[rPD,pPD]=corr(sFracPD(notLELC),pFracPD(notLELC));

%root mean square error
o=FracLP(did);
p=pFracPD(did);
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);

o=FracLP(did2);
p=pFracPD(did2);
n = length(o);
num=nansum((p-o).^2);
rmse2 = sqrt(num/n);

o=sFracPD(notLELC);
p=pFracPD(notLELC);
n = length(o);
num=nansum((p-o).^2);
rmsePD = sqrt(num/n);

%Fmed
Fall=10^(median(FracLP(did)-pFracPD(did)));
Fall2=10^(median(FracLP(did2)-pFracPD(did2)));
FPD=10^(median(sFracPD(notLELC)-pFracPD(notLELC)));


% Table
fish_stat(1,1) = rall;
fish_stat(2,1) = rmse;
fish_stat(3,1) = Fall;
fish_stat(1,2) = rall2;
fish_stat(2,2) = rmse2;
fish_stat(3,2) = Fall2;
fish_stat(1,3) = rPD;
fish_stat(2,3) = rmsePD;
fish_stat(3,3) = FPD;

Fstat = array2table(fish_stat,'RowNames',{'r','RMSE','Fmed'},...
    'VariableNames',{'DvDAllLMEs','DvDnoLELC','SAUnoLELC'});
writetable(Fstat,[fpath 'core_',harv,'_LME_DvD_SAU_stats_' cfile '.csv'],'Delimiter',',','WriteRowNames',true)
save([fpath 'core_',harv,'_LME_DvD_SAU_stats_' cfile '.mat'],'fish_stat')

%% Plot info
[ni,nj]=size(geolon_t);
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90;
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

cmYOR=cbrewer('seq','YlOrRd',28,'PCHIP');
cmRP=cbrewer('seq','RdPu',28,'PCHIP');
cmPR=cbrewer('seq','PuRd',28,'PCHIP');

x=0:0.1:1;

load coastlines

%% Subplot with maps and corr and pval
figure(1)
%SAU
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffS)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.25 0.56 0.5 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
title('FEISTY - SAU difference')

%DvD
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffD)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('FEISTY - vanD difference')

%SAU corr
subplot('Position',[0.1 0.16 0.35 0.35])
plot(x,x,'--k');hold on;
scatter(sFracPD(notLELC),pFracPD(notLELC),20,lme_ptemp(notLELC,1),'filled'); hold on;
cmocean('thermal');
text(0.725,0.55,['r = ' sprintf('%2.2f',rPD)])
text(0.725,0.49,['(p = ' sprintf('%2.2f',pPD) ')'])
text(0.725,0.41,['RMSE = ' sprintf('%2.2f',rmsePD)])
axis([0 1.05 0 1.05])
xlabel('SAU')
ylabel('FEISTY')
%title('Fraction Large Pelagics')

%DvD Corr
subplot('Position',[0.575 0.16 0.35 0.35])
plot(x,x,'--k');hold on;
scatter(FracLP(did),pFracPD(did),20,lme_ptemp(did,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.25 0.05 0.5 0.025],'orientation','horizontal')
text(0.725,0.55,['r = ' sprintf('%2.2f',rall)])
text(0.725,0.49,['(p = ' sprintf('%2.2f',pall) ')'])
text(0.725,0.41,['RMSE = ' sprintf('%2.2f',rmse)])
axis([0 1.05 0 1.05])
xlabel('vanD')
ylabel('FEISTY')
%title('Fraction Large Pelagics')
%stamp(cfile)
print('-dpng',[ppath 'CORE_' harv '_LME_fracPD_catch_SAUP_DvD_comp_subplot.png'])


