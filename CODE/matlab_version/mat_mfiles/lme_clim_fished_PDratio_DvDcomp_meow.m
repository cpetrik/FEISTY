% P:D ratio by LME 
% Climatology
% 150 years
% Saved as mat files
% Compare to Daniel's model results MEOW mapped to LME

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
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

%% POEM
cfile = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';
ppath = [pp cfile '/'];
dpath = [dp cfile '/'];
load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat'],...
    'lme_mcatch','lme_area','rPD_biom','rPD_catch','rPD_catch_mtkm2');

lme_area_km2 = lme_area * 1e-6;

% POEM LME biomass in MT
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
% MT/km2
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;

plme_rPDcatch = plme_Pmcatch ./ (plme_Pmcatch+plme_Dmcatch);

%% DvD on grid
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/DanielVD_PelDem/DvD_results_LME.mat')
dlme_Pfrac = NaN*ones(180,360);
dlme_Sfrac = NaN*ones(180,360);
for L=1:63
    lid = find(tlme==L);
    dlme_Pfrac(lid) = mod_frac(L);
    dlme_Sfrac(lid) = obs_frac(L);
end

%% Comparison stats
keep=[1:61,63];
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
notLELC = notLELC(notLELC<=63);

%r
rall=corr(mod_frac(keep),plme_rPDcatch(keep));
rall2=corr(mod_frac(notLELC),plme_rPDcatch(notLELC));

robs=corr(obs_frac(keep),plme_rPDcatch(keep));
robs2=corr(obs_frac(notLELC),plme_rPDcatch(notLELC));

%root mean square error
o=mod_frac(keep);
p=plme_rPDcatch(keep);
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);

o=mod_frac(notLELC);
p=plme_rPDcatch(notLELC);
n = length(o);
num=nansum((p-o).^2);
rmse2 = sqrt(num/n);

o=obs_frac(keep);
p=plme_rPDcatch(keep);
n = length(o);
num=nansum((p-o).^2);
rmse_obs = sqrt(num/n);

o=obs_frac(notLELC);
p=plme_rPDcatch(notLELC);
n = length(o);
num=nansum((p-o).^2);
rmse_obs2 = sqrt(num/n);

%Fmed
Fall=10^(median(mod_frac(keep)-plme_rPDcatch(keep)));
Fall2=10^(median(mod_frac(notLELC)-plme_rPDcatch(notLELC)));
Fobs=10^(median(obs_frac(keep)-plme_rPDcatch(keep)));
Fobs2=10^(median(obs_frac(notLELC)-plme_rPDcatch(notLELC)));

% Table
fish_stat(1,1) = rall;
fish_stat(2,1) = rmse;
fish_stat(3,1) = Fall;
fish_stat(1,2) = rall2;
fish_stat(2,2) = rmse2;
fish_stat(3,2) = Fall2;
fish_stat(1,3) = robs;
fish_stat(2,3) = rmse_obs;
fish_stat(3,3) = Fobs;
fish_stat(1,4) = robs2;
fish_stat(2,4) = rmse_obs2;
fish_stat(3,4) = Fobs2;

Fstat = array2table(fish_stat,'RowNames',{'r','RMSE','Fmed'},...
    'VariableNames',{'mLMEs','mLELC','oLMEs','oLELC'});
writetable(Fstat,[dpath 'MEOW_LME_DvD_stats_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true)
save([dpath 'MEOW_LME_DvD_stats_' cfile '.mat'],'fish_stat')

%% Plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac
% ENTER -100 TO MAP ORIGIN LONG

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

revamp=colormap(cmocean('amp'));
revamp=flipud(revamp);
% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);
close all

x=0:0.1:1;

%% Figures
% Correlation
%Daniel's food web model
figure(1)
subplot(1,2,1)
plot(x,x,'--k');hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(mod_frac(lme),plme_rPDcatch(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.05,-0.1,['r = ' sprintf('%2.2f',rall)])
text(0.05,-0.2,['RMSE = ' sprintf('%2.2f',rmse)])
text(0.05,-0.3,['Fmed = ' sprintf('%2.2f',Fall)])
%axis([-0.1 1.1 -0.1 1.1])
axis([0 1 0 1])
axis equal
xlabel('DvD model predict')
ylabel('POEM Climatology')
title('Fraction Large Pelagics All LME')

subplot(1,2,2)
plot(x,x,'--k');hold on;
for i=1:length(notLELC)
    lme=notLELC(i);
    plot(mod_frac(lme),plme_rPDcatch(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.05,-0.1,['r = ' sprintf('%2.2f',rall2)])
text(0.05,-0.2,['RMSE = ' sprintf('%2.2f',rmse2)])
text(0.05,-0.3,['Fmed = ' sprintf('%2.2f',Fall2)])
%axis([-0.1 1.1 -0.1 1.1])
axis([0 1 0 1])
axis equal
xlabel('DvD model predict')
ylabel('POEM Climatology')
title('Fraction Large Pelagics no LELC')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_MEOW_DvDmod_comp_LELC_subplot.png'])

%Daniel's catch data
figure(2)
subplot(1,2,1)
plot(x,x,'--k');hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(obs_frac(lme),plme_rPDcatch(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.05,-0.1,['r = ' sprintf('%2.2f',robs)])
text(0.05,-0.2,['RMSE = ' sprintf('%2.2f',rmse_obs)])
text(0.05,-0.3,['Fmed = ' sprintf('%2.2f',Fobs)])
%axis([-0.1 1.1 -0.1 1.1])
axis([0 1 0 1])
axis equal
xlabel('DvD catch data')
ylabel('POEM Climatology')
title('Fraction Large Pelagics All LME')

subplot(1,2,2)
plot(x,x,'--k');hold on;
for i=1:length(notLELC)
    lme=notLELC(i);
    plot(obs_frac(lme),plme_rPDcatch(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.05,-0.1,['r = ' sprintf('%2.2f',robs2)])
text(0.05,-0.2,['RMSE = ' sprintf('%2.2f',rmse_obs2)])
text(0.05,-0.3,['Fmed = ' sprintf('%2.2f',Fobs2)])
%axis([-0.1 1.1 -0.1 1.1])
axis([0 1 0 1])
axis equal
xlabel('DvD catch data')
ylabel('POEM Climatology')
title('Fraction Large Pelagics no LELC')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_MEOW_DvDobs_comp_LELC_subplot.png'])


%% Catch
diffP = rPD_catch - dlme_Pfrac;
diffS = rPD_catch - dlme_Sfrac;

% All 4 on subplots
figure(3)
% POEM
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,rPD_catch)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology Fraction LP catch')

% Diff
subplot('Position',[0.25 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('POEM - DvD difference')

% DvD
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,dlme_Pfrac)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('DvD Model Fraction LP catch')
stamp(cfile)
print('-dpng',[ppath 'Clim_' harv '_MEOW_LME_fracPD_catch_modDvD_comp.png'])

% All 4 on subplots
figure(4)
% POEM
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,rPD_catch)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology Fraction LP catch')

% Diff
subplot('Position',[0.25 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffS)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('POEM - DvD difference')

% DvD
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,dlme_Sfrac)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('DvD Catch Obs Fraction LP catch')
stamp(cfile)
print('-dpng',[ppath 'Clim_' harv '_MEOW_LME_fracPD_catch_obsDvD_comp.png'])


