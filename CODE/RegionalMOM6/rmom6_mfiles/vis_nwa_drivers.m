% View time series and maps
% Reanalysis-forced runs 1993-2023, v20250715
% MOM6-NWA12

clear
close all

cpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-NWA12/';
ppath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/NWA12/';

%% Units
%poc flux: mol N m-2 s-1
%zoo: mol N m-2
%tp: degC
%tb: degC

load([cpath 'temp_100.nwa.full.hcast.monthly.raw.r20250715.199301-202312.mat'],'temp_100');
load([cpath 'btm_temp.nwa.full.hcast.monthly.raw.r20250715.199301-202312.mat'],'btm_temp');
load([cpath 'fntot_btm.nwa.full.hcast.monthly.raw.r20250715.199301-202312.mat'],'fntot_btm','fntot_btm_long_name','fntot_btm_units')
load([cpath 'nmdz_nlgz_100.nwa.full.hcast.monthly.raw.r20250715.199301-202312.mat']);
load([cpath 'jhploss_mdz_lgz_100.nwa.full.hcast.monthly.raw.r20250715.199301-202312.mat']);

%% nans
temp_100(temp_100 > 1.0e19) = nan;
btm_temp(btm_temp > 1.0e19) = nan;
fntot_btm(fntot_btm > 1.0e19) = nan;
nmdz_100(nmdz_100 > 1.0e19) = nan;
nlgz_100(nlgz_100 > 1.0e19) = nan;
jhploss_nmdz_100(jhploss_nmdz_100 > 1.0e19) = nan;
jhploss_nlgz_100(jhploss_nlgz_100 > 1.0e19) = nan;

%% units
nmdz_100 = nmdz_100 * (106/16) * 12.01  * 9.0;
nlgz_100 = nlgz_100 * (106/16) * 12.01  * 9.0;
jhploss_nmdz_100 = jhploss_nmdz_100 * (106/16) * 12.01 * 9.0 * 60 * 60 * 24;
jhploss_nlgz_100 = jhploss_nlgz_100 * (106/16) * 12.01 * 9.0 * 60 * 60 * 24;
fntot_btm = fntot_btm * (106/16) * 12.01 * 9.0 * 60 * 60 * 24;

%%
load([cpath 'nwa_raw_ocean_static_gridspec.mat'],'geolon','geolat');
load([cpath 'Data_grid_mom6_nwa12.mat'],'GRD');

[ni,nj]=size(geolon);
geolon = double(geolon);
geolat = double(geolat);
%NWAtl
plotminlat=5; 
plotmaxlat=60;
plotminlon=-100;
plotmaxlon=-30;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;   %decent looking coastlines

%% colors
cm10=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0.5 0.5 0.5; ...    %med grey
    0 0 0];...      %black
    

set(groot,'defaultAxesColorOrder',cm10);

%% time
%units = 'days since 1993-01-01 00:00:00'
mos = length(time);
nyrs = mos/12;

yrs = floor(yr(1)):yr(end);

%% annual means

st=1:12:mos;
en=12:12:mos;

tp_amean = nan*ones(nyrs,1);
tb_amean = tp_amean;
mz_amean = tp_amean;
lz_amean = tp_amean;
mzhp_amean = tp_amean;
lzhp_amean = tp_amean;
det_amean = tp_amean;

for n=1:nyrs
    % mean 
    tp_amean(n) = squeeze(mean(temp_100(:,:,st(n):en(n)),[1,2,3],'omitnan'));
    tb_amean(n) = squeeze(mean(btm_temp(:,:,st(n):en(n)),[1,2,3],'omitnan'));
    mz_amean(n) = squeeze(mean(nmdz_100(:,:,st(n):en(n)),[1,2,3],'omitnan'));
    lz_amean(n) = squeeze(mean(nlgz_100(:,:,st(n):en(n)),[1,2,3],'omitnan'));
    mzhp_amean(n) = squeeze(mean(jhploss_nmdz_100(:,:,st(n):en(n)),[1,2,3],'omitnan'));
    lzhp_amean(n) = squeeze(mean(jhploss_nlgz_100(:,:,st(n):en(n)),[1,2,3],'omitnan'));
    det_amean(n) = squeeze(mean(fntot_btm(:,:,st(n):en(n)),[1,2,3],'omitnan'));

end

%% seasonal cycle
tp_smean = nan*ones(12,1);
tb_smean = tp_smean;
mz_smean = tp_smean;
lz_smean = tp_smean;
mzhp_smean = tp_smean;
lzhp_smean = tp_smean;
det_smean = tp_smean;

for n=1:12
    s = n:12:mos;
    % mean 
    tp_smean(n)=mean(temp_100(:,:,s),[1,2,3],'omitnan');
    tb_smean(n)=mean(btm_temp(:,:,st(n):en(n)),[1,2,3],'omitnan');
    mz_smean(n)=mean(nmdz_100(:,:,st(n):en(n)),[1,2,3],'omitnan');
    lz_smean(n)=mean(nlgz_100(:,:,st(n):en(n)),[1,2,3],'omitnan');
    mzhp_smean(n)=mean(jhploss_nmdz_100(:,:,st(n):en(n)),[1,2,3],'omitnan');
    lzhp_smean(n)=mean(jhploss_nlgz_100(:,:,st(n):en(n)),[1,2,3],'omitnan');
    det_smean(n)=mean(fntot_btm(:,:,st(n):en(n)),[1,2,3],'omitnan');

end

%% spatial means
tp_mmean = squeeze(mean(temp_100,3,'omitnan'));
tb_mmean = squeeze(mean(btm_temp,3,'omitnan'));
mz_mmean = squeeze(mean(nmdz_100,3,'omitnan'));
lz_mmean = squeeze(mean(nlgz_100,3,'omitnan'));
mzhp_mmean = squeeze(mean(jhploss_nmdz_100,3,'omitnan'));
lzhp_mmean = squeeze(mean(jhploss_nlgz_100,3,'omitnan'));
det_mmean = squeeze(mean(fntot_btm,3,'omitnan'));

%% Annual means, full time series
figure(1)
subplot(3,3,1)
plot(yrs,tp_amean,'r','Linewidth',1); hold on;
title('TP')
ylabel('^oC')

subplot(3,3,2)
plot(yrs,tb_amean,'b','Linewidth',1); hold on;
title('TB')
ylabel('^oC')

subplot(3,3,6)
plot(yrs,(det_amean),'k','Linewidth',1); hold on;
title('Det')
ylabel('gWW/m2/d')

subplot(3,3,4)
plot(yrs,(mz_amean),'m','Linewidth',1); hold on;
title('MZ')
ylabel('gWW/m2')

subplot(3,3,5)
plot(yrs,(lz_amean),'color',[0 0.5 0.75],'Linewidth',1); hold on;
title('LZ')
ylabel('gWW/m2')

subplot(3,3,7)
plot(yrs,(mzhp_amean),'m','Linewidth',1); hold on;
xlabel('Year')
title('MZ HPloss')
ylabel('gWW/m2/d')

subplot(3,3,8)
plot(yrs,(lzhp_amean),'color',[0 0.5 0.75],'Linewidth',1); hold on;
xlabel('Year')
title('LZ HPloss')
ylabel('gWW/m2/d')

stamp('')
print('-dpng',[ppath 'nwa.full.hcast.monthly.raw.r20250715_all_drivers_annmean.png'])

%% Seasonal cycle
figure(2)
subplot(3,3,1)
plot(1:12,tp_smean,'r','Linewidth',1); hold on;
title('TP')
ylabel('^oC')

subplot(3,3,2)
plot(1:12,tb_smean,'b','Linewidth',1); hold on;
title('TB')
ylabel('^oC')

subplot(3,3,6)
plot(1:12,(det_smean),'k','Linewidth',1); hold on;
title('Det')
ylabel('gWW/m2/d')

subplot(3,3,4)
plot(1:12,(mz_smean),'m','Linewidth',1); hold on;
title('MZ')
ylabel('gWW/m2')

subplot(3,3,5)
plot(1:12,(lz_smean),'color',[0 0.5 0.75],'Linewidth',1); hold on;
title('LZ')
ylabel('gWW/m2')

subplot(3,3,7)
plot(1:12,(mzhp_smean),'m','Linewidth',1); hold on;
xlabel('Mo')
title('MZ HPloss')
ylabel('gWW/m2/d')

subplot(3,3,8)
plot(1:12,(lzhp_smean),'color',[0 0.5 0.75],'Linewidth',1); hold on;
xlabel('Mo')
title('LZ HPloss')
ylabel('gWW/m2/d')

stamp('')
print('-dpng',[ppath 'nwa.full.hcast.monthly.raw.r20250715_all_drivers_seasonal_cycle.png'])


%% 8plot by fn type and size
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - TP
subplot('Position',[0.015 0.75 0.44 0.23])
axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,tp_mmean)
cmocean('thermal')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 30])
set(gcf,'renderer','painters')
text(0,1.35,'TP','HorizontalAlignment','center')

%B - MZ
subplot('Position',[0.015 0.5 0.44 0.23])
axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,log10(mz_mmean))
cmocean('dense')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2])
set(gcf,'renderer','painters')
text(0,1.35,'log_1_0 MZ','HorizontalAlignment','center')

%C - MZloss
subplot('Position',[0.015 0.25 0.44 0.23])
axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,log10(mzhp_mmean))
cmocean('dense')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3 0])
set(gcf,'renderer','painters')
text(0,1.35,'log_1_0 MZ HPloss','HorizontalAlignment','center')

%D - Det
subplot('Position',[0.015 0.0 0.44 0.23])
axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,log10(det_mmean))
cmocean('dense')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2.5 0.5])
set(gcf,'renderer','painters')
text(0,1.35,'log_1_0 Det','HorizontalAlignment','center')

%E - TB
subplot('Position',[0.47 0.75 0.44 0.23])
axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,tb_mmean)
cmocean('thermal')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 30])
set(gcf,'renderer','painters')
text(0,1.35,'TB','HorizontalAlignment','center')

%F - LZ
subplot('Position',[0.47 0.5 0.44 0.23])
axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,log10(lz_mmean))
cmocean('dense')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2])
set(gcf,'renderer','painters')
text(0,1.35,'log_1_0 LZ','HorizontalAlignment','center')

%G - LZloss
subplot('Position',[0.47 0.25 0.44 0.23])
axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,log10(lzhp_mmean))
cmocean('dense')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3 0])
set(gcf,'renderer','painters')
text(0,1.35,'log_1_0 LZ HPloss','HorizontalAlignment','center')

%H - all
% subplot('Position',[0.47 0.0 0.44 0.25])
% axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
%    'Grid','off','FLineWidth',1)
% surfm(geolat,geolon,log10(Eld))
% cmocean('dense')
% %colorbar
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-2 2])
% set(gcf,'renderer','painters')
% text(0,1.75,'LD','HorizontalAlignment','center')

stamp('')
print('-dpng',[ppath 'nwa.full.hcast.monthly.raw.r20250715_map_means.png'])


%% SAVE
save([cpath 'nwa.full.hcast.monthly.raw.r20250715.199301-202312_time_space_means.mat'],...
    'tp_mmean','tb_mmean','mz_mmean','lz_mmean','mzhp_mmean','lzhp_mmean','det_mmean',...
    'tp_amean','tb_amean','mz_amean','lz_amean','mzhp_amean','lzhp_amean','det_amean',...
    'tp_smean','tb_smean','mz_smean','lz_smean','mzhp_smean','lzhp_smean','det_smean')