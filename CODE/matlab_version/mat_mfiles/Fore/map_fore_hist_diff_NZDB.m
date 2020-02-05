% Visualize difference between
% ESM2M Hindcast of 1951-2000 
% and Forecast of 2051-2100

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

% colors
cmG=cbrewer('seq','Greens',50,'PCHIP');
cmB=cbrewer('seq','Blues',50,'PCHIP');
cmP=cbrewer('seq','Purples',50,'PCHIP');
cmR=cbrewer('seq','Reds',50,'PCHIP');
cmBW=cbrewer('seq','Greys',50,'PCHIP');

%% NPP and zoop
load([bpath 'cobalt_det_temp_zoop_npp_means.mat']);

% molN/m2/s --> g/m2/d
mzloss_hist = mzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_hist = lzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mzprod_hist = mzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_hist = lzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_hist = npp_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_hist = det_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
ptemp_hist = ptemp_mean_hist - 273;

mzloss_fore = mzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_fore = lzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mzprod_fore = mzprod_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_fore = lzprod_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_fore = npp_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_fore = det_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
ptemp_fore = ptemp_mean_fore - 273;

zloss_hist = mzloss_hist + lzloss_hist;
zloss_fore = mzloss_fore + lzloss_fore;

zprod_hist = mzprod_hist + lzprod_hist;
zprod_fore = mzprod_fore + lzprod_fore;

l10ZlDet_hist = log10(zloss_hist./det_hist);
l10ZlDet_fore = log10(zloss_fore./det_fore);

l10ZpDet_hist = log10(zprod_hist./det_hist);
l10ZpDet_fore = log10(zprod_fore./det_fore);

ZlDet_hist = (zloss_hist./det_hist);
ZlDet_fore = (zloss_fore./det_fore);

ZpDet_hist = (zprod_hist./det_hist);
ZpDet_fore = (zprod_fore./det_fore);

%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

% Hindcast
load([fpath 'Means_Historic_' harv '_' cfile '.mat'],...
    'b_mean50');

[hi,hj]=size(geolon_t);
Hb =NaN*ones(hi,hj);
Hb(grid(:,1)) =b_mean50;

clear b_mean50


% Forecast
load([fpath 'Means_fore_',harv,'_' cfile '.mat'],...
    'b_mean50');

[ni,nj]=size(geolon_t);
Cb =NaN*ones(ni,nj);
Cb(ID) =b_mean50;

clear b_mean50


%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

%%
diffZD = (ZpDet_fore-ZpDet_hist);

pdiffN = (npp_fore-npp_hist) ./ npp_hist;
pdiffDet = (det_fore-det_hist) ./ det_hist;
pdiffMZ = (mzprod_fore-mzprod_hist) ./ mzprod_hist;
pdiffLZ = (lzprod_fore-lzprod_hist) ./ lzprod_hist;
pdiffZ = (zprod_fore-zprod_hist) ./ zprod_hist;
pdiffZD = (l10ZpDet_fore-l10ZpDet_hist) ./ l10ZpDet_hist;
pdiffB = (Cb-Hb) ./ Hb;

% diffB(Hb(:)<1e-6) = nan;

%% Maps
% Individual Zprod:Det
figure(1)
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,ZpDet_hist)
colormap(cmP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 20]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast ZP:Det');

subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,ZpDet_fore)
colormap(cmP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 20]);
colorbar
set(gcf,'renderer','painters')
title('Forecast ZP:Det');

% Zl:Det diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffZD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 10]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast ZP:Det');
print('-dpng',[pp 'Hist_Fore_' harv '_global_diffZpDet_3plot.png'])

%% Individual Zprod, Det
figure(2)
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffZ)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('Percent change ZP');

subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffDet)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('Percent change Det');

% Zl:Det diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffZD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 10]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast ZP:Det');
print('-dpng',[pp 'Hist_Fore_' harv '_global_pdiff_ZpDet_3plot.png'])
