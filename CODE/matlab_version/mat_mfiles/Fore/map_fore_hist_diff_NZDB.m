% Visualize difference between
% ESM2M Hindcast of 1951-2000 
% and Forecast of 2051-2100

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];

% colors
cmG=cbrewer('seq','Greens',50,'PCHIP');
cmB=cbrewer('seq','Blues',50,'PCHIP');
cmP=cbrewer('seq','Purples',50,'PCHIP');
cmR=cbrewer('seq','Reds',50,'PCHIP');
cmBW=cbrewer('seq','Greys',50,'PCHIP');

%% NPP and zoop
load([bpath 'cobalt_det_temp_zoop_npp_means.mat']);
load([fpath 'ts_Hist_Fore_Zp_D_ZpDet_intA.mat']);

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

%time series of 5-yr means as difference from 1951
dtD = tD - tD(19);
dtZ = tZ - tZ(19);
dtZD = tZD - tZD(19);

pdtD = (tD - tD(19)) / tD(19);
pdtZ = (tZ - tZ(19)) / tZ(19);
pdtZD = (tZD - tZD(19)) / tZD(19);

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

%% ts colors
cm=[0.4 0.2 0;... %brown
    0 0 0;...      %black
    0.5 0.5 0.5;...    %med grey
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255;... %peach
    1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.5 0;...    %dk green
    0 0.7 0;...   %g
    127/255 255/255 0;... %lime green
    0 1 1;...     %c
    0 0 0.75;...
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...
    0.5 0 1];%purple

cmap = [...
     0.57255      0.58431      0.56863   %grey
           1       0.8431            0   %yellow
     0.97647         0.19            0   %red
           0            0         0.65   %blue 
         0.4          0.2            0   %brown
         0.1         0.65      0.10196]; %green 

set(groot,'defaultAxesColorOrder',cm);

%% Individual Zprod, Det w/ ts as subplot diffs
figure(3)
subplot('Position',[0 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffZ)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar('Position',[0.05 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Percent change ZP');

subplot('Position',[0 0.05 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffDet)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
title('Percent change Det');

% Zp:Det diff
subplot('Position',[0.5 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffZD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-20 20]);
colorbar('Position',[0.55 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Forecast - Hindcast ZP:Det');

%ts 
subplot('Position',[0.575 0.075 0.35 0.4])
yyaxis left
line(y(19:end),dtD(19:end),'Linewidth',2); hold on;
line(y(19:end),dtZ(19:end),'color',[1 0.5 0],'Linewidth',2); hold on;
ylabel('Change in biomass');

yyaxis right
line(y(19:end),dtZD(19:end),'Linewidth',2); hold on;
xlabel('Year')
ylabel('Change in Z:D');
print('-dpng',[pp 'Hist_Fore_' harv '_global_pdiff_ZpDet_4plot_v1.png'])

%% Individual Zprod, Det w/ ts as subplot pdiff Z, D
figure(4)
subplot('Position',[0 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffZ)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar('Position',[0.05 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Percent change ZP');

subplot('Position',[0 0.05 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffDet)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
title('Percent change Det');

% Zp:Det diff
subplot('Position',[0.5 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffZD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-20 20]);
colorbar('Position',[0.55 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Forecast - Hindcast ZP:Det');

%ts 
subplot('Position',[0.575 0.075 0.35 0.4])
yyaxis left
line(y(19:end),pdtD(19:end),'Linewidth',2); hold on;
line(y(19:end),pdtZ(19:end),'color',[1 0.5 0],'Linewidth',2); hold on;
ylabel('Percent change in biomass');

yyaxis right
line(y(19:end),dtZD(19:end),'Linewidth',2); hold on;
xlabel('Year')
ylabel('Change in Z:D');
print('-dpng',[pp 'Hist_Fore_' harv '_global_pdiff_ZpDet_4plot_v2.png'])


%% Individual Zprod, Det w/ ts as subplot pdiff Z, D
figure(5)
subplot('Position',[0 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffZ)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar('Position',[0.05 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Percent change ZP');

subplot('Position',[0 0.05 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffDet)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
title('Percent change Det');

% Zp:Det diff
subplot('Position',[0.5 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffZD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-20 20]);
colorbar('Position',[0.55 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Forecast - Hindcast ZP:Det');

%ts 
subplot('Position',[0.575 0.075 0.35 0.4])
yyaxis left
line(y(19:end),pdtD(19:end),'Linewidth',2); hold on;
line(y(19:end),pdtZ(19:end),'color',[1 0.5 0],'Linewidth',2); hold on;
ylabel('Percent change in biomass');

yyaxis right
line(y(19:end),pdtZD(19:end),'Linewidth',2); hold on;
xlabel('Year')
ylabel('Percent change in Z:D');
print('-dpng',[pp 'Hist_Fore_' harv '_global_pdiff_ZpDet_4plot_v3.png'])


%% Individual Zprod, Det w/ ts as subplot pdiff Z, D
% USE RED AND BLUE
cmap = [...
     0            0         0.65   %blue
           0       0            0   %black
     0.97647         0.19            0]; %red

set(groot,'defaultAxesColorOrder',cmap);

figure(6)
subplot('Position',[0 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffZ)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar('Position',[0.05 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Percent change Zoo');
text(-2.75,1.75,'A')

subplot('Position',[0 0.05 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffDet)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
title('Percent change Det');
text(-2.75,2,'B')

% Zp:Det diff
subplot('Position',[0.5 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffZD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-20 20]);
colorbar('Position',[0.55 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Change in Zoo:Det');
text(-2.75,1.75,'C')

%ts 
subplot('Position',[0.575 0.075 0.35 0.4])
yyaxis left
line(y(19:end),100*pdtD(19:end),'Linewidth',2); hold on;
line(y(19:end),100*pdtZ(19:end),'color',[0.65 0 0],'Linewidth',2); hold on;
ylabel('Percent change in production');

yyaxis right
line(y(19:end),dtZD(19:end),'Linewidth',2); hold on;
xlabel('Year')
ylabel('Change in Zoo:Det');
text(1925,0.4,'D')
print('-dpng',[pp 'Hist_Fore_' harv '_global_pdiff_ZpDet_4plot_v4.png'])

