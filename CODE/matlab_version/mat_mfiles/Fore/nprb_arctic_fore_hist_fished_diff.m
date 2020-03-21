% Maps of FEISTY in Arctic for NPRB proposal

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
epath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';

%% save hist and fore together
load([fpath 'Means_hist_fore_',harv,'_cobalt_' cfile '.mat']);

zprod_hist = mzprod_hist + lzprod_hist;
zprod_fore = mzprod_fore + lzprod_fore;

ZpDet_hist = (zprod_hist./det_hist);
ZpDet_fore = (zprod_fore./det_fore);

%% plot info
plotminlat=55; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

%%
cAll = cF+cP+cD;
cFracD = cD ./ cAll;
cFracDP = cD ./ (cP+cD);
cFracDF = cD ./ (cD+cF);

hAll = hF+hP+hD;
hFracD = hD ./ hAll;
hFracDP = hD ./ (hP+hD);
hFracDF = hD ./ (hD+hF);

pdiffDet = (det_fore-det_hist) ./ det_hist;
pdiffZ = (zprod_fore-zprod_hist) ./ zprod_hist;
pdiffZD = (ZpDet_fore-ZpDet_hist) ./ ZpDet_hist;

pdiffF = (cF-hF) ./ hF;
pdiffP = (cP-hP) ./ hP;
pdiffD = (cD-hD) ./ hD;
pdiffB = (Cb-Hb) ./ Hb;
pdiffAll = (cAll-hAll) ./ hAll;

diffD = (cFracD-hFracD);
diffDP = (cFracDP-hFracDP);
diffDF = (cFracDF-hFracDF);
diffZD = (ZpDet_fore-ZpDet_hist);

%% Maps
figure(1)
% 6 figure subplot Zp:Det, Benthos, D, D:P, D:F, D:all
%1 - Zp:Det
subplot('Position',[0 0.67 0.5 0.275])
axesm ('ortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,(diffZD))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-20 20]);
colorbar
set(gcf,'renderer','painters')
title('Zooplankton to Detritus')


%2 - Benthos
subplot('Position',[0 0.335 0.5 0.275])
axesm ('ortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*pdiffB)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('Benthos');

%3 - D
subplot('Position',[0 0.02 0.5 0.275])
axesm ('ortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*pdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('Demersals');

%4 - DP
subplot('Position',[0.5 0.67 0.5 0.275])
axesm ('ortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,(diffDP))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('Demersals vs. Large Pelagics')

%5 - DF
subplot('Position',[0.5 0.335 0.5 0.275])
axesm ('ortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,(diffDF))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('Demersals vs. Forage')

%6 - D/all
subplot('Position',[0.5 0.02 0.5 0.275])
axesm ('ortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,(diffD))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('D / All')
%stamp(cfile)
print('-dpng',[pp 'Hist_Fore_',harv,'_arctic_6plot_Dfracs_ZpDet.png'])


%% plot info
plotminlat=55; %Set these bounds for your data
plotmaxlat=77.5;
plotminlon=-180;
plotmaxlon=-125;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

figure(2)
% 6 figure subplot Zp:Det, Benthos, D, D:P, D:F, D:all
%1 - Zp:Det
subplot('Position',[0 0.67 0.5 0.275])
axesm ('murdoch1','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,(diffZD))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-20 20]);
colorbar
set(gcf,'renderer','painters')
title('Zooplankton to Detritus')


%2 - Benthos
subplot('Position',[0 0.335 0.5 0.275])
axesm ('murdoch1','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*pdiffB)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('Benthos');

%3 - D
subplot('Position',[0 0.02 0.5 0.275])
axesm ('murdoch1','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*pdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('Demersals');

%4 - DP
subplot('Position',[0.5 0.67 0.5 0.275])
axesm ('murdoch1','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,(diffDP))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('Demersals vs. Large Pelagics')

%5 - DF
subplot('Position',[0.5 0.335 0.5 0.275])
axesm ('murdoch1','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,(diffDF))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('Demersals vs. Forage')

%6 - D/all
subplot('Position',[0.5 0.02 0.5 0.275])
axesm ('murdoch1','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,(diffD))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('D / All')
%stamp(cfile)
print('-dpng',[pp 'Hist_Fore_',harv,'_Chukchi_6plot_Dfracs_ZpDet.png'])
