% Visualize difference between
% ESM2M Hindcast of 1951-2000 
% and Forecast of 2051-2100
% production of phyto, zoop, and fish

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% NPP and zoop
load([bpath 'cobalt_npp_means.mat']);
load([bpath 'cobalt_zoop_biom_means.mat']);
load([bpath 'cobalt_det_biom_means.mat']);

% molN/m2/s --> g/m2/d
mzloss_hist = mzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_hist = lzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_hist = npp_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_hist = det_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

mzloss_fore = mzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_fore = lzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_fore = npp_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_fore = det_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

zloss_hist = mzloss_hist + lzloss_hist;
zloss_fore = mzloss_fore + lzloss_fore;

%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

% Hindcast
load([fpath 'Means_Historic_' harv '_prod_' cfile '.mat'],...
    'sf_prod50','sp_prod50','sd_prod50',...
    'mf_prod50','mp_prod50','md_prod50',...
    'lp_prod50','ld_prod50','b_mean50');

[hi,hj]=size(geolon_t);
Hsf=NaN*ones(hi,hj);
Hsp=NaN*ones(hi,hj);
Hsd=NaN*ones(hi,hj);
Hmf=NaN*ones(hi,hj);
Hmp=NaN*ones(hi,hj);
Hmd=NaN*ones(hi,hj);
Hlp=NaN*ones(hi,hj);
Hld=NaN*ones(hi,hj);
Hb =NaN*ones(hi,hj);

Hsf(grid(:,1))=sf_prod50;
Hsp(grid(:,1))=sp_prod50;
Hsd(grid(:,1))=sd_prod50;
Hmf(grid(:,1))=mf_prod50;
Hmp(grid(:,1))=mp_prod50;
Hmd(grid(:,1))=md_prod50;
Hlp(grid(:,1))=lp_prod50;
Hld(grid(:,1))=ld_prod50;
Hb(grid(:,1)) =b_mean50;

Hsf(Hsf(:)<0)=0;
Hsp(Hsp(:)<0)=0;
Hsd(Hsd(:)<0)=0;
Hmf(Hmf(:)<0)=0;
Hmp(Hmp(:)<0)=0;
Hmd(Hmd(:)<0)=0;
Hlp(Hlp(:)<0)=0;
Hld(Hld(:)<0)=0;

clear sf_prod50 sp_prod50 sd_prod50 mf_prod50 mp_prod50 md_prod50 lp_prod50 ld_prod50 b_mean50

%% Forecast
load([fpath 'Means_fore_',harv,'_' cfile '.mat'],...
    'sf_prod50','sp_prod50','sd_prod50',...
    'mf_prod50','mp_prod50','md_prod50',...
    'lp_prod50','ld_prod50','b_mean50');

[ni,nj]=size(geolon_t);
Csf=NaN*ones(ni,nj);
Csp=NaN*ones(ni,nj);
Csd=NaN*ones(ni,nj);
Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);
Cb =NaN*ones(ni,nj);

Csf(ID)=sf_prod50;
Csp(ID)=sp_prod50;
Csd(ID)=sd_prod50;
Cmf(ID)=mf_prod50;
Cmp(ID)=mp_prod50;
Cmd(ID)=md_prod50;
Clp(ID)=lp_prod50;
Cld(ID)=ld_prod50;
Cb(ID) =b_mean50;

Csf(Csf(:)<0)=0;
Csp(Csp(:)<0)=0;
Csd(Csd(:)<0)=0;
Cmf(Cmf(:)<0)=0;
Cmp(Cmp(:)<0)=0;
Cmd(Cmd(:)<0)=0;
Clp(Clp(:)<0)=0;
Cld(Cld(:)<0)=0;

clear sf_prod50 sp_prod50 sd_prod50 mf_prod50 mp_prod50 md_prod50 lp_prod50 ld_prod50 b_mean50

%% Groups
cF = Csf+Cmf;
cP = Csp+Cmp+Clp;
cD = Csd+Cmd+Cld;
cS = Csp+Csf+Csd;
cM = Cmp+Cmf+Cmd;
cL = Clp+Cld;
cAll = cF+cP+cD;

hF = Hsf+Hmf;
hP = Hsp+Hmp+Hlp;
hD = Hsd+Hmd+Hld;
hS = Hsp+Hsf+Hsd;
hM = Hmp+Hmf+Hmd;
hL = Hlp+Hld;
hAll = hF+hP+hD;

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

%% Percent differences
diffN = (npp_fore-npp_hist) ./ npp_hist;
diffMZ = (mzloss_fore-mzloss_hist) ./ mzloss_hist;
diffLZ = (lzloss_fore-lzloss_hist) ./ lzloss_hist;
diffZ = (zloss_fore-zloss_hist) ./ zloss_hist;
diffDet = (det_fore-det_hist) ./ det_hist;

diffL = (cL-hL) ./ hL;
diffM = (cM-hM) ./ hM;
diffF = (cF-hF) ./ hF;
diffP = (cP-hP) ./ hP;
diffD = (cD-hD) ./ hD;
diffB = (Cb-Hb) ./ Hb;
diffAll = (cAll-hAll) ./ hAll;


% diffF(hF(:)<1e-6) = nan;
% diffP(hP(:)<1e-6) = nan;
% diffD(hD(:)<1e-6) = nan;
% diffB(Hb(:)<1e-6) = nan;
% diffAll(hAll(:)<1e-6) = nan;

%% Maps
%F
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast F');
print('-dpng',[pp 'Hist_Fore_' harv '_global_prod_diffF.png'])

% P
figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast P');
print('-dpng',[pp 'Hist_Fore_' harv '_global_prod_diffP.png'])

figure(3)
%D
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast D');
print('-dpng',[pp 'Hist_Fore_' harv '_global_prod_diffD.png'])

%diff all
figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast All');
print('-dpng',[pp 'Hist_Fore_' harv '_global_prod_diffAll.png'])

%% All 4 on subplots
figure(5)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffF))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('F')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffD))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('D')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffP))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('P')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffAll))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('All fishes')
%stamp(cfile)
print('-dpng',[pp 'Hist_Fore_',harv,'_global_prod_diff_subplot.png'])

%% Sizes
%M
figure(6)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffM)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast Medium');

%L
subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffL)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast Large');
print('-dpng',[pp 'Hist_Fore_' harv '_global_prod_diff_size.png'])

%% All 4 on subplots
figure(7)
% NPP
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffN))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Forecast - Hindcast NPP')

% Zoop loss
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffZ))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Forecast - Hindcast Z')

% Medium
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffM))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Forecast - Hindcast M Fish')

% Large
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffL))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Forecast - Hindcast L Fish')
%stamp(cfile)
print('-dpng',[pp 'Hist_Fore_',harv,'_global_prod_diff_plankton_fish_subplot.png'])

%% All 4 on subplots
figure(16)
% NPP
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffN))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Forecast - Hindcast NPP')

% Det
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffDet))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Forecast - Hindcast Detritus')

% Zoop
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffZ))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Forecast - Hindcast Z')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffAll))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Forecast - Hindcast Fish')
%stamp(cfile)
print('-dpng',[pp 'Hist_Fore_',harv,'_global_prod_diff_food_fish_subplot.png'])




