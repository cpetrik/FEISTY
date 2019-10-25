% Visualize difference between
% ESM2M Hindcast of 1951-2000 
% and Forecast of 2051-2100

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

% colors
cmPG=cbrewer('div','PRGn',50,'PCHIP');
cmPG = flipud(cmPG);
cmPO=cbrewer('div','PuOr',50,'PCHIP');
cmBG=cbrewer('div','BrBG',50,'PCHIP');

cmPBG=cbrewer('seq','PuBuGn',50,'PCHIP');
cmBW=cbrewer('seq','Greys',50,'PCHIP');

%% NPP and zoop
load([bpath 'cobalt_det_temp_zoop_npp_means.mat']);

% molN/m2/s --> g/m2/d
mzloss_hist = mzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_hist = lzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_hist = npp_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_hist = det_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
ptemp_hist = ptemp_mean_hist - 273;

mzloss_fore = mzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_fore = lzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_fore = npp_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_fore = det_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
ptemp_fore = ptemp_mean_fore - 273;

zloss_hist = mzloss_hist + lzloss_hist;
zloss_fore = mzloss_fore + lzloss_fore;

l10ZlDet_hist = log10(zloss_hist./det_hist);
l10ZlDet_fore = log10(zloss_fore./det_fore);

ZlDet_hist = (zloss_hist./det_hist);
ZlDet_fore = (zloss_fore./det_fore);

%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/US/'];

harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

% Hindcast
load([fpath 'Means_Historic_' harv '_' cfile '.mat'],...
    'sf_mean50','sp_mean50','sd_mean50',...
    'mf_mean50','mp_mean50','md_mean50',...
    'lp_mean50','ld_mean50','b_mean50');

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
Hsf(grid(:,1))=sf_mean50;
Hsp(grid(:,1))=sp_mean50;
Hsd(grid(:,1))=sd_mean50;
Hmf(grid(:,1))=mf_mean50;
Hmp(grid(:,1))=mp_mean50;
Hmd(grid(:,1))=md_mean50;
Hlp(grid(:,1))=lp_mean50;
Hld(grid(:,1))=ld_mean50;
Hb(grid(:,1)) =b_mean50;

clear sf_mean50 sp_mean50 sd_mean50 mf_mean50 mp_mean50 md_mean50 lp_mean50 ld_mean50 b_mean50

% Forecast
load([fpath 'Means_fore_',harv,'_' cfile '.mat'],...
    'sf_mean50','sp_mean50','sd_mean50',...
    'mf_mean50','mp_mean50','md_mean50',...
    'lp_mean50','ld_mean50','b_mean50');

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
Csf(ID)=sf_mean50;
Csp(ID)=sp_mean50;
Csd(ID)=sd_mean50;
Cmf(ID)=mf_mean50;
Cmp(ID)=mp_mean50;
Cmd(ID)=md_mean50;
Clp(ID)=lp_mean50;
Cld(ID)=ld_mean50;
Cb(ID) =b_mean50;

clear sf_mean50 sp_mean50 sd_mean50 mf_mean50 mp_mean50 md_mean50 lp_mean50 ld_mean50 b_mean50

cF = Csf+Cmf;
cP = Csp+Cmp+Clp;
cD = Csd+Cmd+Cld;
cS = Csp+Csf+Csd;
cM = Cmp+Cmf+Cmd;
cL = Clp+Cld;

hF = Hsf+Hmf;
hP = Hsp+Hmp+Hlp;
hD = Hsd+Hmd+Hld;
hS = Hsp+Hsf+Hsd;
hM = Hmp+Hmf+Hmd;
hL = Hlp+Hld;

%% plot info
plotminlat=10; %Set these bounds for your data
plotmaxlat=70;
plotminlon=-180;
plotmaxlon=-50;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

%%
diffZD = (ZlDet_fore-ZlDet_hist);

cAll = cF+cP+cD;
cFracPD = cP ./ (cP+cD);
cFracPF = cP ./ (cP+cF);
cFracDF = cD ./ (cD+cF);
cFracFD = cF ./ (cF+cD);
cFracFP = cF ./ (cF+cP);
cFracLM = cL ./ (cL+cM);
cFracML = cM ./ (cM+cL);

hAll = hF+hP+hD;
hFracPD = hP ./ (hP+hD);
hFracPF = hP ./ (hP+hF);
hFracDF = hD ./ (hD+hF);
hFracFD = hF ./ (hF+hD);
hFracFP = hF ./ (hF+hP);
hFracLM = hL ./ (hL+hM);
hFracML = hM ./ (hM+hL);

pdiffN = (npp_fore-npp_hist) ./ npp_hist;
pdiffDet = (det_fore-det_hist) ./ det_hist;
pdiffMZ = (mzloss_fore-mzloss_hist) ./ mzloss_hist;
pdiffLZ = (lzloss_fore-lzloss_hist) ./ lzloss_hist;
pdiffZ = (zloss_fore-zloss_hist) ./ zloss_hist;
pdiffZD = (l10ZlDet_fore-l10ZlDet_hist) ./ l10ZlDet_hist;
pdiffPT = (ptemp_fore-ptemp_hist) ./ ptemp_hist;
dPT = (ptemp_fore-ptemp_hist);

pdiffL = (cL-hL) ./ hL;
pdiffM = (cM-hM) ./ hM;
pdiffF = (cF-hF) ./ hF;
pdiffP = (cP-hP) ./ hP;
pdiffD = (cD-hD) ./ hD;
pdiffB = (Cb-Hb) ./ Hb;
pdiffAll = (cAll-hAll) ./ hAll;
pdiffPD = (cFracPD-hFracPD) ./ hFracPD;
pdiffPF = (cFracPF-hFracPF) ./ hFracPF;
pdiffDF = (cFracDF-hFracDF) ./ hFracDF;
pdiffFD = (cFracFD-hFracFD) ./ hFracFD;
pdiffFP = (cFracFP-hFracFP) ./ hFracFP;
pdiffLM = (cFracLM-hFracLM) ./ hFracLM;

diffPD = (cFracPD-hFracPD);
diffPF = (cFracPF-hFracPF);
diffDF = (cFracDF-hFracDF);
diffFD = (cFracFD-hFracFD);
diffFP = (cFracFP-hFracFP);
diffLM = (cFracLM-hFracLM);
diffML = (cFracML-hFracML);

% diffF(hF(:)<1e-6) = nan;
% diffP(hP(:)<1e-6) = nan;
% diffD(hD(:)<1e-6) = nan;
% diffB(Hb(:)<1e-6) = nan;
% diffAll(hAll(:)<1e-6) = nan;

%% Maps
% Individual Zl:Det
figure(1)
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,ZlDet_hist)
colormap(cmPBG)
% cmocean('dense')
% cmocean('matter')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 15]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast Zl:Det');
set(gca,'XColor','none','YColor','none')

subplot('Position',[0.0 0.03 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,ZlDet_fore)
colormap(cmPBG)
% cmocean('dense')
% cmocean('matter')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 15]);
colorbar
set(gcf,'renderer','painters')
title('Forecast Zl:Det');
set(gca,'XColor','none','YColor','none')

% Zl:Det diff
subplot('Position',[0.5 0.28 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffZD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 10]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast Zl:Det');
set(gca,'XColor','none','YColor','none')
print('-dpng',[pp 'Hist_Fore_' harv '_US_diffZlDet.png'])

%% Individual F:P
figure(2)
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,hFracFP)
colormap(cmPG)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
title('Hindcast F / (F+P)');

subplot('Position',[0.0 0.03 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,cFracFP)
colormap(cmPG)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
title('Forecast F / (F+P)');

% F:P Ratio
subplot('Position',[0.5 0.28 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffFP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
title('Forecast - Hindcast F / (F+P)')
print('-dpng',[pp 'Hist_Fore_',harv,'_US_FP_diff.png'])

%% Individual F:D
figure(3)
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,hFracFD)
colormap(cmPO)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
title('Hindcast F / (F+D)');

subplot('Position',[0.0 0.03 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,cFracFD)
colormap(cmPO)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
title('Forecast F / (F+D)');

% F:D Ratio
subplot('Position',[0.5 0.28 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffFD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
title('Forecast - Hindcast F / (F+D)')
print('-dpng',[pp 'Hist_Fore_',harv,'_US_FD_diff.png'])


%% Individual P:D
figure(4)
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,hFracPD)
colormap(cmBG)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
title('Hindcast P / (P+D)');

subplot('Position',[0.0 0.03 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,cFracPD)
colormap(cmBG)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
title('Forecast P / (P+D)');

% F:D Ratio
subplot('Position',[0.5 0.28 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
title('Forecast - Hindcast P / (P+D)')
print('-dpng',[pp 'Hist_Fore_',harv,'_US_PD_diff.png'])



%% Individual M:L
cmP=cbrewer('div','PiYG',50,'PCHIP');
cmP = flipud(cmP);

figure(5)
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,hFracML)
colormap(cmP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
title('Hindcast M / (M+L)');

subplot('Position',[0.0 0.03 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,cFracML)
colormap(cmP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
title('Forecast M / (M+L)');

% M:L Ratio
subplot('Position',[0.5 0.28 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffML)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast M / (M+L)')
set(gca,'XColor','none','YColor','none')
print('-dpng',[pp 'Hist_Fore_',harv,'_US_ML_diff.png'])



