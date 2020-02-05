% Visualize difference between
% ESM2M Hindcast of 1951-2000 
% and Forecast of 2051-2100

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

% colors
cmPG=cbrewer('div','PiYG',50,'PCHIP');
cmPO=cbrewer('div','PuOr',50,'PCHIP');
cmG=cbrewer('seq','Greens',50,'PCHIP');
cmB=cbrewer('seq','Blues',50,'PCHIP');
cmP=cbrewer('seq','Purples',50,'PCHIP');
cmR=cbrewer('seq','Reds',50,'PCHIP');
cmBW=cbrewer('seq','Greys',50,'PCHIP');

cmapD = cmR(35,:);
cmapD(2,:) = cmB(35,:);
cmapD(3,:) = cmG(35,:);

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
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

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
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

%%
diffZD = (ZlDet_fore-ZlDet_hist);

cAll = cF+cP+cD;
cFracF = cF ./ cAll;
cFracP = cP ./ cAll;
cFracD = cD ./ cAll;
cFracPD = cP ./ (cP+cD);
cFracPF = cP ./ (cP+cF);
cFracDF = cD ./ (cD+cF);
cFracFD = cF ./ (cF+cD);
cFracFP = cF ./ (cF+cP);
cFracLM = cL ./ (cL+cM);
cFracML = cM ./ (cM+cL);

hAll = hF+hP+hD;
hFracF = hF ./ cAll;
hFracP = hP ./ cAll;
hFracD = hD ./ cAll;
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
pdiffPf = (cFracP-hFracP) ./ hFracP;
pdiffFf = (cFracF-hFracF) ./ hFracF;
pdiffDf = (cFracD-hFracD) ./ hFracD;
pdiffPD = (cFracPD-hFracPD) ./ hFracPD;
pdiffPF = (cFracPF-hFracPF) ./ hFracPF;
pdiffDF = (cFracDF-hFracDF) ./ hFracDF;
pdiffFD = (cFracFD-hFracFD) ./ hFracFD;
pdiffFP = (cFracFP-hFracFP) ./ hFracFP;
pdiffLM = (cFracLM-hFracLM) ./ hFracLM;

diffPf = (cFracP-hFracP);
diffFf = (cFracF-hFracF);
diffDf = (cFracD-hFracD);
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

%% Dominant type
hdF = find(hFracF >= 0.5);
hdP = find(hFracP >= 0.5);
hdD = find(hFracD >= 0.5);

cdF = find(cFracF >= 0.5);
cdP = find(cFracP >= 0.5);
cdD = find(cFracD >= 0.5);

hdom = nan(ni,nj);
hdom(hdF) = 1;
hdom(hdP) = 2;
hdom(hdD) = 3;

cdom = nan(ni,nj);
cdom(cdF) = 1;
cdom(cdP) = 2;
cdom(cdD) = 3;

hdomF = zeros(ni,nj);
hdomF(hdF) = 1;
cdomF = zeros(ni,nj);
cdomF(cdF) = 1;

hdomP = zeros(ni,nj);
hdomP(hdP) = 1;
cdomP = zeros(ni,nj);
cdomP(cdP) = 1;

hdomD = zeros(ni,nj);
hdomD(hdD) = 1;
cdomD = zeros(ni,nj);
cdomD(cdD) = 1;

%% Maps
% Individual Zl:Det
figure(1)
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,ZlDet_hist)
colormap(cmP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 20]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast Zl:Det');

subplot('Position',[0.5 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,ZlDet_fore)
colormap(cmP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 20]);
colorbar
set(gcf,'renderer','painters')
title('Forecast Zl:Det');

% Zl:Det diff
subplot('Position',[0.25 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffZD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 10]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast Zl:Det');
print('-dpng',[pp 'Hist_Fore_' harv '_global_diffZlDet_3plot.png'])

%% Individual Zl:Det
figure(2)
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,ZlDet_hist)
colormap(cmP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 20]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast Zl:Det');

subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,ZlDet_fore)
colormap(cmP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 20]);
colorbar
set(gcf,'renderer','painters')
title('Forecast Zl:Det');

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
title('Forecast - Hindcast Zl:Det');
print('-dpng',[pp 'Hist_Fore_' harv '_global_diffZlDet_3plot_v2.png'])

%% Frac F
figure(3)
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,hFracF)
colormap(cmR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast F / All');

subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,cFracF)
colormap(cmR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
title('Forecast F / All');

% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffFf)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast F / All');
print('-dpng',[pp 'Hist_Fore_' harv '_global_fracF_cmR.png'])


%% Frac P
figure(4)
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,hFracP)
colormap(cmB)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast P / All');

subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,cFracP)
colormap(cmB)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
title('Forecast P / All');

% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffPf)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast P / All');
print('-dpng',[pp 'Hist_Fore_',harv,'_global_fracP_cmB.png'])

%% Frac D
figure(5)
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,hFracD)
colormap(cmG)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast D / All');

subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,cFracD)
colormap(cmG)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
title('Forecast D / All');

% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffDf)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast D / All');
print('-dpng',[pp 'Hist_Fore_',harv,'_global_fracD_cmG.png'])


%% Dominance
figure(6)
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,hdom)
colormap(cmapD)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 3]);
%colorbar('southoutside','Ticks',[1 2 3],'TickLabels',{'F','P','D'})
colorbar('Position',[0.05 0.53 0.4 0.025],'orientation','horizontal','Ticks',[1 2 3],'TickLabels',{'F','P','D'})
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
title('Hindcast');

subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,cdom)
colormap(cmapD)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 3]);
% colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
title('Forecast');

% Diff
subplot('Position',[0.5 0.66 0.48 0.33])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,cdomF-hdomF)
cmocean('balance')
%colormap(cmR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.25 1.25]);
%colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
text(-0.1,1.5,'Forage fish','FontWeight','bold','HorizontalAlignment','center');

subplot('Position',[0.5 0.33 0.475 0.33])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,cdomP-hdomP)
cmocean('balance')
%colormap(cmR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.25 1.25]);
colorbar('Position',[0.95 0.34 0.02 0.33],'orientation','vertical','Ticks',[-1 0 1])
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
text(-0.1,1.5,'Large pelagics','FontWeight','bold','HorizontalAlignment','center');

subplot('Position',[0.5 0.0 0.475 0.33])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,cdomD-hdomD)
cmocean('balance')
%colormap(cmR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.25 1.25]);
%colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
text(-0.1,1.5,'Demersals','FontWeight','bold','HorizontalAlignment','center');
print('-dpng',[pp 'Hist_Fore_',harv,'_global_dominance.png'])

