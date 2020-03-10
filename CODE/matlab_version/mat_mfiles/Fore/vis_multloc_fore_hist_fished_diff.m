% Visualize difference between
% ESM2M Hindcast of 1951-2000 
% and Forecast of 2051-2100

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

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

zprod_hist = mzprod_hist + lzprod_hist;
zprod_fore = mzprod_fore + lzprod_fore;

l10ZpDet_hist = log10(zprod_hist./det_hist);
l10ZpDet_fore = log10(zprod_fore./det_fore);

ZpDet_hist = (zprod_hist./det_hist);
ZpDet_fore = (zprod_fore./det_fore);

%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
%fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
fpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];
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

%% save hist and fore together
save([fpath 'Means_hist_fore_',harv,'_cobalt_' cfile '.mat'],'mzloss_hist',...
    'mzprod_hist','lzloss_hist','lzprod_hist','npp_hist','det_hist',...
    'ptemp_hist','mzloss_fore','mzprod_fore','lzloss_fore','lzprod_fore',...
    'npp_fore','det_fore','ptemp_fore','cF','cP','cD','cS','cM','cL','hF',...
    'hP','hD','hS','hM','hL');

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

%%
diffZD = (ZpDet_fore-ZpDet_hist);

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
pdiffMZ = (mzprod_fore-mzprod_hist) ./ mzprod_hist;
pdiffLZ = (lzprod_fore-lzprod_hist) ./ lzprod_hist;
pdiffZ = (zprod_fore-zprod_hist) ./ zprod_hist;
pdiffZD = (l10ZpDet_fore-l10ZpDet_hist) ./ l10ZpDet_hist;
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
% Individual Hist vs Fore
figure(1)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(hF)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast F');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(cF)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Forecast F');
print('-dpng',[pp 'Hist_Fore_' harv '_global_F.png'])

%P
figure(2)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(hP)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast P');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(cP)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Forecast P');
print('-dpng',[pp 'Hist_Fore_' harv '_global_P.png'])

% D
figure(3)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(hD)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast D');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(cD)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Forecast D');
print('-dpng',[pp 'Hist_Fore_' harv '_global_D.png'])

%4
figure(4)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(hAll)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 2]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast All');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(cAll)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 2]);
colorbar
set(gcf,'renderer','painters')
title('Forecast All');
print('-dpng',[pp 'Hist_Fore_' harv '_global_All.png'])

%% B
figure(5)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(Hb)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast Benthos');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(Cb)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Forecast Benthos');
print('-dpng',[pp 'Hist_Fore_' harv '_global_Bent.png'])

%% Zp:Det
figure(6)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,ZpDet_hist)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 20]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast Zp:Det');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,ZpDet_fore)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 20]);
colorbar
set(gcf,'renderer','painters')
title('Forecast Zp:Det');
print('-dpng',[pp 'Hist_Fore_' harv '_global_ZpDet.png'])

%% pdiffs individual
%F
figure(7)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast F');
print('-dpng',[pp 'Hist_Fore_' harv '_global_pdiffF.png'])

% P
figure(8)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast P');
print('-dpng',[pp 'Hist_Fore_' harv '_global_pdiffP.png'])

figure(9)
%D
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast D');
print('-dpng',[pp 'Hist_Fore_' harv '_global_pdiffD.png'])

%pdiff all
figure(10)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast All');
print('-dpng',[pp 'Hist_Fore_' harv '_global_pdiffAll.png'])

%% B pdiff
figure(11)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffB)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast B');
print('-dpng',[pp 'Hist_Fore_' harv '_global_pdiffB.png'])

%% Zp:Det diff
figure(12)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffZD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-25 25]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast Zp:Det');
print('-dpng',[pp 'Hist_Fore_' harv '_global_diffZpDet.png'])

%% All 4 on subplots 0.5
figure(13)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffF))
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
surfm(geolat_t,geolon_t,(pdiffD))
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
surfm(geolat_t,geolon_t,(pdiffP))
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
surfm(geolat_t,geolon_t,(pdiffAll))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('All fishes')
%stamp(cfile)
print('-dpng',[pp 'Hist_Fore_',harv,'_global_pdiff_subplot50.png'])

%% All 4 on subplots 1
figure(23)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffF))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('F')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffD))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('D')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffP))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('P')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffAll))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('All fishes')
%stamp(cfile)
print('-dpng',[pp 'Hist_Fore_',harv,'_global_pdiff_subplot100.png'])

%% All 4 on subplots 100%
figure(33)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*(pdiffF))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('F')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*(pdiffD))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
title('D')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*(pdiffP))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
title('P')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*(pdiffAll))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
title('All fishes')
%stamp(cfile)
print('-dpng',[pp 'Hist_Fore_',harv,'_global_pdiff_subplot100p.png'])


%% All 4 on subplots 100%
figure(43)
% all M
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*(pdiffM))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Med')

% all L
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*(pdiffL))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
title('Lrg')
%stamp(cfile)
print('-dpng',[pp 'Hist_Fore_',harv,'_global_sizes_pdiff_subplot100p.png'])

%% Ratios on subplots red-white-blue
% 3 figure subplot P:D, P:F, L:M
figure(14)
subplot('Position',[0 0.53 0.5 0.5])
%P:D
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Large Pelagics vs. Demersals')

%P:F
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffPF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Large Pelagics vs. Forage Fishes')

%L:M
subplot('Position',[0.25 0.0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffLM)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Large vs. Medium')
%stamp(cfile)
print('-dpng',[pp 'Hist_Fore_',harv,'_global_ratios_diff_subplot.png'])

%% F Ratios on subplots red-white-blue w/ZpDet
% 3 figure subplot F:D, F:P, Z:Det
figure(24)
subplot('Position',[0 0.53 0.5 0.5])
%F:D
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffFD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Forage Fishes vs. Demersals')

%F:P
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffFP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
colorbar('Position',[0.2 0.56 0.6 0.025],'orientation','horizontal')
title('Forage Fishes vs. Large Pelagics')

%Zoo:Det
subplot('Position',[0.25 0.02 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffZD))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 10]);
set(gcf,'renderer','painters')
title('Zooplankton to Detritus')
colorbar('Position',[0.2 0.05 0.6 0.025],'orientation','horizontal')
print('-dpng',[pp 'Hist_Fore_',harv,'_global_Fratios_diff_subplot_ZpDet.png'])

%% Ratios on subplots with ZpDet
% 4 figure subplot P:D, P:F, L:M, Zp:Det
figure(15)
% PD
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffPD))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Large Pelagics vs. Demersals')

% all ZD
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffZD))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 10]);
set(gcf,'renderer','painters')
title('Zooplankton to Detritus')
colorbar('Position',[0.035 0.05 0.425 0.025],'orientation','horizontal')

% PF
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffPF))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Large Pelagics vs. Forage Fishes')

% LM
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffLM))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Large vs. Medium Fishes')
%stamp(cfile)
print('-dpng',[pp 'Hist_Fore_',harv,'_global_ratios_subplot_ZpDet.png'])

%% F Ratios on subplots with ZpDet
% 4 figure subplot F:D, F:P, M:L, Zp:Det
figure(25)
% FD
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffFD))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Forage Fishes vs.  Demersals')

% all ZD
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffZD))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 10]);
set(gcf,'renderer','painters')
title('Zooplankton to Detritus')
colorbar('Position',[0.035 0.05 0.425 0.025],'orientation','horizontal')

% FP
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffFP))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Forage Fishes vs. Large Pelagics')

% ML
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffML))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Medium vs. Large Fishes')
%stamp(cfile)
print('-dpng',[pp 'Hist_Fore_',harv,'_global_FMLratios_subplot_ZpDet.png'])

%% Ratios on subplots with temp
figure(16)
% all PelT
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(dPT))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
title('Pelagic Temperature')
colorbar('Position',[0.035 0.05 0.425 0.025],'orientation','horizontal')


% FD
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffFD))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Forage Fishes vs.  Demersals')

% FP
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffFP))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Forage Fishes vs. Large Pelagics')

% ML
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffML))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Medium vs. Large Fishes')
%stamp(cfile)
print('-dpng',[pp 'Hist_Fore_',harv,'_global_FMLratios_subplot_PelT.png'])

%% F Ratios on subplots with temp & Zp:Det
figure(26)
% Zp:Det
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffZD))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 10]);
set(gcf,'renderer','painters')
title('Zooplankton to Detritus')
colorbar('Position',[0.035 0.05 0.425 0.025],'orientation','horizontal')

% FD
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffFD))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar('Position',[0.25 0.55 0.5 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Forage Fishes vs.  Demersals')

% FP
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffFP))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Forage Fishes vs. Large Pelagics')

% Pel Temp
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(dPT))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
title('Pelagic Temperature')
colorbar('Position',[0.535 0.05 0.425 0.025],'orientation','horizontal')
print('-dpng',[pp 'Hist_Fore_',harv,'_global_Fratios_subplot_ZpDet_PelT.png'])

%% Size
%M
figure(17)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,pdiffM)
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
surfm(geolat_t,geolon_t,pdiffL)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('Forecast - Hindcast Large');
print('-dpng',[pp 'Hist_Fore_' harv '_global_pdiff_size.png'])

%% All 4 on subplots
figure(18)
% NPP
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffN))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Forecast - Hindcast NPP')

% MZ
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffMZ))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Forecast - Hindcast MZ')

% All Z
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffZ))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Forecast - Hindcast Z')

% LZ
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffLZ))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Forecast - Hindcast LZ')
%stamp(cfile)
print('-dpng',[pp 'Hist_Fore_',harv,'_global_pdiff_plankton_subplot.png'])


%% All 4 on subplots
figure(19)
% NPP
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffN))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Forecast - Hindcast NPP')

% M
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffM))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('Forecast - Hindcast M Fish')

% Z
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffZ))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('Forecast - Hindcast Z')

% L
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffL))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('Forecast - Hindcast L Fish')
%stamp(cfile)
print('-dpng',[pp 'Hist_Fore_',harv,'_global_pdiff_plankton_fish_subplot.png'])

%% All 4 on subplots
figure(20)
% NPP
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffN))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Forecast - Hindcast NPP')

% Det
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffDet))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('Forecast - Hindcast Detritus')

% Zoop
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffZ))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('Forecast - Hindcast Z')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffAll))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('Forecast - Hindcast Fish')
%stamp(cfile)
print('-dpng',[pp 'Hist_Fore_',harv,'_global_pdiff_food_fish_subplot.png'])



