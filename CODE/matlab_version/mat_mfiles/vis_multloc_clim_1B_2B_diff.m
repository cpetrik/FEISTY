% Visualize difference between
% ESM2.6 Climatology of 5 yrs (w/1990 ICs) 
% with 1 or 2 benthos groups

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';

%% Climatology grid
Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']); %depth,ID,lat,lmask,lon
load([cpath 'esm26_area_1deg.mat']); %area
AREA_OCN = max(area,1);

%% FEISTY Output
cfile1 = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath1=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile1 '/'];
cfile2 = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_2B_BE08_noCC_RE00100';
fpath2=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile2 '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile2 '/'];

harv = 'All_fish03';

%% 1B
load([fpath1 '/Climatology/Means_Climatol_' harv '_' cfile1 '.mat'],...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean');

% Put biomass on grid
[ni,nj]=size(lon);
Csf=NaN*ones(ni,nj);
Csp=NaN*ones(ni,nj);
Csd=NaN*ones(ni,nj);
Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);
Cb =NaN*ones(ni,nj);
Csf(ID)=sf_mean;
Csp(ID)=sp_mean;
Csd(ID)=sd_mean;
Cmf(ID)=mf_mean;
Cmp(ID)=mp_mean;
Cmd(ID)=md_mean;
Clp(ID)=lp_mean;
Cld(ID)=ld_mean;
Cb(ID) =b_mean;

clear sf_mean sp_mean sd_mean mf_mean mp_mean md_mean lp_mean ld_mean

%% 2B
load([fpath2 'Means_Climatol_',harv,'_' cfile2 '.mat'],...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','sb_mean','mb_mean');

% Put biomass on grid
Hsf=NaN*ones(ni,nj);
Hsp=NaN*ones(ni,nj);
Hsd=NaN*ones(ni,nj);
Hmf=NaN*ones(ni,nj);
Hmp=NaN*ones(ni,nj);
Hmd=NaN*ones(ni,nj);
Hlp=NaN*ones(ni,nj);
Hld=NaN*ones(ni,nj);
Hsb =NaN*ones(ni,nj);
Hmb =NaN*ones(ni,nj);
Hsf(ID)=sf_mean;
Hsp(ID)=sp_mean;
Hsd(ID)=sd_mean;
Hmf(ID)=mf_mean;
Hmp(ID)=mp_mean;
Hmd(ID)=md_mean;
Hlp(ID)=lp_mean;
Hld(ID)=ld_mean;
Hsb(ID)=sb_mean;
Hmb(ID)=mb_mean;

clear sf_mean sp_mean sd_mean mf_mean mp_mean md_mean lp_mean ld_mean

%%
CF = Csf+Cmf;
CP = Csp+Cmp+Clp;
CD = Csd+Cmd+Cld;
CS = Csp+Csf+Csd;
CM = Cmp+Cmf+Cmd;
CL = Clp+Cld;
CB = Cb;
CAll = CF + CP + CD;

HF = Hsf+Hmf;
HP = Hsp+Hmp+Hlp;
HD = Hsd+Hmd+Hld;
HS = Hsp+Hsf+Hsd;
HM = Hmp+Hmf+Hmd;
HL = Hlp+Hld;
HB = Hsb+Hmb;
HAll = HF + HP + HD;

%% plot info
geolat = lat;
geolon = lon;
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

%%
pdiffF = (HF-CF) ./ CF;
pdiffP = (HP-CP) ./ CP;
pdiffD = (HD-CD) ./ CD;
pdiffB = (HB-CB) ./ CB;
pdiffAll = (HAll-CAll) ./ CAll;

%% Maps
figure(1)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat,geolon,real(log10(HF)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Clim2B F');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat,geolon,real(log10(CF)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Clim1B F');
print('-dpng',[pp 'Climatol_1B_2B_' harv '_global_F.png'])

%% P
figure(2)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat,geolon,real(log10(HP)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Clim2B P');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat,geolon,real(log10(CP)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Clim1B P');
print('-dpng',[pp 'Climatol_1B_2B_' harv '_global_P.png'])

% D
figure(3)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat,geolon,real(log10(HD)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Clim2B D');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat,geolon,real(log10(CD)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Clim1B D');
print('-dpng',[pp 'Climatol_1B_2B_' harv '_global_D.png'])

%4
figure(4)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat,geolon,real(log10(HAll)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 2]);
colorbar
set(gcf,'renderer','painters')
title('Clim2B All');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat,geolon,real(log10(CAll)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 2]);
colorbar
set(gcf,'renderer','painters')
title('Clim1B All');
print('-dpng',[pp 'Climatol_1B_2B_' harv '_global_All.png'])

%% B
figure(5)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat,geolon,real(log10(HB)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Clim2B Benthos');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat,geolon,real(log10(CB)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Clim1B Benthos');
print('-dpng',[pp 'Climatol_1B_2B_' harv '_global_Bent.png'])

%% 2B
figure(6)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat,geolon,real(log10(Hsb)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Clim2B S Benthos');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat,geolon,real(log10(Hmb)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Clim2B M Benthos');
print('-dpng',[pp 'Climatol_2B_' harv '_global_Bent.png'])

%% pdiffs
%F
figure(7)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat,geolon,pdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Clim2B - Clim1B F');
print('-dpng',[pp 'Climatol_1B_2B_' harv '_global_pdiffF.png'])

%P
figure(8)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat,geolon,pdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Clim2B - Clim1B P');
print('-dpng',[pp 'Climatol_1B_2B_' harv '_global_pdiffP.png'])

figure(9)
%D
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat,geolon,pdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Clim2B - Clim1B D');
print('-dpng',[pp 'Climatol_1B_2B_' harv '_global_pdiffD.png'])

%diif all
figure(10)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat,geolon,pdiffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Clim2B - Clim1B All');
print('-dpng',[pp 'Climatol_1B_2B_' harv '_global_pdiffAll.png'])

%% B
figure(11)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat,geolon,pdiffB)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Clim2B - Clim1B B');
print('-dpng',[pp 'Climatol_1B_2B_' harv '_global_pdiffB.png'])

%% Calc differences in total biomass


