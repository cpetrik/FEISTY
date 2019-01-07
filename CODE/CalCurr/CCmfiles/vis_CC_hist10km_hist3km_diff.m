% Visualize difference between
% 10km and 3km Cal Current models

clear all
close all

%% 10km grid
hpath = '/Volumes/GFDL/NEMURO/10km/';
load([hpath 'gridspec_10km.mat'],'LON','LAT');
load([hpath 'Data_grid_10km_hist.mat']);

geolon_t = LON;
geolat_t = LAT;
GRD10 = GRD;

% plotminlat=30; %Set these bounds for your data
% plotmaxlat=48;
% plotminlon=-134;
% plotmaxlon=-115;

clear LON LAT GRD

%% 3km grid
cpath = '/Volumes/GFDL/NEMURO/3km/';
load([cpath 'gridspec_3km.mat'],'LON','LAT');
load([cpath 'Data_grid_3km_hist.mat']);

GRD3 = GRD;
clear GRD

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
pp = ['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/CalCurr/'];
ppath = [pp cfile '/CC/'];

harv = 'All_fish01';
tharv = 'Harvest all fish 0.1 yr^-^1';

load([fpath 'Means_Historic10km_' harv '_' cfile '.mat'],...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean');
sf_mean10 = sf_mean; 
sp_mean10 = sp_mean;
sd_mean10 = sd_mean;
mf_mean10 = mf_mean;
mp_mean10 = mp_mean;
md_mean10 = md_mean;
lp_mean10 = lp_mean;
ld_mean10 = ld_mean;
b_mean10  = b_mean;

load([fpath 'Means_Historic3km_' harv '_' cfile '.mat'],...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean');
sf_mean3 = sf_mean; 
sp_mean3 = sp_mean;
sd_mean3 = sd_mean;
mf_mean3 = mf_mean;
mp_mean3 = mp_mean;
md_mean3 = md_mean;
lp_mean3 = lp_mean;
ld_mean3 = ld_mean;
b_mean3  = b_mean;


%% Put biomass on grid
%10km
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
Hsf(GRD10.ID)=sf_mean10;
Hsp(GRD10.ID)=sp_mean10;
Hsd(GRD10.ID)=sd_mean10;
Hmf(GRD10.ID)=mf_mean10;
Hmp(GRD10.ID)=mp_mean10;
Hmd(GRD10.ID)=md_mean10;
Hlp(GRD10.ID)=lp_mean10;
Hld(GRD10.ID)=ld_mean10;
Hb(GRD10.ID) =b_mean10;
%3km
[ni,nj]=size(LON);
Csf=NaN*ones(ni,nj);
Csp=NaN*ones(ni,nj);
Csd=NaN*ones(ni,nj);
Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);
Cb =NaN*ones(ni,nj);
Csf(GRD3.ID)=sf_mean3;
Csp(GRD3.ID)=sp_mean3;
Csd(GRD3.ID)=sd_mean3;
Cmf(GRD3.ID)=mf_mean3;
Cmp(GRD3.ID)=mp_mean3;
Cmd(GRD3.ID)=md_mean3;
Clp(GRD3.ID)=lp_mean3;
Cld(GRD3.ID)=ld_mean3;
Cb(GRD3.ID) =b_mean3;

CF = Csf+Cmf;
CP = Csp+Cmp+Clp;
CD = Csd+Cmd+Cld;
CS = Csp+Csf+Csd;
CM = Cmp+Cmf+Cmd;
CL = Clp+Cld;

HF = Hsf+Hmf;
HP = Hsp+Hmp+Hlp;
HD = Hsd+Hmd+Hld;
HS = Hsp+Hsf+Hsd;
HM = Hmp+Hmf+Hmd;
HL = Hlp+Hld;

%% Interpolate to same grid
% geolat = [30 48];
% geolon = [-134 -115]; %1/10 degrees
% LAT = [32 44];
% LON = [-129 -116];    %1/30 degrees

lats = 32:0.05:44;      %1/20 degrees
lons = -129:0.05:-116; 
[glon,glat] = meshgrid(lons,lats);

hF = griddata(geolat_t,geolon_t,HF,glat,glon);
hP = griddata(geolat_t,geolon_t,HP,glat,glon);
hD = griddata(geolat_t,geolon_t,HD,glat,glon);
hS = griddata(geolat_t,geolon_t,HS,glat,glon);
hM = griddata(geolat_t,geolon_t,HM,glat,glon);
hL = griddata(geolat_t,geolon_t,HL,glat,glon);

cF = griddata(LAT,LON,CF,glat,glon);
cP = griddata(LAT,LON,CP,glat,glon);
cD = griddata(LAT,LON,CD,glat,glon);
cS = griddata(LAT,LON,CS,glat,glon);
cM = griddata(LAT,LON,CM,glat,glon);
cL = griddata(LAT,LON,CL,glat,glon);

%% plot info
plotminlat=32; %Set these bounds for your data
plotmaxlat=44;
plotminlon=-129;
plotmaxlon=-116;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%%
cAll = cF+cP+cD;
cFracPD = cP ./ (cP+cD);
cFracPF = cP ./ (cP+cF);
cFracLM = cL ./ (cL+cM);

hAll = hF+hP+hD;
hFracPD = hP ./ (hP+hD);
hFracPF = hP ./ (hP+hF);
hFracLM = hL ./ (hL+hM);

pdiffF = (cF-hF) ./ cF;
pdiffP = (cP-hP) ./ cP;
pdiffD = (cD-hD) ./ cD;
pdiffAll = (cAll-hAll) ./ cAll;

diffF = (cF-hF);
diffP = (cP-hP);
diffD = (cD-hD);
diffAll = (cAll-hAll);

ldiffF = (log10(cF)-log10(hF));
ldiffP = (log10(cP)-log10(hP));
ldiffD = (log10(cD)-log10(hD));
ldiffAll = (log10(cAll)-log10(hAll));

%% Maps
figure(1)
subplot(2,1,1)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(hF)))
colormap('jet')
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('10km F');

subplot(2,1,2)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(cF)))
colormap('jet')
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('3km F');
print('-dpng',[ppath 'Hist_10km_3km_' harv '_global_F.png'])

%P
figure(2)
subplot(2,1,1)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(hP)))
colormap('jet')
caxis([0 2]);
colorbar
set(gcf,'renderer','painters')
title('10km P');

subplot(2,1,2)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(cP)))
colormap('jet')
caxis([0 2]);
colorbar
set(gcf,'renderer','painters')
title('3km P');
print('-dpng',[ppath 'Hist_10km_3km_' harv '_global_P.png'])

% D
figure(3)
subplot(2,1,1)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(hD)))
colormap('jet')
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('10km D');

subplot(2,1,2)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(cD)))
colormap('jet')
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('3km D');
print('-dpng',[ppath 'Hist_10km_3km_' harv '_global_D.png'])

%4
figure(4)
subplot(2,1,1)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(hAll)))
colormap('jet')
caxis([0 2]);
colorbar
set(gcf,'renderer','painters')
title('10km All');

subplot(2,1,2)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(cAll)))
colormap('jet')
caxis([0 2]);
colorbar
set(gcf,'renderer','painters')
title('3km All');
print('-dpng',[ppath 'Hist_10km_3km_' harv '_global_All.png'])

%% 5
figure(5)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,pdiffF)
cmocean('balance')
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Percent difference 3km - 10km F');
print('-dpng',[ppath 'Hist_10km_3km_' harv '_global_pdiffF.png'])

%6
figure(6)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,pdiffP)
cmocean('balance')
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Percent difference 3km - 10km P');
print('-dpng',[ppath 'Hist_10km_3km_' harv '_global_pdiffP.png'])

figure(7)
%7
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,pdiffD)
cmocean('balance')
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Percent difference 3km - 10km D');
print('-dpng',[ppath 'Hist_10km_3km_' harv '_global_pdiffD.png'])

%8
figure(8)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,pdiffAll)
cmocean('balance')
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Percent difference 3km - 10km All');
print('-dpng',[ppath 'Hist_10km_3km_' harv '_global_pdiffAll.png'])

%% Subplot
figure(9)
subplot('Position',[0 0.5 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,pdiffF)
cmocean('balance')
caxis([-1 1]);
set(gcf,'renderer','painters')
title('Percent difference 3km - 10km F');

subplot('Position',[0.5 0.5 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,pdiffP)
cmocean('balance')
caxis([-1 1]);
set(gcf,'renderer','painters')
title('Percent difference 3km - 10km P');

subplot('Position',[0 0 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,pdiffD)
cmocean('balance')
caxis([-1 1]);
colorbar('Position',[0.45 0.25 0.03 0.5],'orientation','vertical')
set(gcf,'renderer','painters')
title('Percent difference 3km - 10km D');

subplot('Position',[0.5 0 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,pdiffAll)
cmocean('balance')
caxis([-1 1]);
set(gcf,'renderer','painters')
title('Percent difference 3km - 10km All');
print('-dpng',[ppath 'Hist_10km_3km_' harv '_global_pdiff_subplot.png'])


%% Subplot
figure(10)
subplot('Position',[0 0.5 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,diffF)
cmocean('balance')
caxis([-30 30]);
set(gcf,'renderer','painters')
title('Difference 3km - 10km F');

subplot('Position',[0.5 0.5 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,diffP)
cmocean('balance')
caxis([-30 30]);
set(gcf,'renderer','painters')
title('Difference 3km - 10km P');

subplot('Position',[0 0 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,diffD)
cmocean('balance')
caxis([-30 30]);
colorbar('Position',[0.45 0.25 0.03 0.5],'orientation','vertical')
set(gcf,'renderer','painters')
title('Difference 3km - 10km D');

subplot('Position',[0.5 0 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,diffAll)
cmocean('balance')
caxis([-30 30]);
set(gcf,'renderer','painters')
title('Difference 3km - 10km All');
print('-dpng',[ppath 'Hist_10km_3km_' harv '_global_diff_subplot.png'])


figure(11)
subplot('Position',[0 0.5 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,(ldiffF))
cmocean('balance')
caxis([-2 2]);
set(gcf,'renderer','painters')
title('log_1_0 Difference 3km - 10km F');

subplot('Position',[0.5 0.5 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,(ldiffP))
cmocean('balance')
caxis([-2 2]);
set(gcf,'renderer','painters')
title('log_1_0 Difference 3km - 10km P');

subplot('Position',[0 0 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,real(ldiffD))
cmocean('balance')
caxis([-2 2]);
colorbar('Position',[0.45 0.25 0.03 0.5],'orientation','vertical')
set(gcf,'renderer','painters')
title('log_1_0 Difference 3km - 10km D');

subplot('Position',[0.5 0 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,(ldiffAll))
cmocean('balance')
caxis([-2 2]);
set(gcf,'renderer','painters')
title('log_1_0 Difference 3km - 10km All');
print('-dpng',[ppath 'Hist_10km_3km_' harv '_global_logdiff_subplot.png'])


