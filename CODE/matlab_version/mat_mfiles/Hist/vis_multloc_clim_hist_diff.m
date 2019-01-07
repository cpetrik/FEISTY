% Visualize difference between
% ESM2.6 Climatology of 5 yrs (w/1990 ICs) and
% ESM2M Hindcast of 1990-1994

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';

%% Climatology grid
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']); %depth,ID,lat,lmask,lon
load([cpath 'esm26_area_1deg.mat']); %area
AREA_OCN = max(area,1);

%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

load([fpath 'Means_Climatol_' harv '_' cfile '.mat'],...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean');

load([fpath 'Means_Historic_',harv,'_' cfile '.mat'],...
    'sf_mean5','sp_mean5','sd_mean5',...
    'mf_mean5','mp_mean5','md_mean5',...
    'lp_mean5','ld_mean5','b_mean5');

%% Put biomass on grid
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
Hsf(grid(:,1))=sf_mean5;
Hsp(grid(:,1))=sp_mean5;
Hsd(grid(:,1))=sd_mean5;
Hmf(grid(:,1))=mf_mean5;
Hmp(grid(:,1))=mp_mean5;
Hmd(grid(:,1))=md_mean5;
Hlp(grid(:,1))=lp_mean5;
Hld(grid(:,1))=ld_mean5;
Hb(grid(:,1)) =b_mean5;

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
%lat        [-89.5 89.5]
%geolat_t   [-81.5 89.4879]
%lon        [0.5 359.5]
%geolon_t   [-279.9803 79.9803]

% Need to fix both longitudes
test = lon-360;
id=find(test<-180);
test(id)=test(id)+360;
lon = test;

test2=geolon_t;
id=find(test2<-180);
test2(id)=test2(id)+360;
geolon_t = test2;

geolat = geolat_t';
geolon = geolon_t';

lats = -89.5:89.5;
lons = -179.5:179.5;
[glon,glat] = meshgrid(lons,lats);

clat = lat;
clon = lon;

hF = griddata(geolat,geolon,HF',glat,glon);
hP = griddata(geolat,geolon,HP',glat,glon);
hD = griddata(geolat,geolon,HD',glat,glon);
hS = griddata(geolat,geolon,HS',glat,glon);
hM = griddata(geolat,geolon,HM',glat,glon);
hL = griddata(geolat,geolon,HL',glat,glon);

cF = griddata(clat,clon,CF,glat,glon);
cP = griddata(clat,clon,CP,glat,glon);
cD = griddata(clat,clon,CD,glat,glon);
cS = griddata(clat,clon,CS,glat,glon);
cM = griddata(clat,clon,CM,glat,glon);
cL = griddata(clat,clon,CL,glat,glon);

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

%%
cAll = cF+cP+cD;
cFracPD = cP ./ (cP+cD);
cFracPF = cP ./ (cP+cF);
cFracLM = cL ./ (cL+cM);

hAll = hF+hP+hD;
hFrahPD = hP ./ (hP+hD);
hFrahPF = hP ./ (hP+hF);
hFrahLM = hL ./ (hL+hM);

diffF = (cF-hF) ./ cF;
diffP = (cP-hP) ./ cP;
diffD = (cD-hD) ./ cD;
diffAll = (cAll-hAll) ./ cAll;

%% Maps
figure(1)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(hF)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast F');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(cF)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Climatology F');
print('-dpng',[pp 'Hist_Climatol_' harv '_global_F.png'])

%P
figure(2)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(hP)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast P');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(cP)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Climatology P');
print('-dpng',[pp 'Hist_Climatol_' harv '_global_P.png'])

% D
figure(3)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(hD)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast D');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(cD)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Climatology D');
print('-dpng',[pp 'Hist_Climatol_' harv '_global_D.png'])

%4
figure(4)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(hAll)))
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
surfm(glat,glon,real(log10(cAll)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 2]);
colorbar
set(gcf,'renderer','painters')
title('Climatology All');
print('-dpng',[pp 'Hist_Climatol_' harv '_global_All.png'])

%% 5
figure(5)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,diffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Climatology - Hindcast F');
print('-dpng',[pp 'Hist_Climatol_' harv '_global_diffF.png'])

%6
figure(6)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,diffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Climatology - Hindcast P');
print('-dpng',[pp 'Hist_Climatol_' harv '_global_diffP.png'])

figure(7)
%7
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,diffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Climatology - Hindcast D');
print('-dpng',[pp 'Hist_Climatol_' harv '_global_diffD.png'])

%8
figure(8)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,diffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Climatology - Hindcast All');
print('-dpng',[pp 'Hist_Climatol_' harv '_global_diffAll.png'])

%% Calc differences in total biomass


