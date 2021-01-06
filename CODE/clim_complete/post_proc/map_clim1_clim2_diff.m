% Visualize difference between
% ESM2.6 Climatology of 5 yrs (w/1990 ICs) and
% Climatology with NoNuUpdate

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';

%% Climatology grid
Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']); %depth,ID,lat,lmask,lon
load([cpath 'esm26_area_1deg.mat']); %area
AREA_OCN = max(area,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/Climatology/'];
dpath=['/Volumes/FEISTY/NC/Clim_comp_tests/' cfile '/NoNuUpdate_'];
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/clim_complete/post_proc/pp_figs/NoNuUpdate_';

harv = 'All_fish03';

%% Orig climatol
load([fpath 'Means_Climatol_' harv '_' cfile '.mat'],...
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

CF = Csf+Cmf;
CP = Csp+Cmp+Clp;
CD = Csd+Cmd+Cld;
CS = Csp+Csf+Csd;
CM = Cmp+Cmf+Cmd;
CL = Clp+Cld;

%% NoNuUpdate
load([dpath 'Means_Climatol_',harv,'_' cfile '.mat'],...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean');

% Put biomass on grid
[hi,hj]=size(lon);
Hsf=NaN*ones(hi,hj);
Hsp=NaN*ones(hi,hj);
Hsd=NaN*ones(hi,hj);
Hmf=NaN*ones(hi,hj);
Hmp=NaN*ones(hi,hj);
Hmd=NaN*ones(hi,hj);
Hlp=NaN*ones(hi,hj);
Hld=NaN*ones(hi,hj);
Hb =NaN*ones(hi,hj);
Hsf(ID)=sf_mean;
Hsp(ID)=sp_mean;
Hsd(ID)=sd_mean;
Hmf(ID)=mf_mean;
Hmp(ID)=mp_mean;
Hmd(ID)=md_mean;
Hlp(ID)=lp_mean;
Hld(ID)=ld_mean;
Hb(ID) =b_mean;

HF = Hsf+Hmf;
HP = Hsp+Hmp+Hlp;
HD = Hsd+Hmd+Hld;
HS = Hsp+Hsf+Hsd;
HM = Hmp+Hmf+Hmd;
HL = Hlp+Hld;

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

glat = lat;
glon = lon;

%%
CAll = CF+CP+CD;
CFracPD = CP ./ (CP+CD);
CFracPF = CP ./ (CP+CF);
CFracLM = CL ./ (CL+CM);

HAll = HF+HP+HD;
HFracPD = HP ./ (HP+HD);
HFracPF = HP ./ (HP+HF);
HFracLM = HL ./ (HL+HM);

diffF = (CF-HF) ./ CF;
diffP = (CP-HP) ./ CP;
diffD = (CD-HD) ./ CD;
diffB = (Cb-Hb) ./ Cb;
diffAll = (CAll-HAll) ./ CAll;

dFracPD = CFracPD - HFracPD;
dFracPF = CFracPF - HFracPF;
dFracLM = CFracLM - HFracLM;

%% Maps
% F
figure(1)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(HF)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('NoNuUpdate F');

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(CF)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Climatology F');

% diffs
subplot('Position',[0.5 0.25 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,100*diffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('Climatology - NoNuUpdate F');
print('-dpng',[pp 'Climatol_' harv '_global_F.png'])

%% P
figure(2)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(HP)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('NoNuUpdate P');

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(CP)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Climatology P');

% diffs
subplot('Position',[0.5 0.25 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,100*diffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('Climatology - NoNuUpdate P');
print('-dpng',[pp 'Climatol_' harv '_global_P.png'])

%% D
figure(3)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(HD)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('NoNuUpdate D');

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(CD)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Climatology D');

% diffs
subplot('Position',[0.5 0.25 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,100*diffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('Climatology - NoNuUpdate D');
print('-dpng',[pp 'Climatol_' harv '_global_D.png'])

%% 4
figure(4)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(HAll)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('NoNuUpdate All');

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(CAll)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Climatology All');

% diffs
subplot('Position',[0.5 0.25 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,100*diffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('Climatology - NoNuUpdate All');
print('-dpng',[pp 'Climatol_' harv '_global_All.png'])

%% B
figure(5)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(Hb)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('NoNuUpdate Benthos');

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(Cb)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Climatology Benthos');

% diffs
subplot('Position',[0.5 0.25 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,diffB)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Climatology - NoNuUpdate B');
print('-dpng',[pp 'Climatol_' harv '_global_B.png'])

%% Diffs & Fracs -------------------------------------------------------------
%F
figure(6)
subplot('Position',[0.01 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,100*diffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('F Climatology - NoNuUpdate ');

%P
subplot('Position',[0.475 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,100*diffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('P Climatology - NoNuUpdate ');

%D
subplot('Position',[0.01 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,100*diffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('D Climatology - NoNuUpdate ');

%diff all
subplot('Position',[0.475 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,100*diffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('All Climatology - NoNuUpdate ');
print('-dpng',[pp 'Climatol_' harv '_global_pdiff_fish.png'])

%% Fracs
figure(7)
%P:D
subplot('Position',[0.01 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,dFracPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('P/(P+D)');

%P:F
subplot('Position',[0.475 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,dFracPF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('P/(P+F)');

%L:M
subplot('Position',[0.01 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,dFracLM)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('L/(L+M)');
print('-dpng',[pp 'Climatol_' harv '_global_diff_fracs.png'])








