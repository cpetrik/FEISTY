function map_hist_1meso_overcon_diff_loop(sF,sP,sD,sB,sM,sL,ppath)

% Visualize difference between
% ESM2M Hindcast w/1 meso, 
% with and without overconsump from HPloss

close all

%% Hindcast grid
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% 1 meso structurally FEISTY
cfile1 = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
fpath1 = ['/Volumes/FEISTY/NC/Matlab_new_size/' cfile1 '/'];
harv = 'All_fish03';

load([fpath1 'Means_Historic_1meso_',harv,'_' cfile1 '.mat'],...
    'sf_mean50','sp_mean50','sd_mean50',...
    'mf_mean50','mp_mean50','md_mean50',...
    'lp_mean50','ld_mean50','b_mean50');

[ni,nj]=size(geolon_t);
CF=NaN*ones(ni,nj);
CP=NaN*ones(ni,nj);
CD=NaN*ones(ni,nj);
CB=NaN*ones(ni,nj);
CL=NaN*ones(ni,nj);
CM=NaN*ones(ni,nj);
CF(ID) = sf_mean50 + mf_mean50;
CP(ID) = sp_mean50 + mp_mean50 + lp_mean50;
CD(ID) = sd_mean50 + md_mean50 + ld_mean50;
CB(ID) = b_mean50;
CM(ID) = mf_mean50 + mp_mean50 + md_mean50;
CL(ID) = lp_mean50 + ld_mean50;

CAll = CF+CP+CD;
cFracPD = CP ./ (CP+CD);
cFracPF = CP ./ (CP+CF);
cFracLM = CL ./ (CL+CM);

%% 1 meso without overcon
ZF=NaN*ones(ni,nj);
ZP=NaN*ones(ni,nj);
ZD=NaN*ones(ni,nj);
ZM=NaN*ones(ni,nj);
ZL=NaN*ones(ni,nj);
ZB=NaN*ones(ni,nj);

ZF(grid(:,1))=sF;
ZP(grid(:,1))=sP;
ZD(grid(:,1))=sD;
ZM(grid(:,1))=sM;
ZL(grid(:,1))=sL;
ZB(grid(:,1))=sB;

ZAll = ZF+ZP+ZD;
zFracPD = ZP ./ (ZP+ZD);
zFracPF = ZP ./ (ZP+ZF);
zFracLM = ZL ./ (ZL+ZM);

%%
pdiffF = (CF-ZF) ./ CF;
pdiffP = (CP-ZP) ./ CP;
pdiffD = (CD-ZD) ./ CD;
pdiffB = (CB-ZB) ./ CB;
pdiffAll = (CAll-ZAll) ./ CAll;

dFracPD = cFracPD - zFracPD;
dFracPF = cFracPF - zFracPF;
dFracLM = cFracLM - zFracLM;

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%% Maps
% diffs -------------------------------------------------------------
%F
figure(1)
subplot('Position',[0.01 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('F HPloss - no HPloss');

%P
subplot('Position',[0.475 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('P HPloss - no HPloss');

%D
subplot('Position',[0.01 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('D HPloss - no HPloss');

%diff all
subplot('Position',[0.475 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('All HPloss - no HPloss');
print('-dpng',[ppath 'Meso1_HPloss_noHP_' harv '_global_pdiff_fish.png'])

%% B
figure(2)
%P:D
subplot('Position',[0.01 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,dFracPD)
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
surfm(geolat_t,geolon_t,dFracPF)
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
surfm(geolat_t,geolon_t,dFracLM)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('L/(L+M)');

%B
subplot('Position',[0.475 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffB)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('B HPloss - no HPloss');
print('-dpng',[ppath 'Meso1_HPloss_noHP_' harv '_global_pdiff_fracs.png'])



