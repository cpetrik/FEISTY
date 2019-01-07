% Plot mean biomasses of zoop, det, temp

clear all
close all

cpath = '/Volumes/GFDL/NEMURO/3km/';
load([cpath 'nemuro_3km_means.mat'])

% Figures
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = [pp cfile '/CC/'];

%% Convert units
%from gWW m-2 to mg C m-2
%from gWW m-2 d-1 to mg C m-2 d-1
% 106/16 mol C in 1 mol N
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
gc_m2 = (1.0/9.0) * 1e3;

% det_mean_3km = det_mean_3km *gc_m2;
% mz_mean_3km = mz_mean_3km *gc_m2;
% lz_mean_3km = lz_mean_3km *gc_m2;
% mzloss_mean_3km = mz_mean_3km * 0.06;
% lzloss_mean_3km = lz_mean_3km * 0.08;

%% Map data
load([cpath 'gridspec_3km.mat'],'LON','LAT');
load([cpath 'Data_grid_3km_hist.mat']);
[ni,nj]=size(LON);
plotminlat=32; 
plotmaxlat=44;
plotminlon=-129;
plotmaxlon=-116;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%% Plots in space
MZ=NaN*ones(ni,nj);
LZ=NaN*ones(ni,nj);
MZloss=NaN*ones(ni,nj);
LZloss=NaN*ones(ni,nj);
Tp=NaN*ones(ni,nj);
Tb=NaN*ones(ni,nj);
Det=NaN*ones(ni,nj);

MZ(GRD.ID)=mz_mean_3km;
LZ(GRD.ID)=lz_mean_3km;
MZloss(GRD.ID)=mzloss_mean_3km;
LZloss(GRD.ID)=lzloss_mean_3km;
Tp(GRD.ID)=ptemp_mean_3km;
Tb(GRD.ID)=btemp_mean_3km;
Det(GRD.ID)=det_mean_3km;

%%
figure(1)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(Tp))
colormap('jet')
caxis([1 17]);
colorbar('h');
set(gcf,'renderer','painters')
title('mean Top 100 m Temp (^oC)')
print('-dpng',[ppath 't100_nemuro3km.png'])

figure(2)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(Tb))
colormap('jet')
caxis([1 17]);
colorbar('h');
set(gcf,'renderer','painters')
title('mean Bottom Temp (^oC)')
print('-dpng',[ppath 'btemp_nemuro3km.png'])

figure(3)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(Det))
colormap('jet')
caxis([-1 1]);
colorbar('h');
set(gcf,'renderer','painters')
title('log_1_0 mean Detritus flux (g WW m^-^2 d^-^1)')
print('-dpng',[ppath 'LZ_loss_nemuro3km.png'])
%%
figure(4)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(MZ))
colormap('jet')
caxis([0 2]);
colorbar('h');
set(gcf,'renderer','painters')
title('log_1_0 mean MZ biomass (g WW m^-^2)')
print('-dpng',[ppath 'MZ_nemuro3km.png'])

figure(5)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(LZ))
colormap('jet')
caxis([0 2]);
colorbar('h');
set(gcf,'renderer','painters')
title('log_1_0 mean LZ biomass (g WW m^-^2)')
print('-dpng',[ppath 'LZ_nemuro3km.png'])

figure(8)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(MZ+LZ))
colormap('jet')
caxis([0 2]);
colorbar('h');
set(gcf,'renderer','painters')
title('log_1_0 mean Zoop biomass (g WW m^-^2)')
print('-dpng',[ppath 'Z_nemuro3km.png'])

%%
figure(6)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(MZloss))
colormap('jet')
caxis([-1 1]);
colorbar('h');
set(gcf,'renderer','painters')
title('log_1_0 mean MZ biomass loss (g WW m^-^2 d^-^1)')
print('-dpng',[ppath 'MZ_loss_nemuro3km.png'])

figure(7)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(LZloss))
colormap('jet')
caxis([-1 1]);
colorbar('h');
set(gcf,'renderer','painters')
title('log_1_0 mean LZ biomass loss (g WW m^-^2 d^-^1)')
print('-dpng',[ppath 'LZ_loss_nemuro3km.png'])

figure(9)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(MZloss+LZloss))
colormap('jet')
caxis([-1 1]);
colorbar('h');
set(gcf,'renderer','painters')
%title('log_1_0 mean Zoop biomass loss (mg C m^-^2 d^-^1)')
title('log_1_0 mean Zoop biomass loss (g WW m^-^2 d^-^1)')
print('-dpng',[ppath 'Z_loss_nemuro3km.png'])


