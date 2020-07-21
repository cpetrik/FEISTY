% Visualize output of FEISTY
% Projection time period (2006-2100) at all locations
% Fraction of zoop hp loss consumed 
% And number of times overconsumed

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/Forecast_RCP85_ESM2M/'];
%fpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

ppath = [pp cfile '/'];

load([fpath 'Means_fore_' harv '_' cfile '.mat'],...
    'mz_tmfrac','mz_mfrac50','mz_mfrac','mz_ttf',...
    'mz_mtf50','mz_mtf',...
    'lz_tmfrac','lz_mfrac50','lz_mfrac','lz_ttf',...
    'lz_mtf50','lz_mtf');

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

[ni,nj]=size(geolon_t);
ID = grid(:,1);
nx = length(ID);

y = 2005+(1/12):(1/12):2100;

%% Pick which time period mean
% 2051-2100
%Mean fraction
mz_smfrac = mz_mfrac50;
lz_smfrac = lz_mfrac50;

%Total times it happens over time
Fmz_ttover = mz_ttf/nx;
Flz_ttover = lz_ttf/nx;

%Total times it happens in 50 yrs * 12 mos in space
mz_stover = mz_mtf50 ./ (50*12);
lz_stover = lz_mtf50 ./ (50*12);

%% Plot in time
figure(1)
plot(y, Fmz_ttover,'r','LineWidth',2); hold on;
plot(y, Flz_ttover,'b','LineWidth',2); hold on;
xlim([2006 2100])
legend('M','L')
xlabel('Year')
ylabel('Fraction of grid cells over-consumed')
print('-dpng',[ppath 'Fore_' harv '_timeseries_zoop_overcon.png'])

%% Plots in space

FFmz=NaN*ones(ni,nj);
FOmz=NaN*ones(ni,nj);
FFlz=NaN*ones(ni,nj);
FOlz=NaN*ones(ni,nj);

FFmz(ID)=mz_smfrac;
FFlz(ID)=lz_smfrac;

FOmz(ID)=mz_stover;
FOlz(ID)=lz_stover;

bpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/ESM2M_Hist_Fore/'];
save([bpath 'Hist_Fore_' harv '_ts_map_zoop_overcon.mat'],...
    'Fmz_ttover','Flz_ttover','FFmz','FOmz','FFlz','FOlz','-append');

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

cmBP=cbrewer('seq','BuPu',10,'PCHIP');
cmBP2 = cmBP;
cmBP2(11,:) = [0 0 0];

%% 4 plot of all maps
figure(2)
%1 - m frac
subplot('Position',[0.01 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FFmz)
colormap(cmBP2)
%cmocean('oxy')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1.1]);
colorbar
%colorbar('Position',[0.385 0.695 0.025 0.275],'orientation','vertical','AxisLocation','out')
%set(gcf,'renderer','painters')
text(0,1.75,'Mean fraction MZ hploss consumed','HorizontalAlignment','center')
text(-2.75,1.75,'A')

%2 - m over
subplot('Position',[0.475 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FOmz)
colormap(cmBP2)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
%colorbar('Position',[0.825 0.695 0.025 0.275],'orientation','vertical','AxisLocation','out')
%set(gcf,'renderer','painters')
text(0,1.75,'Mean times MZ overconsumed','HorizontalAlignment','center')
text(-2.75,1.75,'B')

%3 - l frac
subplot('Position',[0.01 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FFlz)
colormap(cmBP2)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1.1]);
colorbar
%colorbar('Position',[0.385 0.075 0.025 0.275],'orientation','vertical','AxisLocation','out')
%set(gcf,'renderer','painters')
text(0,1.75,'Mean fraction LZ hploss consumed','HorizontalAlignment','center')
text(-2.75,1.75,'C')

%4 - l over
subplot('Position',[0.475 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FOlz)
colormap(cmBP2)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
%colorbar('Position',[0.825 0.075 0.025 0.275],'orientation','vertical','AxisLocation','out')
%set(gcf,'renderer','painters')
text(0,1.75,'Mean times MZ overconsumed','HorizontalAlignment','center')
text(-2.75,1.75,'D')
%print('-dpng',[ppath 'Fore_' harv '_global_zoop_overcon.png'])

%% 2 plots of both maps
figure(3)
%1 - m frac
subplot('Position',[0.01 0.5 0.9 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FFmz)
colormap(cmBP2)
%cmocean('oxy')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1.1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'Mean fraction MZ hploss consumed','HorizontalAlignment','center')
%text(-2.75,1.75,'A')

%3 - l frac
subplot('Position',[0.01 0.01 0.9 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FFlz)
colormap(cmBP2)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1.1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'Mean fraction LZ hploss consumed','HorizontalAlignment','center')
%text(-2.75,1.75,'C')
print('-dpng',[ppath 'Fore_' harv '_global_zoop_fraccon.png'])

%%
figure(4)
%2 - m over
subplot('Position',[0.01 0.5 0.9 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FOmz)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'Mean times MZ overconsumed','HorizontalAlignment','center')
%text(-2.75,1.75,'B')

%4 - l over
subplot('Position',[0.01 0.01 0.9 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FOlz)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'Mean times LZ overconsumed','HorizontalAlignment','center')
%text(-2.75,1.75,'D')
print('-dpng',[ppath 'Fore_' harv '_global_zoop_overcon.png'])







