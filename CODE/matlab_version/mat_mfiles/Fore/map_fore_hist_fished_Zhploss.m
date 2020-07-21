% Visualize difference between
% ESM2M Hindcast of 1951-2000 and Forecast of 2051-2100
% Zoop hploss overconsumption

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/ESM2M_Hist_Fore/'];
%fpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';

% saved hist and fore together
load([fpath 'Hist_Fore_' harv '_ts_map_zoop_overcon.mat']);

hy = 1860+(1/12):(1/12):2005;
fy = 2005+(1/12):(1/12):2100;
y = 1860+(1/12):(1/12):2100;

%% 
%full ts
mz_ttover = [Hmz_ttover Fmz_ttover];
lz_ttover = [Hlz_ttover Flz_ttover];
mmMZ = movmean(mz_ttover,12);
mmLZ = movmean(lz_ttover,12);

%diffs
pdiffFmz = (FFmz-HFmz) ./ HFmz;
pdiffFlz = (FFlz-HFlz) ./ HFlz;
diffOmz = (FOmz-HOmz);
diffOlz = (FOlz-HOlz);

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

%% Time series
figure(1)
% plot(y, mz_ttover,'r','LineWidth',2); hold on;
% plot(y, lz_ttover,'b','LineWidth',2); hold on;
plot(y, mmMZ,'r','LineWidth',2); hold on;
plot(y, mmLZ,'b','LineWidth',2); hold on;
xlim([1951 2100])
legend('M','L')
xlabel('Year')
ylabel('Fraction of grid cells over-consumed')
print('-dpng',[pp 'Hist_Fore_' harv '_timeseries_zoop_overcon.png'])

figure(2)
subplot(2,1,1)
plot(y, mmMZ,'r','LineWidth',2); hold on;
xlim([1951 2100])
ylim([0.2 0.3])
ylabel('Fraction of grid cells over-consumed')

subplot(2,1,2)
plot(y, mmLZ,'b','LineWidth',2); hold on;
xlim([1951 2100])
ylim([0.61 0.71])
xlabel('Year')
ylabel('Fraction of grid cells over-consumed')
print('-dpng',[pp 'Hist_Fore_' harv '_timeseries_zoop_overcon_sub.png'])

%% Maps
% 
figure(3)
%A - F
subplot('Position',[0.01 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,pdiffFmz)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.75,'% change frac mz','HorizontalAlignment','center')

%B - B
subplot('Position',[0.5 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,pdiffFlz)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.75,'% change frac lz','HorizontalAlignment','center')

%C - P
subplot('Position',[0.01 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffOmz)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
text(0,1.75,'Diff over mz','HorizontalAlignment','center')

%D - D
subplot('Position',[0.5 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffOlz)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
text(0,1.75,'Diff over lz','HorizontalAlignment','center')
print('-dpng',[pp 'Hist_Fore_',harv,'_global_diff_zoop_overcon.png'])

