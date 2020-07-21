% Visualize output of FEISTY
% Climatology at all locations
% Fraction of zoop hp loss consumed 
% And number of times overconsumed

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/Climatology/'];
%fpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

ppath = [pp cfile '/'];

load([fpath 'Means_Climatol_' harv '_' cfile '.mat'],...
    'mz_tmfrac','mz_mfrac5','mz_ttf','mz_mtf5','mz_mtf',...
    'lz_tmfrac','lz_mfrac5','lz_ttf','lz_mtf5','lz_mtf');

load('/Volumes/FEISTY/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_gridspec.mat',...
    'lon','lat','ID');

[ni,nj]=size(lon);
nx = length(ID);

%% Pick which time period mean
% 1951-2000
%Mean fraction
Cmz_smfrac = mz_mfrac5;
Clz_smfrac = lz_mfrac5;

%Total times it happens over time
Cmz_ttover = mz_ttf/nx;
Clz_ttover = lz_ttf/nx;

%Total times it happens in last year (12 mos) in space
Cmz_stover5 = mz_mtf5 ./ (12);
Clz_stover5 = lz_mtf5 ./ (12);

%%
%happens whole year
test=floor(Clz_stover5);
%histogram(test)
sum(test)/length(test) % = 0.3154

test2=floor(Cmz_stover5);
%histogram(test2)
sum(test2)/length(test2) % = 4.4322e-05 (2 grid cells total)

%happens >=50% of year
test3=round(Clz_stover5);
%histogram(test3)
sum(test3)/length(test3) % = 0.6006

test4=round(Cmz_stover5);
%histogram(test4)
sum(test4)/length(test4) % = 0.0703

%% Plot in time
y=(1/12):(1/12):150;
figure(1)
plot(y, Cmz_ttover,'r','LineWidth',2); hold on;
plot(y, Clz_ttover,'b','LineWidth',2); hold on;
%xlim([1951 2005])
legend('M','L')
xlabel('Years')
ylabel('Fraction of grid cells over-consumed')
print('-dpng',[ppath 'Clim_' harv '_timeseries_zoop_overcon.png'])

%% Plots in space

CFmz=NaN*ones(ni,nj);
COmz=NaN*ones(ni,nj);
CFlz=NaN*ones(ni,nj);
COlz=NaN*ones(ni,nj);
COmz5=NaN*ones(ni,nj);
COlz5=NaN*ones(ni,nj);

CFmz(ID)=Cmz_smfrac;
CFlz(ID)=Clz_smfrac;

COmz(ID)=Cmz_stover5;
COlz(ID)=Clz_stover5;

bpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/ESM2M_Hist_Fore/'];
save([bpath 'Hist_Fore_' harv '_ts_map_zoop_overcon.mat'],...
    'Cmz_ttover','Clz_ttover','CFmz','COmz','CFlz','COlz','-append');

%% plot info
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

cmBP=cbrewer('seq','BuPu',10,'PCHIP');
cmBP2 = cmBP;
cmBP2(11,:) = [0 0 0];

%% 2 plots of both maps
figure(3)
%1 - m frac
subplot('Position',[0.01 0.5 0.9 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,CFmz)
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
surfm(geolat_t,geolon_t,CFlz)
colormap(cmBP2)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1.1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'Mean fraction LZ hploss consumed','HorizontalAlignment','center')
%text(-2.75,1.75,'C')
print('-dpng',[ppath 'Clim_' harv '_global_zoop_fraccon.png'])

figure(4)
%2 - m over
subplot('Position',[0.01 0.5 0.9 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,COmz)
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
surfm(geolat_t,geolon_t,COlz)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'Mean times LZ overconsumed','HorizontalAlignment','center')
%text(-2.75,1.75,'D')
print('-dpng',[ppath 'Clim_' harv '_global_zoop_overcon.png'])








