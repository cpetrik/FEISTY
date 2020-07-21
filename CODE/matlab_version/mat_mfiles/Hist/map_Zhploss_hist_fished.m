% Visualize output of FEISTY
% Historic time period (1861-2005) at all locations
% Fraction of zoop hp loss consumed 
% And number of times overconsumed

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/Historic_ESM2M/'];
%fpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

ppath = [pp cfile '/'];

load([fpath 'Means_Historic_' harv '_' cfile '.mat'],'y',...
    'mz_tmfrac','mz_mfrac50','mz_mfrac5','mz_mfrac90','mz_mfrac','mz_ttf',...
    'mz_mtf50','mz_mtf5','mz_mtf90','mz_mtf',...
    'lz_tmfrac','lz_mfrac50','lz_mfrac5','lz_mfrac90','lz_mfrac','lz_ttf',...
    'lz_mtf50','lz_mtf5','lz_mtf90','lz_mtf');

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

[ni,nj]=size(geolon_t);
ID = grid(:,1);
nx = length(ID);

%% Pick which time period mean
% 1951-2000
%Mean fraction
Hmz_smfrac = mz_mfrac50;
Hlz_smfrac = lz_mfrac50;

%Total times it happens over time
Hmz_ttover = mz_ttf/nx;
Hlz_ttover = lz_ttf/nx;

%Total times it happens in 50 yrs * 12 mos in space
Hmz_stover = mz_mtf50 ./ (50*12);
Hlz_stover = lz_mtf50 ./ (50*12);

%Total times it happens in 1990 (12 mos) in space
Hmz_stover90 = mz_mtf90 ./ (12);
Hlz_stover90 = lz_mtf90 ./ (12);

%%
%happens whole year
test=floor(Hlz_stover90);
%histogram(test)
sum(test)/length(test) % = 0.3228

test2=floor(Hmz_stover90);
%histogram(test2)
sum(test2)/length(test2) % = 1.4550e-04 (7 grid cells total)

%happens >=50% of year
test3=round(Hlz_stover90);
%histogram(test3)
sum(test3)/length(test3) % = 0.6820

test4=round(Hmz_stover90);
%histogram(test4)
sum(test4)/length(test4) % = 0.0822

%% Plot in time
figure(1)
plot(y, Hmz_ttover,'r','LineWidth',2); hold on;
plot(y, Hlz_ttover,'b','LineWidth',2); hold on;
xlim([1951 2005])
legend('M','L')
xlabel('Year')
ylabel('Fraction of grid cells over-consumed')
print('-dpng',[ppath 'Hist_' harv '_timeseries_zoop_overcon.png'])

%% Plots in space

HFmz=NaN*ones(ni,nj);
HOmz=NaN*ones(ni,nj);
HFlz=NaN*ones(ni,nj);
HOlz=NaN*ones(ni,nj);
HOmz90=NaN*ones(ni,nj);
HOlz90=NaN*ones(ni,nj);

HFmz(ID)=Hmz_smfrac;
HFlz(ID)=Hlz_smfrac;

HOmz(ID)=Hmz_stover;
HOlz(ID)=Hlz_stover;

HOmz90(ID)=Hmz_stover90;
HOlz90(ID)=Hlz_stover90;

bpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/ESM2M_Hist_Fore/'];
save([bpath 'Hist_Fore_' harv '_ts_map_zoop_overcon.mat'],...
    'Hmz_ttover','Hlz_ttover','HFmz','HOmz','HFlz','HOlz','-append');

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
surfm(geolat_t,geolon_t,HFmz)
colormap(cmBP2)
%cmocean('oxy')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1.2]);
colorbar
%colorbar('Position',[0.385 0.695 0.025 0.275],'orientation','vertical','AxisLocation','out')
%set(gcf,'renderer','painters')
text(0,1.75,'Mean fraction MZ hploss consumed','HorizontalAlignment','center')
text(-2.75,1.75,'A')

%2 - m over
subplot('Position',[0.475 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,HOmz)
colormap(cmBP)
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
surfm(geolat_t,geolon_t,HFlz)
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
surfm(geolat_t,geolon_t,HOlz)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
%colorbar('Position',[0.825 0.075 0.025 0.275],'orientation','vertical','AxisLocation','out')
%set(gcf,'renderer','painters')
text(0,1.75,'Mean times MZ overconsumed','HorizontalAlignment','center')
text(-2.75,1.75,'D')
%print('-dpng',[ppath 'Hist_' harv '_global_zoop_overcon.png'])

%% 2 plots of both maps
figure(3)
%1 - m frac
subplot('Position',[0.01 0.5 0.9 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,HFmz)
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
surfm(geolat_t,geolon_t,HFlz)
colormap(cmBP2)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1.1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'Mean fraction LZ hploss consumed','HorizontalAlignment','center')
%text(-2.75,1.75,'C')
print('-dpng',[ppath 'Hist_' harv '_global_zoop_fraccon.png'])

figure(4)
%2 - m over
subplot('Position',[0.01 0.5 0.9 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,HOmz)
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
surfm(geolat_t,geolon_t,HOlz)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'Mean times LZ overconsumed','HorizontalAlignment','center')
%text(-2.75,1.75,'D')
print('-dpng',[ppath 'Hist_' harv '_global_zoop_overcon.png'])

%% 1990 frac grid cells overcon
figure(5)
%2 - m over
subplot('Position',[0.01 0.5 0.9 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,HOmz90)
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
surfm(geolat_t,geolon_t,HOlz90)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'Mean times LZ overconsumed','HorizontalAlignment','center')
%text(-2.75,1.75,'D')
print('-dpng',[ppath 'Hist_1990_' harv '_global_zoop_overcon.png'])







