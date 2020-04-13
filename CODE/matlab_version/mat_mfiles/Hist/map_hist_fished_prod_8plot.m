% Visualize output of FEISTY
% Historic time period (1861-2005) at all locations
% 145 years
% Saved as mat files
% Production

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

%fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
fpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

ppath = [pp cfile '/'];

load([fpath 'Means_Historic_' harv '_prod_' cfile '.mat']);

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

%% Pick which time period mean
% 1951-2000
sp_smean=sp_prod50;
sf_smean=sf_prod50;
sd_smean=sd_prod50;
mp_smean=mp_prod50;
mf_smean=mf_prod50;
md_smean=md_prod50;
lp_smean=lp_prod50;
ld_smean=ld_prod50;
b_smean=b_mean50;

%% Plots in space
[ni,nj]=size(geolon_t);
ID = grid(:,1);

Psf=NaN*ones(ni,nj);
Psp=NaN*ones(ni,nj);
Psd=NaN*ones(ni,nj);
Pmf=NaN*ones(ni,nj);
Pmp=NaN*ones(ni,nj);
Pmd=NaN*ones(ni,nj);
Plp=NaN*ones(ni,nj);
Pld=NaN*ones(ni,nj);
Pb=NaN*ones(ni,nj);

Psf(ID)=sf_prod50;
Psp(ID)=sp_prod50;
Psd(ID)=sd_prod50;
Pmf(ID)=mf_prod50;
Pmp(ID)=mp_prod50;
Pmd(ID)=md_prod50;
Plp(ID)=lp_prod50;
Pld(ID)=ld_prod50;
Pb(ID)=b_mean50;

Psf(Psf(:)<0) = 0;
Psp(Psp(:)<0) = 0;
Psd(Psd(:)<0) = 0;
Pmf(Pmf(:)<0) = 0;
Pmp(Pmp(:)<0) = 0;
Pmd(Pmd(:)<0) = 0;
Plp(Plp(:)<0) = 0;
Pld(Pld(:)<0) = 0;
Pb(Pb(:)<0) = 0;

All = Psp+Psf+Psd+Pmp+Pmf+Pmd+Plp+Pld;
AllF = Psf+Pmf;
AllP = Psp+Pmp+Plp;
AllD = Psd+Pmd+Pld;
AllS = Psp+Psf+Psd;
AllM = Pmp+Pmf+Pmd;
AllL = Plp+Pld;

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

cmBP=cbrewer('seq','BuPu',50,'PCHIP');

%% COBALT
gpath='/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
load([gpath 'cobalt_zoop_biom_means.mat'],'mzprod_mean_hist','lzprod_mean_hist'); 
load([gpath 'cobalt_det_biom_means.mat'],'det_mean_hist');
load([gpath 'cobalt_npp_means.mat'],'npp_mean_hist');

% Zoop and det and npp 
%ESM2M in mmol N m-2 or mmol N m-2 d-1
%molN/m2 --> g/m2
%106/16 mol C in 1 mol N
%12.01 g C in 1 mol C
%1 g dry W in 9 g wet W
mmz_prod = mzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mlz_prod = lzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_hist = det_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_hist = npp_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

zprod = mmz_prod + mlz_prod;

%% calculate: 
%benthic production = benthic biomass * detritus flux * benthic efficiency
prodB = Pb .* det_hist .* 0.075;

%%
edges = -5:3;

figure(10)
subplot(3,3,1)
histogram(log10(npp_hist(:)),edges)

subplot(3,3,2)
histogram(log10(zprod(:)),edges)

subplot(3,3,3)
histogram(log10(det_hist(:)),edges)

subplot(3,3,4)
histogram(log10(prodB(:)),edges)

subplot(3,3,5)
histogram(log10(AllF(:)),edges)

subplot(3,3,6)
histogram(log10(AllP(:)),edges)

subplot(3,3,7)
histogram(log10(AllD(:)),edges)

subplot(3,3,8)
histogram(log10(All(:)),edges)

%% figure info
f1 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - npp
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(npp_hist))
colormap(cmBP)
colorbar('Position',[0.43 0.76 0.025 0.225],'orientation','vertical','AxisLocation','out')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 1.5]);
set(gcf,'renderer','painters')
text(0,1.75,'Net primary','HorizontalAlignment','center')
text(-2.5,1.75,'A')

% - mesozoo
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(zprod))
colormap(cmBP)
colorbar('Position',[0.43 0.51 0.025 0.225],'orientation','vertical','AxisLocation','out')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(0,1.75,'Mesozoo','HorizontalAlignment','center')
text(-2.5,1.75,'B')

% - det
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(det_hist))
colormap(cmBP)
colorbar('Position',[0.43 0.26 0.025 0.225],'orientation','vertical','AxisLocation','out')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2.0 0.5]);
set(gcf,'renderer','painters')
text(0,1.75,'Detritus','HorizontalAlignment','center')
text(-2.5,1.75,'C')

% - benthos
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(prodB))
colormap(cmBP)
colorbar('Position',[0.43 0.01 0.025 0.225],'orientation','vertical','AxisLocation','out')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 0]);
set(gcf,'renderer','painters')
text(0,1.75,'Benthos','HorizontalAlignment','center')
text(-2.5,1.75,'D')

% - forage
subplot('Position',[0.525 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllF))
colormap(cmBP)
colorbar('Position',[0.93 0.76 0.025 0.225],'orientation','vertical','AxisLocation','out')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 0]);
set(gcf,'renderer','painters')
text(0,1.75,'Forage','HorizontalAlignment','center')
text(-2.5,1.75,'E')

% - LP
subplot('Position',[0.525 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 0]);
colorbar('Position',[0.93 0.51 0.025 0.225],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Large pelagic','HorizontalAlignment','center')
text(-2.5,1.75,'F')

%F - dem
subplot('Position',[0.525 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllD))
colormap(cmBP)
colorbar('Position',[0.93 0.26 0.025 0.225],'orientation','vertical','AxisLocation','out')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 0]);
set(gcf,'renderer','painters')
text(0,1.75,'Demersal','HorizontalAlignment','center')
text(-2.5,1.75,'G')

%- All
subplot('Position',[0.525 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All))
colormap(cmBP)
colorbar('Position',[0.93 0.01 0.025 0.225],'orientation','vertical','AxisLocation','out')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 0]);
set(gcf,'renderer','painters')
text(0,1.75,'All','HorizontalAlignment','center')
text(-2.5,1.75,'H')

print('-dpng',[pp 'Hist_',harv,'_global_prod_8plot.png'])

