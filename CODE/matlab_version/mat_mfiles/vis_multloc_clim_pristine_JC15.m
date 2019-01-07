% Visualize output of POEM
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files
% Map same quantities and scales as J&C15

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end

load([fpath 'Means_bio_prod_Climatol_pristine_' cfile '.mat']);

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

close all

% plot info
[ni,nj]=size(lon);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

geolat_t=lat;
geolon_t=lon;

cmS=cbrewer('div','Spectral',11);
cmS2=cmS(6:end,:);
cmS2=flipud(cmS2);

%% Plots in space
Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zb=NaN*ones(ni,nj);

%mean
Zsf(ID)=sf_mean;
Zsp(ID)=sp_mean;
Zsd(ID)=sd_mean;
Zmf(ID)=mf_mean;
Zmp(ID)=mp_mean;
Zmd(ID)=md_mean;
Zlp(ID)=lp_mean;
Zld(ID)=ld_mean;
Zb(ID)=b_mean;

Psf=NaN*ones(ni,nj);
Psp=NaN*ones(ni,nj);
Psd=NaN*ones(ni,nj);
Pmf=NaN*ones(ni,nj);
Pmp=NaN*ones(ni,nj);
Pmd=NaN*ones(ni,nj);
Plp=NaN*ones(ni,nj);
Pld=NaN*ones(ni,nj);

%sum over year to be per year instead of per day
Psf(ID)=sf_tprod;
Psp(ID)=sp_tprod;
Psd(ID)=sd_tprod;
Pmf(ID)=mf_tprod;
Pmp(ID)=mp_tprod;
Pmd(ID)=md_tprod;
Plp(ID)=lp_tprod;
Pld(ID)=ld_tprod;

ocean=NaN*ones(ni,nj);
ocean(ID)=ones(size(sf_mean));


%% Diff maps of all fish
%JC 1-10^6g
%JC 10^2-10^4g
%Production should be per day

All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;

PAll = Psp+Psf+Psd+Pmp+Pmf+Pmd+Plp+Pld;

%Prod:Biom
PB=PAll./All;

%% All Biom
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 2]);
hcb = colorbar('h');
ylim(hcb,[-1 2])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology log10 Biomass All fishes (g m^-^2)')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_pristine_global_All_biomass_JCscale.png'])

%% All Prod
figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(PAll+eps))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1.5]);
hcb = colorbar('h');
ylim(hcb,[-2 1.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology log10 Prodution All fishes (g m^-^2 y^-^1)')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_pristine_global_All_prod_JCscale.png'])

% P:B
figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,PB)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.25 1.5]);
hcb = colorbar('h');
ylim(hcb,[0.25 1.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology Production:Biomass')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_pristine_global_All_PBratio_JCscale.png'])

%% All Biom
figure(4)
subplot('Position',[0 0.67 1 0.32])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 2]);
colorbar('Position',[0.69 0.71 0.025 0.25],'orientation','vertical')
text(0.2,1.75,'log_1_0 Biomass All fishes (g m^-^2)',...
    'HorizontalAlignment','center');
set(gcf,'renderer','painters')

% All Prod
subplot('Position',[0 0.34 1 0.32])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(PAll+eps))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1.5]);
colorbar('Position',[0.69 0.38 0.025 0.25],'orientation','vertical')
text(0.2,1.75,'log_1_0 Prodution All fishes (g m^-^2 y^-^1)',...
    'HorizontalAlignment','center');
set(gcf,'renderer','painters')

% P:B
subplot('Position',[0 0.01 1 0.32])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,PB)
colormap('jet')
%colormap(cmS2)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.25 1.5]);
colorbar('Position',[0.69 0.05 0.025 0.25],'orientation','vertical')
text(0.2,1.75,'Production:Biomass',...
    'HorizontalAlignment','center');
set(gcf,'renderer','painters')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_pristine_global_All_JCscale.png'])

%% All Biom
figure(5)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',[-180 180])
surfm(geolat_t,geolon_t,log10(All));
%xlim([-180 180])
%ylim([-90 90])
colormap('jet')
caxis([-1 2]);
colorbar('h')
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('Climatology log10 Biomass All fishes (g m^-^2)')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_pristine_global_All_biomass_JCscale_axes.png'])

