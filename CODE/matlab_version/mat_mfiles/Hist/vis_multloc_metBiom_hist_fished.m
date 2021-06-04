% Visualize output of FEISTY
% Historic time period (1861-2005) at all locations
% Metab*Biomass
% Use stoich to calc equiv O2 demand

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/Historic_ESM2M/'];
%fpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

ppath = [pp cfile '/'];

load([fpath 'Means_Historic_' harv '_metBiom_' cfile '.mat']);

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

%% Pick which time period mean
% 1951-2000
sp_smean=sp_met50;
sf_smean=sf_met50;
sd_smean=sd_met50;
mp_smean=mp_met50;
mf_smean=mf_met50;
md_smean=md_met50;
lp_smean=lp_met50;
ld_smean=ld_met50;

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

Psf(ID)=sf_met50;
Psp(ID)=sp_met50;
Psd(ID)=sd_met50;
Pmf(ID)=mf_met50;
Pmp(ID)=mp_met50;
Pmd(ID)=md_met50;
Plp(ID)=lp_met50;
Pld(ID)=ld_met50;

Psf(Psf(:)<0) = 0;
Psp(Psp(:)<0) = 0;
Psd(Psd(:)<0) = 0;
Pmf(Pmf(:)<0) = 0;
Pmp(Pmp(:)<0) = 0;
Pmd(Pmd(:)<0) = 0;
Plp(Plp(:)<0) = 0;
Pld(Pld(:)<0) = 0;

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
load coastlines;                     %decent looking coastlines
cmBP=cbrewer('seq','YlOrRd',45,'PCHIP');

%% WOA O2 data
gpath='/Volumes/MIP/Obs_Data/WOA/O2/';
load([gpath 'woa18_all_o00_01_m100.mat']); %Climatology in ug kg-1

[mlat,mlon] = meshgrid(double(lat),double(lon));

% ug/kg --> mol/L
%1 kg = 1 L water
%10^-6 ug in 1 g
%32.00 g O2 in 1 mol O2

%O2 umol/kg --> mol/L
%10^6 umol in 1 mol
o2_100 = o_100 * 10e-6;

%% FEISTY units
%1 kg = 1 L water
%1 L = 10^-3 m^3
%10^-6 ug in 1 g
%32.00 g O2 in 1 mol O2
%106/138 mol C in 1 mol O2 (Redfield ratio)
%12.01 g C in 1 mol C
%1 g dry W in 9 g wet

%g WW --> mol/L
%1 g dry in 9 g wet; 12.01 g C in 1 mol C; 100 m depth -> m^3; 10^-3 m^3 in 1 L; 106/138 mol C in 1 mol O2
Omet_sf = Psf * (1/9.0) * (1/12.01) * (1/100) * 1e-3 * (106/138);
Omet_sp = Psp * (1/9.0) * (1/12.01) * (1/100) * 1e-3 * (106/138);
Omet_sd = Psd * (1/9.0) * (1/12.01) * (1/100) * 1e-3 * (106/138);
Omet_mf = Pmf * (1/9.0) * (1/12.01) * (1/100) * 1e-3 * (106/138);
Omet_mp = Pmp * (1/9.0) * (1/12.01) * (1/100) * 1e-3 * (106/138);
Omet_md = Pmd * (1/9.0) * (1/12.01) * (1/100) * 1e-3 * (106/138);
Omet_lp = Plp * (1/9.0) * (1/12.01) * (1/100) * 1e-3 * (106/138);
Omet_ld = Pld * (1/9.0) * (1/12.01) * (1/100) * 1e-3 * (106/138);
Omet_S = AllS * (1/9.0) * (1/12.01) * (1/100) * 1e-3 * (106/138);
Omet_M = AllM * (1/9.0) * (1/12.01) * (1/100) * 1e-3 * (106/138);
Omet_L = AllL * (1/9.0) * (1/12.01) * (1/100) * 1e-3 * (106/138);
Omet_F = AllF * (1/9.0) * (1/12.01) * (1/100) * 1e-3 * (106/138);
Omet_P = AllP * (1/9.0) * (1/12.01) * (1/100) * 1e-3 * (106/138);
Omet_D = AllD * (1/9.0) * (1/12.01) * (1/100) * 1e-3 * (106/138);

%% Map of O2
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(mlat,mlon,(o2_100))
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2e-5 12e-5]);
colorbar('orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
title('WOA O_2 100 m mean (mol L^-^1)')
print('-dpng',[ppath 'WOA_O2_global_molL_YOR.png'])

%% 6 plot of types
figure(2)
%A - Small
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Omet_S))
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([2e-5 12e-5]);
colorbar('Position',[0.385 0.695 0.025 0.275],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Small','HorizontalAlignment','center')

%D - forage
subplot('Position',[0.45 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Omet_F))
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([2e-5 12e-5]);
colorbar('Position',[0.825 0.695 0.025 0.275],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Forage','HorizontalAlignment','center')

%B - med
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Omet_M))
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([2e-5 12e-5]);
colorbar('Position',[0.385 0.385 0.025 0.275],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Medium','HorizontalAlignment','center')

%E - P
subplot('Position',[0.45 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Omet_P))
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([2e-5 12e-5]);
colorbar('Position',[0.825 0.385 0.025 0.275],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Large pelagic','HorizontalAlignment','center')

%C - large
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Omet_L))
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([2e-5 12e-5]);
colorbar('Position',[0.385 0.075 0.025 0.275],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Large','HorizontalAlignment','center')

%F - dem
subplot('Position',[0.45 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Omet_D))
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([2e-5 12e-5]);
colorbar('Position',[0.825 0.075 0.025 0.275],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Demersal','HorizontalAlignment','center')

print('-dpng',[ppath 'Hist_O2_metBiom_' harv '_global_types_6plot_YOR.png'])

%% 8plot by stage and type
% Mollweid proj
f3 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - sf
subplot('Position',[0.025 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omet_sf)
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2e-9 12e-9]);
set(gcf,'renderer','painters')
text(0,1.75,'SF','HorizontalAlignment','center')

%B - sp
subplot('Position',[0.025 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omet_sp)
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2e-9 12e-9]);
set(gcf,'renderer','painters')
text(0,1.75,'SP','HorizontalAlignment','center')

%C - mp
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omet_mp)
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2e-9 12e-9]);
set(gcf,'renderer','painters')
text(0,1.75,'MP','HorizontalAlignment','center')

%D - lp
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omet_lp)
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2e-9 12e-9]);
set(gcf,'renderer','painters')
text(0,1.75,'LP','HorizontalAlignment','center')

%E - mf
subplot('Position',[0.475 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omet_mf)
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2e-9 12e-9]);
set(gcf,'renderer','painters')
text(0,1.75,'MF','HorizontalAlignment','center')

%F - sd
subplot('Position',[0.475 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omet_sd)
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2e-9 12e-9]);
colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'SD','HorizontalAlignment','center')

%G - md
subplot('Position',[0.475 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omet_md)
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2e-9 12e-9]);
set(gcf,'renderer','painters')
text(0,1.75,'MD','HorizontalAlignment','center')

%H - ld
subplot('Position',[0.475 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omet_ld)
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2e-9 12e-9]);
set(gcf,'renderer','painters')
text(0,1.75,'LD','HorizontalAlignment','center')

print('-dpng',[ppath 'Hist_O2_metBiom_' harv '_global_stages_6plot_YOR.png'])
