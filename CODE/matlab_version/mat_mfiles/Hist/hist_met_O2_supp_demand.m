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

load([fpath 'Means_Historic_' harv '_met_' cfile '.mat']);

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
cmBP=cbrewer('seq','BuPu',50,'PCHIP');

%% WOA data
gpath='/Volumes/MIP/Obs_Data/WOA/O2/';
tpath='/Volumes/MIP/Obs_Data/WOA/Temp/';
load([gpath 'woa18_all_o00_01_m100.mat']); %Climatology in ug kg-1
load([tpath 'woa18_decav_t00_01_m100.mat']); %Climatology

% ug/kg --> mol/L
%1 kg = 1 L water
%10^-6 ug in 1 g
%32.00 g O2 in 1 mol O2

%O2 umol/kg --> mol/L
%10^6 umol in 1 mol
o2_100 = o_100 * 10e-6;
%Calc pressure (at sea surface, pO2=0.21)
Rgas = 0.082057; %units = L atm K−1 mol−1
T = t_100+273;
Tref = 15+273;
kb = 8.6173e-5; %Boltzmann const
pO2 = o2_100 .* T .* Rgas; %ambient O2 pressure (we use units of atm).

%% plot
[mlat,mlon] = meshgrid(double(lat),double(lon));
load coastlines;                     %decent looking coastlines
cmYR=cbrewer('seq','YlOrRd',45,'PCHIP');

figure(18) 
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(mlat,mlon,o2_100)
colormap(cmYR)                
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([-2 2]);
colorbar('Position',[0.05 0.55 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('O_2 (mol L^-^1)')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(mlat,mlon,pO2)
colormap(cmYR)             
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([-2 2]);
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('pO_2 (atm)')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(mlat,mlon,t_100)
cmocean('thermal')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 30]);
colorbar('Position',[0.55 0.55 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Temp (^oC)')

% All
%subplot('Position',[0.5 0 0.5 0.5])

print('-dpng',[ppath 'WOA_global_O2_pO2_Temp_100m.png'])

%% Calc example supply
load('MI_fntypes3_spp.mat','FBname','alphaS','del','Esup')
DalphaS = alphaS(1); %Atl cod in μmol O2 g−3/4 h−1 atm−1
FalphaS = alphaS(2); %Lanternfish
PalphaS = alphaS(3); %Striped bass
Ddel = del(1);
Fdel = del(2);
Pdel = del(3);
DEsup = Esup(1);
FEsup = Esup(2);
PEsup = Esup(3);

%convert mass to dry g?
M_s = (10^((log10(0.001)+log10(0.5))/2)) /9;
M_m = (10^((log10(0.5)+log10(250))/2)) /9;
M_l = (10^((log10(250)+log10(125000))/2)) /9;

%resulting units in μmol O2
supp_sf = FalphaS .* exp(-FEsup./kb .* ((1./T) - (1./Tref))) .* M_s.^Fdel .* pO2;
supp_sp = PalphaS .* exp(-PEsup./kb .* ((1./T) - (1./Tref))) .* M_s.^Pdel .* pO2;
supp_sd = DalphaS .* exp(-DEsup./kb .* ((1./T) - (1./Tref))) .* M_s.^Ddel .* pO2;
supp_mf = FalphaS .* exp(-FEsup./kb .* ((1./T) - (1./Tref))) .* M_m.^Fdel .* pO2;
supp_mp = PalphaS .* exp(-PEsup./kb .* ((1./T) - (1./Tref))) .* M_m.^Pdel .* pO2;
supp_md = DalphaS .* exp(-DEsup./kb .* ((1./T) - (1./Tref))) .* M_m.^Ddel .* pO2;
supp_lp = PalphaS .* exp(-PEsup./kb .* ((1./T) - (1./Tref))) .* M_l.^Pdel .* pO2;
supp_ld = DalphaS .* exp(-DEsup./kb .* ((1./T) - (1./Tref))) .* M_l.^Ddel .* pO2;

%% FEISTY demand units
%1 kg = 1 L water
%1 L = 10^-3 m^3
%10^-6 ug in 1 g
%32.00 g O2 in 1 mol O2
%106/138 mol C in 1 mol O2 (Redfield ratio)
%12.01 g C in 1 mol C
%1e6 μmol in 1 mol
%1 g dry W in 9 g wet

%g WW C --> μmol O2
%1 g dry in 9 g wet; 12.01 g C in 1 mol C; 106/138 mol C in 1 mol O2; 1e6 μmol in 1 mol
Omet_sf = Psf * (1/9.0) * (1/12.01) * (106/138) * 1e6;
Omet_sp = Psp * (1/9.0) * (1/12.01) * (106/138) * 1e6;
Omet_sd = Psd * (1/9.0) * (1/12.01) * (106/138) * 1e6;
Omet_mf = Pmf * (1/9.0) * (1/12.01) * (106/138) * 1e6;
Omet_mp = Pmp * (1/9.0) * (1/12.01) * (106/138) * 1e6;
Omet_md = Pmd * (1/9.0) * (1/12.01) * (106/138) * 1e6;
Omet_lp = Plp * (1/9.0) * (1/12.01) * (106/138) * 1e6;
Omet_ld = Pld * (1/9.0) * (1/12.01) * (106/138) * 1e6;
Omet_S = AllS * (1/9.0) * (1/12.01) * (106/138) * 1e6;
Omet_M = AllM * (1/9.0) * (1/12.01) * (106/138) * 1e6;
Omet_L = AllL * (1/9.0) * (1/12.01) * (106/138) * 1e6;
Omet_F = AllF * (1/9.0) * (1/12.01) * (106/138) * 1e6;
Omet_P = AllP * (1/9.0) * (1/12.01) * (106/138) * 1e6;
Omet_D = AllD * (1/9.0) * (1/12.01) * (106/138) * 1e6;


%% 8plot by stage and type - Supply
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - sf
subplot('Position',[0.025 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(mlat,mlon,supp_sf)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([10 70])
set(gcf,'renderer','painters')
text(0,1.75,'SF','HorizontalAlignment','center')

%B - sp
subplot('Position',[0.025 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(mlat,mlon,supp_sp)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([10 70])
set(gcf,'renderer','painters')
text(0,1.75,'SP','HorizontalAlignment','center')

%C - mp
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(mlat,mlon,supp_mp)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([5 100])
set(gcf,'renderer','painters')
text(0,1.75,'MP','HorizontalAlignment','center')

%D - lp
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(mlat,mlon,supp_lp)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 100])
set(gcf,'renderer','painters')
text(0,1.75,'LP','HorizontalAlignment','center')

%E - mf
subplot('Position',[0.5 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(mlat,mlon,supp_mf)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([5 100])
set(gcf,'renderer','painters')
text(0,1.75,'MF','HorizontalAlignment','center')

%F - sd
subplot('Position',[0.5 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(mlat,mlon,supp_sd)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([10 70])
%colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'SD','HorizontalAlignment','center')

%G - md
subplot('Position',[0.5 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(mlat,mlon,supp_md)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([5 100])
set(gcf,'renderer','painters')
text(0,1.75,'MD','HorizontalAlignment','center')

%H - ld
subplot('Position',[0.5 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(mlat,mlon,supp_ld)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 100])
set(gcf,'renderer','painters')
text(0,1.75,'LD','HorizontalAlignment','center')

print('-dpng',[ppath 'WOA_O2supp_global_stages_6plot_YOR.png'])

%% 8plot by stage and type - Demand
f3 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - sf
subplot('Position',[0.025 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omet_sf)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([10 70])
set(gcf,'renderer','painters')
text(0,1.75,'SF','HorizontalAlignment','center')

%B - sp
subplot('Position',[0.025 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omet_sp)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([10 70])
set(gcf,'renderer','painters')
text(0,1.75,'SP','HorizontalAlignment','center')

%C - mp
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omet_mp)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([5 100])
set(gcf,'renderer','painters')
text(0,1.75,'MP','HorizontalAlignment','center')

%D - lp
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omet_lp)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 100])
set(gcf,'renderer','painters')
text(0,1.75,'LP','HorizontalAlignment','center')

%E - mf
subplot('Position',[0.5 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omet_mf)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([5 100])
set(gcf,'renderer','painters')
text(0,1.75,'MF','HorizontalAlignment','center')

%F - sd
subplot('Position',[0.5 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omet_sd)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([10 70])
%colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'SD','HorizontalAlignment','center')

%G - md
subplot('Position',[0.5 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omet_md)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([5 100])
set(gcf,'renderer','painters')
text(0,1.75,'MD','HorizontalAlignment','center')

%H - ld
subplot('Position',[0.5 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omet_ld)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 100])
set(gcf,'renderer','painters')
text(0,1.75,'LD','HorizontalAlignment','center')

print('-dpng',[ppath 'Hist_O2met_' harv '_global_stages_6plot_YOR.png'])


