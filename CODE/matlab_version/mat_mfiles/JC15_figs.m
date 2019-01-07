% Compare POEM to J&C15

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70_cmax-metab20_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00400';
harv = '03';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

%% Hist pristine figures -------------------------------------------
load([fpath 'Means_hist_pristine_' cfile '.mat']);

% Pick which time period mean
%1956-2005
sp_smean=sp_mean50;
sf_smean=sf_mean50;
sd_smean=sd_mean50;
mp_smean=mp_mean50;
mf_smean=mf_mean50;
md_smean=md_mean50;
lp_smean=lp_mean50;
ld_smean=ld_mean50;
b_smean=b_mean50;

sp_sprod=sp_prod50;
sf_sprod=sf_prod50;
sd_sprod=sd_prod50;
mp_sprod=mp_prod50;
mf_sprod=mf_prod50;
md_sprod=md_prod50;
lp_sprod=lp_prod50;
ld_sprod=ld_prod50;

%% Plots in space
[ni,nj]=size(geolon_t);

%fix lon shift
id=find(grid(:,2)<-180);
grid(id,2)=grid(id,2)+360;

x=-180:180;
y=-90:90;
[X,Y]=meshgrid(x,y);

all_histf = sf_smean + sp_smean + sd_smean + ...
    mf_smean + mp_smean + md_smean + ...
    lp_smean + ld_smean;

%Global grid
All_histf_gg = griddata(grid(:,2),grid(:,3),all_histf(:,1),X,Y);

Zsf=griddata(grid(:,2),grid(:,3),sf_smean,X,Y);
Zsp=griddata(grid(:,2),grid(:,3),sp_smean,X,Y);
Zsd=griddata(grid(:,2),grid(:,3),sd_smean,X,Y);
Zmf=griddata(grid(:,2),grid(:,3),mf_smean,X,Y);
Zmp=griddata(grid(:,2),grid(:,3),mp_smean,X,Y);
Zmd=griddata(grid(:,2),grid(:,3),md_smean,X,Y);
Zlp=griddata(grid(:,2),grid(:,3),lp_smean,X,Y);
Zld=griddata(grid(:,2),grid(:,3),ld_smean,X,Y);

% Diff maps of all fish
All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
AllF = Zsf+Zmf;
AllP = Zsp+Zmp+Zlp;
AllD = Zsd+Zmd+Zld;
AllS = Zsp+Zsf+Zsd;
AllM = Zmp+Zmf+Zmd;
AllL = Zlp+Zld;
FracPD = AllP ./ (AllP+AllD);
FracPF = AllP ./ (AllP+AllF);
FracPFvD = (AllP+AllF) ./ (AllP+AllF+AllD);

% Production
Psf=griddata(grid(:,2),grid(:,3),sf_sprod*365,X,Y);
Psp=griddata(grid(:,2),grid(:,3),sp_sprod*365,X,Y);
Psd=griddata(grid(:,2),grid(:,3),sd_sprod*365,X,Y);
Pmf=griddata(grid(:,2),grid(:,3),mf_sprod*365,X,Y);
Pmp=griddata(grid(:,2),grid(:,3),mp_sprod*365,X,Y);
Pmd=griddata(grid(:,2),grid(:,3),md_sprod*365,X,Y);
Plp=griddata(grid(:,2),grid(:,3),lp_sprod*365,X,Y);
Pld=griddata(grid(:,2),grid(:,3),ld_sprod*365,X,Y);

Psp(Psp<=0)=NaN;
Psf(Psf<=0)=NaN;
Psd(Psd<=0)=NaN;
Pmp(Pmp<=0)=NaN;
Pmf(Pmf<=0)=NaN;
Pmd(Pmd<=0)=NaN;
Plp(Plp<=0)=NaN;
Pld(Pld<=0)=NaN;

PAll = Psp+Psf+Psd+Pmp+Pmf+Pmd+Plp+Pld;
PAllF = Psf+Pmf;
PAllP = Psp+Pmp+Plp;
PAllD = Psd+Pmd+Pld;
PAllS = Psp+Psf+Psd;
PAllM = Pmp+Pmf+Pmd;
PAllL = Plp+Pld;

%plot info
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%% ALL
% Same colorbar scale as Jennings & Collingridge
figure(21)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(All)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
colormap('jet') 
colorbar('h')
caxis([-1 2])
title('Historic pristine 1956-2005 log10 mean biomass All Fishes (g m^-^2)')
print('-dpng',[ppath 'Hist_pristine_global_All_JC15.png'])


figure(31)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(PAll)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
colormap('jet') 
colorbar('h')
caxis([-2 1.5])
title('Historic pristine 1956-2005 log10 mean production All Fishes (g m^-^2 yr^-^1)')
print('-dpng',[ppath 'Hist_pristine_prod_All_JC15.png'])

figure(41)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(PAll./All)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
colormap('jet') 
colorbar('h')
caxis([0.25 1.5])
title('Historic pristine 1956-2005 mean P:B All Fishes')
print('-dpng',[ppath 'Hist_pristine_pbratio_All_JC15.png'])

