% Compare POEM to J&C15

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70_cmax-metab20_fcrit20_D075_J100_A050_Sm025_nmort5_BE05_CC050_lgRE00100_mdRE00400';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

%% Hist pristine figures -------------------------------------------
load([fpath 'Means_spinup_' cfile '.mat']);

% Pick which time period mean
%1956-2005
sp_smean=sp_mean;
sf_smean=sf_mean;
sd_smean=sd_mean;
mp_smean=mp_mean;
mf_smean=mf_mean;
md_smean=md_mean;
lp_smean=lp_mean;
ld_smean=ld_mean;
b_smean=b_mean;


%% Plots in space
[ni,nj]=size(geolon_t);

%fix lon shift
id=find(grid(:,2)<-180);
grid(id,2)=grid(id,2)+360;

x=-180:180;
y=-90:90;
[X,Y]=meshgrid(x,y);

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
title('Spinup log10 mean biomass All Fishes (g m^-^2)')
print('-dpng',[ppath 'Spinup_global_All_JC15.png'])


