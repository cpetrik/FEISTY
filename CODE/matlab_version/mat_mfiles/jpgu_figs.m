% Visualize output of POEM

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

% Zsf=NaN*ones(ni,nj);
% Zsp=NaN*ones(ni,nj);
% Zsd=NaN*ones(ni,nj);
% Zmf=NaN*ones(ni,nj);
% Zmp=NaN*ones(ni,nj);
% Zmd=NaN*ones(ni,nj);
% Zlp=NaN*ones(ni,nj);
% Zld=NaN*ones(ni,nj);
% Zb=NaN*ones(ni,nj);
% 
% Zsf(grid(:,1))=sf_smean;
% Zsp(grid(:,1))=sp_smean;
% Zsd(grid(:,1))=sd_smean;
% Zmf(grid(:,1))=mf_smean;
% Zmp(grid(:,1))=mp_smean;
% Zmd(grid(:,1))=md_smean;
% Zlp(grid(:,1))=lp_smean;
% Zld(grid(:,1))=ld_smean;
% Zb(grid(:,1))=b_smean;
% 
% % Diff maps of all fish
% All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
% AllF = Zsf+Zmf;
% AllP = Zsp+Zmp+Zlp;
% AllD = Zsd+Zmd+Zld;
% AllS = Zsp+Zsf+Zsd;
% AllM = Zmp+Zmf+Zmd;
% AllL = Zlp+Zld;
% FracPD = AllP ./ (AllP+AllD);
% FracPF = AllP ./ (AllP+AllF);
% FracPFvD = (AllP+AllF) ./ (AllP+AllF+AllD);
% 
% % Production
% Psf=NaN*ones(ni,nj);
% Psp=NaN*ones(ni,nj);
% Psd=NaN*ones(ni,nj);
% Pmf=NaN*ones(ni,nj);
% Pmp=NaN*ones(ni,nj);
% Pmd=NaN*ones(ni,nj);
% Plp=NaN*ones(ni,nj);
% Pld=NaN*ones(ni,nj);
% 
% Psf(grid(:,1))=sf_sprod*365;
% Psp(grid(:,1))=sp_sprod*365;
% Psd(grid(:,1))=sd_sprod*365;
% Pmf(grid(:,1))=mf_sprod*365;
% Pmp(grid(:,1))=mp_sprod*365;
% Pmd(grid(:,1))=md_sprod*365;
% Plp(grid(:,1))=lp_sprod*365;
% Pld(grid(:,1))=ld_sprod*365;
% 
% Psp(Psp<=0)=NaN;
% Psf(Psf<=0)=NaN;
% Psd(Psd<=0)=NaN;
% Pmp(Pmp<=0)=NaN;
% Pmf(Pmf<=0)=NaN;
% Pmd(Pmd<=0)=NaN;
% Plp(Plp<=0)=NaN;
% Pld(Pld<=0)=NaN;
% 
% PAll = Psp+Psf+Psd+Pmp+Pmf+Pmd+Plp+Pld;
% PAllF = Psf+Pmf;
% PAllP = Psp+Pmp+Plp;
% PAllD = Psd+Pmd+Pld;
% PAllS = Psp+Psf+Psd;
% PAllM = Pmp+Pmf+Pmd;
% PAllL = Plp+Pld;

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
% figure(21)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-280 80],'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat_t,geolon_t,log10(All))
% colormap('jet')              
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 2]);
% hcb = colorbar('h');
% ylim(hcb,[-1 2])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Historic pristine 1956-2005 log10 mean biomass All Fishes (g m^-^2)')
% print('-dpng',[ppath 'Hist_pristine_global_All_jpgu.png'])

figure(21)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(All_histf_gg)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
colormap('jet') 
colorbar('h')
caxis([-1 2])
title('Historic pristine 1956-2005 log10 mean biomass All Fishes (g m^-^2)')
print('-dpng',[ppath 'Hist_pristine_global_All_JC15.png'])

%%
figure(31)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-280 80],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,log10(PAll))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1.5]);
hcb = colorbar('h');
ylim(hcb,[-2 1.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic pristine 1956-2005 log10 mean production All Fishes (g m^-^2 yr^-^1)')
print('-dpng',[ppath 'Hist_pristine_prod_All_jpgu.png'])

figure(41)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-280 80],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,PAll./All)
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.25 1.5]);
hcb = colorbar('h');
ylim(hcb,[0.25 1.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic pristine 1956-2005 mean P:B All Fishes')
print('-dpng',[ppath 'Hist_pristine_pbratio_All_jpgu.png'])

%% all F
figure(22)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-280 80],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,log10(AllF))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic pristine 1956-2005 log10 mean biomass All F (g m^-^2)')
print('-dpng',[ppath 'Hist_pristine',harv,'_global_AllF.png'])

% all D
figure(23)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-280 80],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,log10(AllD))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic pristine 1956-2005 log10 mean biomass All D (g m^-^2)')
print('-dpng',[ppath 'Hist_pristine',harv,'_global_AllD.png'])

% All P
figure(24)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-280 80],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,log10(AllP))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic pristine 1956-2005 log10 mean biomass All P (g m^-^2)')
print('-dpng',[ppath 'Hist_pristine',harv,'_global_AllP.png'])


%% Change in catch potential figure --------------------------------

%Historic
load([fpath 'Means_hist_fished',harv,'_' cfile '.mat']);

Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);

Cmf(grid(:,1))=mf_my50;
Cmp(grid(:,1))=mp_my50;
Cmd(grid(:,1))=md_my50;
Clp(grid(:,1))=lp_my50;
Cld(grid(:,1))=ld_my50;

hCAll = Cmp+Cmf+Cmd+Clp+Cld;
hCAllF = Cmf;
hCAllP = Cmp+Clp;
hCAllD = Cmd+Cld;
hCAllM = Cmp+Cmf+Cmd;
hCAllL = Clp+Cld;


%Forecast
load([fpath 'Means_fore_fished',harv,'_' cfile '.mat']);

Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);

Cmf(grid(:,1))=mf_my50;
Cmp(grid(:,1))=mp_my50;
Cmd(grid(:,1))=md_my50;
Clp(grid(:,1))=lp_my50;
Cld(grid(:,1))=ld_my50;

fCAll = Cmp+Cmf+Cmd+Clp+Cld;
fCAllF = Cmf;
fCAllP = Cmp+Clp;
fCAllD = Cmd+Cld;
fCAllM = Cmp+Cmf+Cmd;
fCAllL = Clp+Cld;

dCAll = 100 * (fCAll - hCAll) ./ hCAll;

%% All

rjet = colormap(jet);
rjet = flipud(rjet);

figure(48)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,dCAll)
colormap(rjet)              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 50]);
hcb = colorbar('h');
ylim(hcb,[-100 50])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference from 1956-2005 mean catch All Fishes')
print('-dpng',[ppath 'Fished',harv, '_global_All_catch_foref_histf.png'])
