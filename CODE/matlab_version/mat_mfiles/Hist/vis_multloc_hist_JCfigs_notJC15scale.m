% Visualize output of POEM
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files
% Map same quantities and scales as J&C15

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

load([fpath 'Means_Historic_',harv,'_' cfile '.mat'],...
    'sf_mean5','sp_mean5','sd_mean5',...
    'mf_mean5','mp_mean5','md_mean5',...
     'b_mean5','lp_mean5','ld_mean5');
load([fpath 'Means_Historic_',harv,'_prod_' cfile '.mat'],...
    'sf_ptot','sp_ptot','sd_ptot',...
    'mf_ptot','mp_ptot','md_ptot',...
    'lp_ptot','ld_ptot');


% plot info
[ni,nj]=size(geolon_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

geolat_t=double(geolat_t);
geolon_t=double(geolon_t);

cmS=cbrewer('div','Spectral',11);
cmS2=cmS(6:end,:);
cmS2=flipud(cmS2);
cmBP=cbrewer('seq','BuPu',50);

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
Zsf(ID)=sf_mean5;
Zsp(ID)=sp_mean5;
Zsd(ID)=sd_mean5;
Zmf(ID)=mf_mean5;
Zmp(ID)=mp_mean5;
Zmd(ID)=md_mean5;
Zlp(ID)=lp_mean5;
Zld(ID)=ld_mean5;
Zb(ID)=b_mean5;

Psf=NaN*ones(ni,nj);
Psp=NaN*ones(ni,nj);
Psd=NaN*ones(ni,nj);
Pmf=NaN*ones(ni,nj);
Pmp=NaN*ones(ni,nj);
Pmd=NaN*ones(ni,nj);
Plp=NaN*ones(ni,nj);
Pld=NaN*ones(ni,nj);

%sum over year to be per year instead of per day
Psf(ID)=sf_ptot;
Psp(ID)=sp_ptot;
Psd(ID)=sd_ptot;
Pmf(ID)=mf_ptot;
Pmp(ID)=mp_ptot;
Pmd(ID)=md_ptot;
Plp(ID)=lp_ptot;
Pld(ID)=ld_ptot;

ocean=NaN*ones(ni,nj);
ocean(ID)=ones(size(sf_mean5));


%% Diff maps of all fish
%JC 1-10^6g
%JC 10^2-10^4g
%Production should be per day

All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;

PAll = Psp+Psf+Psd+Pmp+Pmf+Pmd+Plp+Pld;
PAll(PAll<=0) = eps;

%Prod:Biom
PB=PAll./All;

%% Jet
%All Biom
figure(4)
subplot('Position',[0 0.67 1 0.32])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 3]);
colorbar('Position',[0.69 0.71 0.025 0.25],'orientation','vertical')
text(0.2,1.75,'log_1_0 Biomass All fishes (g m^-^2)',...
    'HorizontalAlignment','center');
set(gcf,'renderer','painters')

% All Prod
subplot('Position',[0 0.34 1 0.32])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(PAll))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 3]);
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
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 8]);
colorbar('Position',[0.69 0.05 0.025 0.25],'orientation','vertical')
text(0.2,1.75,'Production:Biomass',...
    'HorizontalAlignment','center');
set(gcf,'renderer','painters')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Hist9095_' harv '_global_All_notJCscale.png'])

%% Colorbrewer
%All Biom
figure(5)
subplot('Position',[0 0.67 1 0.32])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 3]);
colorbar('Position',[0.69 0.71 0.025 0.25],'orientation','vertical')
text(0.2,1.75,'log_1_0 Biomass All fishes (g m^-^2)',...
    'HorizontalAlignment','center');
set(gcf,'renderer','painters')

% All Prod
subplot('Position',[0 0.34 1 0.32])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(PAll+eps))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 3]);
colorbar('Position',[0.69 0.38 0.025 0.25],'orientation','vertical')
text(0.2,1.75,'log_1_0 Prodution All fishes (g m^-^2 y^-^1)',...
    'HorizontalAlignment','center');
set(gcf,'renderer','painters')

% P:B
subplot('Position',[0 0.01 1 0.32])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,PB)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 8]);
colorbar('Position',[0.69 0.05 0.025 0.25],'orientation','vertical')
text(0.2,1.75,'Production:Biomass',...
    'HorizontalAlignment','center');
set(gcf,'renderer','painters')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Hist9095_' harv '_global_All_notJCscale_BP.png'])


