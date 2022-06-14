% Visualize output of FEISTY on same colorbar scale as Watson et al. 2015
% CORE-forced
% Observed effort
% Saved as nc files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%Orig: 
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end

%load([fpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);
load([fpath 'Climatology/Means_Climatol_' harv '_' cfile '.mat']);

close all

% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
% plotminlon=-280;
% plotmaxlon=80;
plotminlon=0;
plotmaxlon=360;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

geolat_t=lat;
geolon_t=lon;

ocean=NaN*ones(ni,nj);
ocean(ID)=ones(size(sf_mean));

land=ones(ni,nj);
land(ID)=nan(size(sf_mean));

% colors
load('MyColormaps.mat')
cjet = colormap('jet');
cjet2 = cjet(8:16:end,:);
cjet3 = [1 1 1;...
    153/255 153/255 255/255;...
    cjet2];

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

Zsf(ID)=sf_mean;
Zsp(ID)=sp_mean;
Zsd(ID)=sd_mean;
Zmf(ID)=mf_mean;
Zmp(ID)=mp_mean;
Zmd(ID)=md_mean;
Zlp(ID)=lp_mean;
Zld(ID)=ld_mean;
Zb(ID)=b_mean;


%% lp
cg = colormap('gray');
cg = flipud(cg);

figure(1)
ax1 = axes;
surf(ax1,geolon_t,geolat_t,log10(Zlp)); view(2); hold on;
shading flat
xlim([0 360])
ylim([-90 90])
title('log10 mean Adult P biomass (g m^-^2)')
colormap(ax1,'jet')
colorbar(ax1,'h')
caxis(ax1,[-10 -1])
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_global_LP.png'])

%% Diff maps of all fish
All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
AllF = Zsf+Zmf;
AllP = Zsp+Zmp+Zlp;
AllD = Zsd+Zmd+Zld;

%% 
% All P
figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Climatology log10 mean All P (g m^-^2)')
stamp([harv '_' cfile])
%print('-dpng',[ppath 'Climatol_' harv '_global_AllP.png'])

%% LP Watson scale
%divide by 100 m to convert g/m2 to g/m3
figure(3)
% axesm ('giso','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
axesm ('eqdcylin','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1);%,'origin',[0 -100 0])
% axesm ('miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zlp/100))
colormap(cjet3)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 -1]);
colorbar('h');
set(gcf,'renderer','painters')
title('Climatology log_1_0 mean LP (g m^-^3)')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_global_LP_Watson.png'])

%% All P Watson scale
%divide by 100 m to convert g/m2 to g/m3
figure(4)
% axesm ('giso','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
axesm ('eqdcylin','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)%,'origin',[0 -100 0])
% axesm ('miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP/100))
colormap(cjet3)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 -1]);
colorbar('h');
set(gcf,'renderer','painters')
title('Climatology log10 mean P (g m^-^3)')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_global_AllP_Watson.png'])

