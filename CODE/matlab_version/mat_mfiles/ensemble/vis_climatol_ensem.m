function vis_climatol_ensem(simname,sf_mean,sp_mean,sd_mean,...
    mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean)

% Visualize output of POEM
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files

Pdrpbx = '/Users/cpetrik/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

pp = [Pdrpbx 'Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/param_ensemble/Dc_D075_Sm025_nmort1_noCC_RE00100/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac


% colors
load('MyColormaps.mat')

close all

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

All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
AllF = Zsf+Zmf;
AllP = Zsp+Zmp+Zlp;
AllD = Zsd+Zmd+Zld;
AllS = Zsp+Zsf+Zsd;
AllM = Zmp+Zmf+Zmd;
AllL = Zlp+Zld;
FracPD = AllP ./ (AllP+AllD);
FracPF = AllP ./ (AllP+AllF);
FracLM = AllL ./ (AllM+AllL);

%% bent
% figure(1)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(Zb))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-2.5 0.5]);
% hcb = colorbar('h');
% ylim(hcb,[-2.5 0.5])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology log10 mean Benthic inverts (g m^-^2)')
% stamp(simname)
% print('-dpng',[pp 'global_BENT_',simname,'.png'])

%% All 4 on subplots
figure(2)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 mean All F (g m^-^2)')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All D (g m^-^2)')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All P (g m^-^2)')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All fishes (g m^-^2)')
stamp(simname)
print('-dpng',[pp 'global_All_subplot_',simname,'.png'])

%% 3 figure subplot P:D, P:F, M:L
% figure(3)
% subplot('Position',[0 0.53 0.5 0.5])
% %P:D
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,FracPD)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0 1]);
% set(gcf,'renderer','painters')
% title('Fraction Large Pelagics vs. Demersals')
% 
% %P:F
% subplot('Position',[0.5 0.53 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,FracPF)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0 1]);
% set(gcf,'renderer','painters')
% title('Fraction Large Pelagics vs. Forage Fishes')
% 
% %L:M
% subplot('Position',[0.25 0.0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,FracLM)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0 1]);
% colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
% set(gcf,'renderer','painters')
% title('Fraction Large vs. Medium')
% stamp(simname)
% print('-dpng',[pp 'global_ratios_subplot_',simname,'.png'])

end
