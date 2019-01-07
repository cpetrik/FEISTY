% Estimate POEM harvest with different efforts scaled from land distance
% Climatology
% 150 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/CSV/Matlab_new_size/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'esm26_1deg_shore_dist.mat']);

AREA_OCN = max(area,1);
tlme = lme_mask_onedeg;
mask = tlme;
mask(isnan(mask))=0;
mask(mask>0)=1;

%% Plots in space
% plot info
[ni,nj]=size(lon);
geolat_t=lat;
geolon_t=lon;
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));



%% Effort scalings
mins = [0;0;0;0.5;0.5;0.5];
maxs = [1;1.5;2;1;1.5;2];

min_dist1 = -1*min_dist;
maxd = max(min_dist1(:));
mind = min(min_dist1(:));

%% apply effort scaling
for n=1%:length(mins)
    eff1 = mins(n) + (maxs(n)-mins(n)).*(min_dist1-mind) ./ (maxd-mind);
end

%%
min_dist2 = -1*min_dist.^2;
maxd = max(min_dist2(:));
mind = min(min_dist2(:));
eff2 = mins(n) + (maxs(n)-mins(n)).*(min_dist2-mind) ./ (maxd-mind);

min_dist3 = -1*min_dist.^(1/2);
maxd = max(min_dist3(:));
mind = min(min_dist3(:));
eff3 = mins(n) + (maxs(n)-mins(n)).*(min_dist3-mind) ./ (maxd-mind);

min_dist4 = -1*min_dist.^(1/3);
maxd = max(min_dist4(:));
mind = min(min_dist4(:));
eff4 = mins(n) + (maxs(n)-mins(n)).*(min_dist4-mind) ./ (maxd-mind);

%The logistic function is given by:
% f(x) = L / (1 + exp(-k*(x-x0)));
%where
%L ? Curve?s maximum value
%k ? Steepness of the curve
%x0 ? x value of Sigmoid?s midpoint
min_dist5 = 1 ./ (1 + exp(-(10/500)*(min_dist-500)));
maxd = max(min_dist5(:));
mind = min(min_dist5(:));
eff5 = mins(n) + (maxs(n)-mins(n)).*(min_dist5-mind) ./ (maxd-mind);

min_dist6 = 1+tanh((10/500)*(min_dist-500));
maxd = max(min_dist6(:));
mind = min(min_dist6(:));
eff6 = mins(n) + (maxs(n)-mins(n)).*(min_dist6-mind) ./ (maxd-mind);

min_dist7 = exp(-500.*exp((-10/500)*min_dist));
maxd = max(min_dist7(:));
mind = min(min_dist7(:));
eff7 = mins(n) + (maxs(n)-mins(n)).*(min_dist7-mind) ./ (maxd-mind);

%%
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,eff1.*mask)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('linear')
print('-dpng',[pp 'Effort_dist_scaling_line.png'])

figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,eff2.*mask)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('quad')
print('-dpng',[pp 'Effort_dist_scaling_quad.png'])

figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,eff3.*mask)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('sqrt')
print('-dpng',[pp 'Effort_dist_scaling_sqrt.png'])

figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,eff4.*mask)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('cube root')
print('-dpng',[pp 'Effort_dist_scaling_cube.png'])

figure(5)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(1-eff5).*mask)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('sigmoid')
print('-dpng',[pp 'Effort_dist_scaling_sig.png'])

figure(6)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(1-eff6).*mask)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('tanh')
print('-dpng',[pp 'Effort_dist_scaling_tanh.png'])

figure(7)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(1-eff7).*mask)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Gompertz')
print('-dpng',[pp 'Effort_dist_scaling_gompertz.png'])


%%
figure(10)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,eff1)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('linear')

figure(11)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,eff2)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('quad')

figure(12)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,eff3)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('sqrt')

figure(13)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,eff4)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('cube root')

figure(14)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,eff5)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('sigmoid')

figure(15)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,1-eff6)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('tanh')

figure(16)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,1-eff7)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Gompertz')

%%
figure(17)
plot(min_dist,eff1,'.k'); hold on;
plot(min_dist,eff2,'.b'); hold on;
plot(min_dist,eff3,'.r'); hold on;
plot(min_dist,eff4,'.g'); hold on;
plot(min_dist,1-eff5,'.c'); hold on;
plot(min_dist,1-eff6,'.m'); hold on;
plot(min_dist,1-eff7,'.','color',[.5 .5 .5]); hold on;
xlim([0 500])

