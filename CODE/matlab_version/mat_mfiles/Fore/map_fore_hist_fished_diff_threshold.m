% Visualize difference between
% ESM2M Hindcast of 1951-2000 and Forecast of 2051-2100
% Only in regions where historic biomass was greater than threshold

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
%fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
fpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

% saved hist and fore together
load([fpath 'Means_hist_fore_',harv,'_cobalt_' cfile '.mat'],'mzloss_hist',...
    'mzprod_hist','lzloss_hist','lzprod_hist','npp_hist','det_hist',...
    'ptemp_hist','mzloss_fore','mzprod_fore','lzloss_fore','lzprod_fore',...
    'npp_fore','det_fore','ptemp_fore','cF','cP','cD','cS','cM','cL','hF',...
    'hP','hD','hS','hM','hL','Hb','Cb');

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

%%
cAll = cF+cP+cD;
hAll = hF+hP+hD;
cPel = cF+cP;
hPel = hF+hP;

pdiffF = (cF-hF) ./ hF;
pdiffP = (cP-hP) ./ hP;
pdiffD = (cD-hD) ./ hD;
pdiffB = (Cb-Hb) ./ Hb;
pdiffAll = (cAll-hAll) ./ hAll;
pdiffPel = (cPel-hPel) ./ hPel;

hmaskF = nan*ones(size(hAll));
hmaskP = nan*ones(size(hAll));
hmaskD = nan*ones(size(hAll));
hmaskB = nan*ones(size(hAll));
hmaskA = nan*ones(size(hAll));
hmaskPel = nan*ones(size(hAll));
hmaskF(hF(:)>1e-6) = 1;
hmaskP(hP(:)>1e-6) = 1;
hmaskD(hD(:)>1e-6) = 1;
hmaskB(Hb(:)>1e-6) = 1;
hmaskA(hAll(:)>1e-6) = 1;
hmaskPel(hPel(:)>1e-6) = 1;

%% saved masks
save([fpath 'Means_hist_fore_',harv,'_biom_masks_' cfile '.mat'],...
    'hmaskF','hmaskP','hmaskD','hmaskB','hmaskA','hmaskPel');

%% Maps
% 6 plot pdiff bio subplot for ms Supp
figure(1)
%A - F
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffF.*hmaskF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
%colorbar('Position',[0.385 0.695 0.025 0.275],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Forage','HorizontalAlignment','center')
text(-2.75,1.75,'A')
%B - B
subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffB.*hmaskB)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
%colorbar('Position',[0.825 0.695 0.025 0.275],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Benthos','HorizontalAlignment','center')
text(-2.75,1.75,'B')
%C - P
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffP.*hmaskP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
%colorbar('Position',[0.385 0.385 0.025 0.275],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Large pelagic','HorizontalAlignment','center')
text(-2.75,1.75,'C')
%D - D
subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffD.*hmaskD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar('Position',[0.8 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Demersal','HorizontalAlignment','center')
text(-2.75,1.75,'D')
%E - Pel
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffPel.*hmaskPel)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
%colorbar('Position',[0.385 0.075 0.025 0.275],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'All pelagic','HorizontalAlignment','center')
text(-2.75,1.75,'E')
%F - All
subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffAll.*hmaskA)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
%colorbar('Position',[0.825 0.075 0.025 0.275],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'All fish','HorizontalAlignment','center')
text(-2.75,1.75,'F')
print('-dpng',[pp 'Hist_Fore_',harv,'_global_pdiff_biom_6plot_mask.png'])



