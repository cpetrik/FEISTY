% Visualize output of encounter and consumption
% Jared noticed consumption after "overcon" subroutine 
% was causing patchwork
% - I do not see this when only running encounter, consump, and offline_zm
% - must be a result of movement?
% - also, consump looks more like temp than prey, so change that behavior

Sf.con_before = Sf.con_zm .* Sf.bio;
Sp.con_before = Sp.con_zm .* Sp.bio;
Sd.con_before = Sd.con_zm .* Sd.bio;
Mf.con_before = Mf.con_zm .* Mf.bio;
Mp.con_before = Mp.con_zm .* Mp.bio;

Sf.con_after = Sf.con_zm .* Sf.bio;
Sp.con_after = Sp.con_zm .* Sp.bio;
Sd.con_after = Sd.con_zm .* Sd.bio;
Mf.con_after = Mf.con_zm .* Mf.bio;
Mp.con_after = Mp.con_zm .* Mp.bio;

%%
load('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');

[ni,nj]=size(geolon_t);
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; 
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];


%% Plots in space
TP=NaN*ones(ni,nj);
TP(GRD.ID)=COBALT.Tp(:,4);

Zm=NaN*ones(ni,nj);
Zm(GRD.ID)=COBALT.Zm(:,4);

dZm=NaN*ones(ni,nj);
dZm(GRD.ID)=COBALT.dZm(:,4);

% Encounter
Esf=NaN*ones(ni,nj);
Esp=NaN*ones(ni,nj);
Esd=NaN*ones(ni,nj);
Emf=NaN*ones(ni,nj);
Emp=NaN*ones(ni,nj);

Esf(GRD.ID)=Sf.enc_zm;
Esp(GRD.ID)=Sp.enc_zm;
Esd(GRD.ID)=Sd.enc_zm;
Emf(GRD.ID)=Mf.enc_zm;
Emp(GRD.ID)=Mp.enc_zm;

% 1st Consumption
Csf=NaN*ones(ni,nj);
Csp=NaN*ones(ni,nj);
Csd=NaN*ones(ni,nj);
Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);

Csf(GRD.ID)=Sf.con_before;
Csp(GRD.ID)=Sp.con_before;
Csd(GRD.ID)=Sd.con_before;
Cmf(GRD.ID)=Mf.con_before;
Cmp(GRD.ID)=Mp.con_before;

% 2nd Consumption (after offline_zm)
Osf=NaN*ones(ni,nj);
Osp=NaN*ones(ni,nj);
Osd=NaN*ones(ni,nj);
Omf=NaN*ones(ni,nj);
Omp=NaN*ones(ni,nj);

Osf(GRD.ID)=Sf.con_after;
Osp(GRD.ID)=Sp.con_after;
Osd(GRD.ID)=Sd.con_after;
Omf(GRD.ID)=Mf.con_after;
Omp(GRD.ID)=Mp.con_after;


%% encounter
figure(1)
% top left
subplot('Position',[0 0.66 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Esf)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
%clim([-2 2]);
set(gcf,'renderer','painters')
title('SF enc')

% mid left
subplot('Position',[0 0.33 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Esp)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
%clim([-2 2]);
set(gcf,'renderer','painters')
title('SP enc')

% btm left
subplot('Position',[0 0 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Esd)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
%clim([-2 2]);
set(gcf,'renderer','painters')
title('SD enc')

% top right
subplot('Position',[0.5 0.66 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Emf)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
%clim([-2 2]);
set(gcf,'renderer','painters')
title('MF enc')

% mid right
subplot('Position',[0.5 0.33 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Emp)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
%clim([-2 2]);
set(gcf,'renderer','painters')
title('MP enc')
%colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')

%btm right
subplot('Position',[0.5 0 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Zm)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
%clim([-2 2]);
set(gcf,'renderer','painters')
title('MZ Biomass')

stamp('')
%print('-dpng',[ppath mod 'global_BENT.png'])

%% consump before
figure(2)
% top left
subplot('Position',[0 0.66 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Csf)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
clim([0 0.2]);
set(gcf,'renderer','painters')
title('SF con')

% mid left
subplot('Position',[0 0.33 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Csp)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
clim([0 0.2]);
set(gcf,'renderer','painters')
title('SP con')

% btm left
subplot('Position',[0 0 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Csd)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
clim([0 0.2]);
set(gcf,'renderer','painters')
title('SD con')

% top right
subplot('Position',[0.5 0.66 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Cmf)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
clim([0 0.2]);
set(gcf,'renderer','painters')
title('MF con')

% mid right
subplot('Position',[0.5 0.33 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Cmp)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
clim([0 0.2]);
set(gcf,'renderer','painters')
title('MP con')
%colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')

%btm right
subplot('Position',[0.5 0 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,TP)
cmocean('thermal')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
%clim([-2 2]);
set(gcf,'renderer','painters')
title('Pel temp')

stamp('')
%print('-dpng',[ppath mod 'global_BENT.png'])

%% consumption after
figure(3)
% top left
subplot('Position',[0 0.66 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Osf)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
clim([0 0.2]);
set(gcf,'renderer','painters')
title('SF overcon')

% mid left
subplot('Position',[0 0.33 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Osp)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
clim([0 0.2]);
set(gcf,'renderer','painters')
title('SP overcon')

% btm left
subplot('Position',[0 0 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Osd)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
clim([0 0.2]);
set(gcf,'renderer','painters')
title('SD overcon')

% top right
subplot('Position',[0.5 0.66 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omf)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
clim([0 0.2]);
set(gcf,'renderer','painters')
title('MF overcon')

% mid right
subplot('Position',[0.5 0.33 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Omp)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
clim([0 0.2]);
set(gcf,'renderer','painters')
title('MP overcon')
%colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')

%btm right
subplot('Position',[0.5 0 0.45 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,dZm)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
clim([0 0.2]);
set(gcf,'renderer','painters')
title('MZ HP loss')

stamp('')
%print('-dpng',[ppath mod 'global_BENT.png'])
