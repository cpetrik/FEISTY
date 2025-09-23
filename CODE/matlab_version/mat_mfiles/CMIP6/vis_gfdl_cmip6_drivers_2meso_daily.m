% Visualize COBALT output 
% CMIP6 Historic global on 1 degree grid
% w/ 2 mesozoo
% Saved as mat files

clear 
close all

%%
gpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/GFDL/';

ppath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/FishMIP/FishMIP6/';

load([gpath 'hist/Time_Means_gfdl_cmip6_hist_2meso_daily.mat']);
load([gpath 'hist/Space_Means_gfdl_cmip6_hist_2meso_daily.mat']);

load([gpath 'Data_grid_gfdl.mat'],'GRD');
load([gpath 'gridspec_gfdl_cmip6.mat']);

%%
[ni,nj]=size(LON);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;   %decent looking coastlines

%% colors
cm10=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0.5 0.5 0.5; ...    %med grey
    0 0 0];...      %black
    

set(groot,'defaultAxesColorOrder',cm10);

%% Plots in space
TP=NaN*ones(ni,nj);
TB=NaN*ones(ni,nj);
Zm=NaN*ones(ni,nj);
Zl=NaN*ones(ni,nj);
dZm=NaN*ones(ni,nj);
dZl=NaN*ones(ni,nj);
Det=NaN*ones(ni,nj);

TP(GRD.ID)=tp_smean;
TB(GRD.ID)=tb_smean;
Zm(GRD.ID)=mz_smean;
Zl(GRD.ID)=lz_smean;
dZm(GRD.ID)=mzhp_smean;
dZl(GRD.ID)=lzhp_smean;
Det(GRD.ID)=det_smean;
	
%% 
figure(1)
subplot(3,3,1)
plot(yrs,tp_tmean,'r','Linewidth',1); hold on;
title('TP')

subplot(3,3,2)
plot(yrs,tb_tmean,'b','Linewidth',1); hold on;
title('TB')

subplot(3,3,3)
plot(yrs,(det_tmean),'k','Linewidth',1); hold on;
title('Det')

subplot(3,3,4)
plot(yrs,(mz_tmean),'m','Linewidth',1); hold on;
title('MZ')

subplot(3,3,5)
plot(yrs,(lz_tmean),'color',[0 0.5 0.75],'Linewidth',1); hold on;
title('LZ')

subplot(3,3,7)
plot(yrs,(mzhp_tmean),'m','Linewidth',1); hold on;
xlabel('Year')
title('MZ HPloss')

subplot(3,3,8)
plot(yrs,(lzhp_tmean),'color',[0 0.5 0.75],'Linewidth',1); hold on;
xlabel('Year')
title('LZ HPloss')

stamp('')
print('-dpng',[ppath 'gfdl_cmip6_hist_2meso_ts_all_drivers_annmean_daily.png'])


%% 8plot by fn type and size
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - TP
subplot('Position',[0.015 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,TP)
cmocean('thermal')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 30])
set(gcf,'renderer','painters')
text(0,1.75,'TP','HorizontalAlignment','center')

%B - MZ
subplot('Position',[0.015 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(Zm))
cmocean('dense')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2])
set(gcf,'renderer','painters')
text(0,1.75,'MZ','HorizontalAlignment','center')

%C - MZloss
subplot('Position',[0.015 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(dZm))
cmocean('dense')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3 0])
set(gcf,'renderer','painters')
text(0,1.75,'MZ HPloss','HorizontalAlignment','center')

%D - Det
subplot('Position',[0.015 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(Det))
cmocean('dense')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2.5 0.5])
set(gcf,'renderer','painters')
text(0,1.75,'Det','HorizontalAlignment','center')

%E - TB
subplot('Position',[0.47 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,TB)
cmocean('thermal')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 30])
set(gcf,'renderer','painters')
text(0,1.75,'TB','HorizontalAlignment','center')

%F - LZ
subplot('Position',[0.47 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(Zl))
cmocean('dense')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2])
set(gcf,'renderer','painters')
text(0,1.75,'LZ','HorizontalAlignment','center')

%G - LZloss
subplot('Position',[0.47 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(dZl))
cmocean('dense')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3 0])
set(gcf,'renderer','painters')
text(0,1.75,'LZ HPloss','HorizontalAlignment','center')

%H - all
% subplot('Position',[0.47 0.0 0.44 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(LAT,LON,log10(Eld))
% cmocean('dense')
% %colorbar
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-2 2])
% set(gcf,'renderer','painters')
% text(0,1.75,'LD','HorizontalAlignment','center')

stamp('')
print('-dpng',[ppath 'gfdl_cmip6_hist_2meso_global_means_daily.png'])
