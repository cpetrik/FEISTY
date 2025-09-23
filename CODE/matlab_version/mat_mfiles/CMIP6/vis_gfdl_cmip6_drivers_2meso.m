% Visualize COBALT output 
% CMIP6 Historic global on 1 degree grid
% w/ 2 mesozoo
% Saved as mat files

clear 
close all

%%
gpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/GFDL/';

ppath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/FishMIP/FishMIP6/';

load([gpath 'hist/Time_Means_gfdl_cmip6_hist_2meso.mat']);
load([gpath 'hist/Space_Means_gfdl_cmip6_hist_2meso.mat']);

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

%% time
%units = 'days since 1850-01-01 00:00:00'
y = (time/365) + 1850;

yr = floor(y(1)):floor(y(end));

%% annual means
nt = length(time);
nyr = nt/12;

st=1:12:nt;
en=12:12:nt;

[ni,nj] = size(LAT);

tp_amean = nan*ones(length(nyr),1);
tb_amean = tp_amean;
mz_amean = tp_amean;
lz_amean = tp_amean;
mzhp_amean = tp_amean;
lzhp_amean = tp_amean;
det_amean = tp_amean;

for n=1:length(st)
    % mean 
    tp_amean(n)=mean(tp_tmean(st(n):en(n)),'omitnan');
    tb_amean(n)=mean(tb_tmean(st(n):en(n)),'omitnan');
    mz_amean(n)=mean(mz_tmean(st(n):en(n)),'omitnan');
    lz_amean(n)=mean(lz_tmean(st(n):en(n)),'omitnan');
    mzhp_amean(n)=mean(mzhp_tmean(st(n):en(n)),'omitnan');
    lzhp_amean(n)=mean(lzhp_tmean(st(n):en(n)),'omitnan');
    det_amean(n)=mean(det_tmean(st(n):en(n)),'omitnan');

end


%% 
figure(1)
subplot(3,3,1)
plot(y(3:12:end),tp_tmean(3:12:end),'r','Linewidth',1); hold on;
title('TP')

subplot(3,3,2)
plot(y(3:12:end),tb_tmean(3:12:end),'b','Linewidth',1); hold on;
title('TB')

subplot(3,3,3)
plot(y(3:12:end),(det_tmean(3:12:end)),'k','Linewidth',1); hold on;
title('Det')

subplot(3,3,4)
plot(y(3:12:end),(mz_tmean(3:12:end)),'m','Linewidth',1); hold on;
title('MZ')

subplot(3,3,5)
plot(y(3:12:end),(lz_tmean(3:12:end)),'color',[0 0.5 0.75],'Linewidth',1); hold on;
title('LZ')

subplot(3,3,7)
plot(y(3:12:end),(mzhp_tmean(3:12:end)),'m','Linewidth',1); hold on;
xlabel('Time (mo)')
title('MZ HPloss')

subplot(3,3,8)
plot(y(3:12:end),(lzhp_tmean(3:12:end)),'color',[0 0.5 0.75],'Linewidth',1); hold on;
xlabel('Time (mo)')
title('LZ HPloss')

stamp('')
print('-dpng',[ppath 'gfdl_cmip6_hist_2meso_ts_all_drivers_March.png'])
	
	
%% 
figure(3)
subplot(3,3,1)
plot(yr,tp_amean,'r','Linewidth',1); hold on;
title('TP')

subplot(3,3,2)
plot(yr,tb_amean,'b','Linewidth',1); hold on;
title('TB')

subplot(3,3,3)
plot(yr,(det_amean),'k','Linewidth',1); hold on;
title('Det')

subplot(3,3,4)
plot(yr,(mz_amean),'m','Linewidth',1); hold on;
title('MZ')

subplot(3,3,5)
plot(yr,(lz_amean),'color',[0 0.5 0.75],'Linewidth',1); hold on;
title('LZ')

subplot(3,3,7)
plot(yr,(mzhp_amean),'m','Linewidth',1); hold on;
xlabel('Year')
title('MZ HPloss')

subplot(3,3,8)
plot(yr,(lzhp_amean),'color',[0 0.5 0.75],'Linewidth',1); hold on;
xlabel('Year')
title('LZ HPloss')

stamp('')
print('-dpng',[ppath 'gfdl_cmip6_hist_2meso_ts_all_drivers_annmean.png'])


%% 8plot by fn type and size
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - TP
subplot('Position',[0.015 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,tp_smean)
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
surfm(LAT,LON,log10(mz_smean))
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
surfm(LAT,LON,log10(mzhp_smean))
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
surfm(LAT,LON,log10(det_smean))
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
surfm(LAT,LON,tb_smean)
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
surfm(LAT,LON,log10(lz_smean))
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
surfm(LAT,LON,log10(lzhp_smean))
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
print('-dpng',[ppath 'gfdl_cmip6_hist_2meso_global_means.png'])
