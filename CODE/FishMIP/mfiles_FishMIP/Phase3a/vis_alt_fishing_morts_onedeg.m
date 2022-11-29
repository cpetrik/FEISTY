% Visualize fishing mort input of FEISTY
% Time series plots and maps

clear 
close all

%% load f/fmsy

%% v1 1841-1960
%/Users/cpetrik/Dropbox/Princeton/FEISTY_other/fishing_ms_ideas/fishing_effort_impl/grid_mortality_guilds_v1
alt1 = 'grid_mortality_guilds_v1'; 
alt2 = '_v1';

gpath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY_other/fishing_ms_ideas/fishing_effort_impl/',alt1,'/'];

% load([gpath 'grid_mortality_all',alt2,'.mat'])

load([gpath 'gfdl-mom6-cobalt2_ctrlclim_onedeg_fmort_ID_annual_1841_1960_tempSc.mat'],'year','WID',...
    'fmD','fmF','fmP');

FmD_v1 = fmD;
FmF_v1 = fmF;
FmP_v1 = fmP;

clear fmD fmF fmP

%% 1961-2010
% load([gpath 'gfdl-mom6-cobalt2_ctrlclim_onedeg_fmort_ID_annual_1961_2010_tempSc.mat'],'year','WID',...
%     'fmD','fmF','fmP');

load([gpath 'gfdl-mom6-cobalt2_obsclim_onedeg_fmort_ID_annual_1961_2010_tempSc.mat'],'year','WID',...
    'fmD','fmF','fmP');

FmD_v1(:,121:170) = fmD;
FmF_v1(:,121:170) = fmF;
FmP_v1(:,121:170) = fmP;

clear fmD fmF fmP

%% v1.2 1841-1960
alt1 = 'pristine_grid_mortality_guilds_v1'; 
alt2 = '_v1_pristine'; 

hpath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY_other/fishing_ms_ideas/fishing_effort_impl/',alt1,'/'];

% load([hpath 'grid_mortality_all',alt2,'.mat'])

load([hpath 'gfdl-mom6-cobalt2_ctrlclim_onedeg_fmort_ID_annual_1841_1960_tempSc',alt2,'.mat'],'year','WID',...
    'fmD','fmF','fmP');

FmD_v12 = fmD;
FmF_v12 = fmF;
FmP_v12 = fmP;

clear fmD fmF fmP

%% 1961-2010
% load([hpath 'gfdl-mom6-cobalt2_ctrlclim_onedeg_fmort_ID_annual_1961_2010_tempSc',alt2,'.mat'],'year','WID',...
%     'fmD','fmF','fmP');

load([hpath 'gfdl-mom6-cobalt2_obsclim_onedeg_fmort_ID_annual_1961_2010_tempSc',alt2,'.mat'],'year','WID',...
    'fmD','fmF','fmP');

FmD_v12(:,121:170) = fmD;
FmF_v12(:,121:170) = fmF;
FmP_v12(:,121:170) = fmP;

clear fmD fmF fmP

%% v2 1841-1960
alt1 = 'grid_mortality_guilds_v2';
alt2 = '_v2'; 

ipath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY_other/fishing_ms_ideas/fishing_effort_impl/',alt1,'/'];

% load([ipath 'grid_mortality_all',alt2,'.mat'])

load([ipath 'gfdl-mom6-cobalt2_ctrlclim_onedeg_fmort_ID_annual_1841_1960_tempSc',alt2,'.mat'],'year','WID',...
    'fmD','fmF','fmP');

FmD_v2 = fmD;
FmF_v2 = fmF;
FmP_v2 = fmP;

clear fmD fmF fmP

%% 1961-2010
% load([ipath 'gfdl-mom6-cobalt2_ctrlclim_onedeg_fmort_ID_annual_1961_2010_tempSc',alt2,'.mat'],'year','WID',...
%     'fmD','fmF','fmP');

load([ipath 'gfdl-mom6-cobalt2_obsclim_onedeg_fmort_ID_annual_1961_2010_tempSc',alt2,'.mat'],'year','WID',...
    'fmD','fmF','fmP');

FmD_v2(:,121:170) = fmD;
FmF_v2(:,121:170) = fmF;
FmP_v2(:,121:170) = fmP;

clear fmD fmF fmP

%%
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP/Phase3a/';
ppath = [pp 'OneDeg/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

% Map data
%cpath = '/Volumes/MIP/Fish-MIP/Phase3/OneDeg/';
cpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';
load([cpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat']);
load([cpath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat'], 'GRD');

[ni,nj]=size(LON);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

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

load coastlines;                     %decent looking coastlines

%% Time means
tfd1 = mean(FmD_v1); %'omitnan'
tff1 = mean(FmF_v1);
tfp1 = mean(FmP_v1);

tfd12 = mean(FmD_v12); %'omitnan'
tff12 = mean(FmF_v12);
tfp12 = mean(FmP_v12);

tfd2 = mean(FmD_v2); %'omitnan'
tff2 = mean(FmF_v2);
tfp2 = mean(FmP_v2);

%% Space means
sfd1 = mean(FmD_v1,2); %'omitnan'
sff1 = mean(FmF_v1,2);
sfp1 = mean(FmP_v1,2);

sfd12 = mean(FmD_v12,2); %'omitnan'
sff12 = mean(FmF_v12,2);
sfp12 = mean(FmP_v12,2);

sfd2 = mean(FmD_v2,2); %'omitnan'
sff2 = mean(FmF_v2,2);
sfp2 = mean(FmP_v2,2);

%% Plots in time
%t = 1:length(sp_tmean); %time;
y = 1841:2010;

figure(1)
plot(y,log10(tff1),'Linewidth',1); hold on;
plot(y,log10(tff12),'Linewidth',1); hold on;
plot(y,log10(tfp1),'Linewidth',1); hold on;
plot(y,log10(tfp12),'Linewidth',1); hold on;
plot(y,log10(tfp2),'Linewidth',1); hold on;
plot(y,log10(tfd1),'Linewidth',1); hold on;
plot(y,log10(tfd12),'Linewidth',1); hold on;
plot(y,log10(tfd2),'Linewidth',1); hold on;
plot(y,log10(tff2),'Linewidth',1); hold on;
legend('F1','F12','P1','P12','P2','D1','D12','D2','F2')
legend('location','eastoutside')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('log_1_0 fmort (g m^-^2 y^-^1)')
%title('Hist')
stamp('')
print('-dpng',[ppath 'onedeg_fmort_annual_tempSc_alts_time_mean_log10.png'])

%%
figure(2)
plot(y,(tff1),'Linewidth',1); hold on;
plot(y,(tff12),'Linewidth',1); hold on;
plot(y,(tfp1),'Linewidth',1); hold on;
plot(y,(tfp12),'Linewidth',1); hold on;
plot(y,(tfp2),'Linewidth',1); hold on;
plot(y,(tfd1),'Linewidth',1); hold on;
plot(y,(tfd12),'Linewidth',1); hold on;
plot(y,(tfd2),'Linewidth',1); hold on;
plot(y,(tff2),'Linewidth',1); hold on;
legend('F1','F12','P1','P12','P2','D1','D12','D2','F2')
legend('location','eastoutside')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('fmort (g m^-^2 y^-^1)')
%title('Hist')
stamp('')
print('-dpng',[ppath 'onedeg_fmort_annual_tempSc_alts_time_mean.png'])

 
%% Plots in space
Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);

Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);

Zlf=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);

Zsf(GRD.ID)=sff1;
Zsp(GRD.ID)=sfp1;
Zsd(GRD.ID)=sfd1;

Zmf(GRD.ID)=sff12;
Zmp(GRD.ID)=sfp12;
Zmd(GRD.ID)=sfd12;

Zlf(GRD.ID)=sff2;
Zlp(GRD.ID)=sfp2;
Zld(GRD.ID)=sfd2;

%% All 3 versions of F on subplots
figure(4)
% v1
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(Zsf))
cmocean('matter')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.5]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('v1     log10 mean fmort F (g m^-^2 y^-^1)')

% v1.2
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(Zmf))
cmocean('matter')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.5]);
set(gcf,'renderer','painters')
title('v1.2')

% v2
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(Zlf))
cmocean('matter')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.5]);
set(gcf,'renderer','painters')
title('v2')

% All
% subplot('Position',[0.5 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,log10(All))
% cmocean('matter')
% h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-2 2]);
% set(gcf,'renderer','painters')
% title('log10 mean All fishes (g m^-^2)')

stamp('')
print('-dpng',[ppath 'onedeg_fmort_annual_tempSc_alts_map_mean_F.png'])


%% All 3 versions of P on subplots
figure(5)
% v1
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(Zsp))
cmocean('matter')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.5]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('v1     log10 mean fmort P (g m^-^2 y^-^1)')

% v1.2
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(Zmp))
cmocean('matter')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.5]);
set(gcf,'renderer','painters')
title('v1.2')

% v2
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(Zlp))
cmocean('matter')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.5]);
set(gcf,'renderer','painters')
title('v2')

stamp('')
print('-dpng',[ppath 'onedeg_fmort_annual_tempSc_alts_map_mean_P.png'])


%% All 3 versions of D on subplots
figure(6)
% v1
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(Zsd))
cmocean('matter')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.5]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('v1     log10 mean fmort D (g m^-^2 y^-^1)')

% v1.2
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(Zmd))
cmocean('matter')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.5]);
set(gcf,'renderer','painters')
title('v1.2')

% v2
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(Zld))
cmocean('matter')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.5]);
set(gcf,'renderer','painters')
title('v2')

stamp('')
print('-dpng',[ppath 'onedeg_fmort_annual_tempSc_alts_map_mean_D.png'])

%% Diffs red-white-blue
figure(7)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(Zmf - Zsf))
cmocean('balance')
caxis([-0.3 0.3])
title('F v1.2 - v1')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(Zlf - Zsf))
cmocean('balance')
caxis([-0.3 0.3])
title('F v2 - v1')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(Zmp - Zsp))
cmocean('balance')
caxis([-0.3 0.3])
title('P v1.2 - v1')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(Zlp - Zsp))
cmocean('balance')
cb = colorbar('Position',[0.85 0.4 0.03 0.4],'orientation','vertical');
xlabel(cb,'fmort')
caxis([-0.3 0.3])
title('P v2 - v1')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(Zmd - Zsd))
cmocean('balance')
caxis([-0.3 0.3])
title('D v1.2 - v1')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(Zld - Zsd))
cmocean('balance')
caxis([-0.3 0.3])
title('D v2 - v1')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);

stamp('')
print('-dpng',[ppath 'onedeg_fmort_annual_tempSc_alts_map_mean_diffs.png'])

%% Hist period
figure(11)
plot(y,log10(tff1),'Linewidth',1); hold on;
plot(y,log10(tff12),'Linewidth',1); hold on;
plot(y,log10(tfp1),'Linewidth',1); hold on;
plot(y,log10(tfp12),'Linewidth',1); hold on;
plot(y,log10(tfp2),'Linewidth',1); hold on;
plot(y,log10(tfd1),'Linewidth',1); hold on;
plot(y,log10(tfd12),'Linewidth',1); hold on;
plot(y,log10(tfd2),'Linewidth',1); hold on;
plot(y,log10(tff2),'Linewidth',1); hold on;
legend('F1','F12','P1','P12','P2','D1','D12','D2','F2')
legend('location','eastoutside')
xlim([1961 2010])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('log_1_0 fmort (g m^-^2 y^-^1)')
%title('Hist')
stamp('')
print('-dpng',[ppath 'Hist_onedeg_fmort_annual_tempSc_alts_time_mean_log10.png'])

%%
figure(12)
plot(y,(tff1),'Linewidth',1); hold on;
plot(y,(tff12),'Linewidth',1); hold on;
plot(y,(tfp1),'Linewidth',1); hold on;
plot(y,(tfp12),'Linewidth',1); hold on;
plot(y,(tfp2),'Linewidth',1); hold on;
plot(y,(tfd1),'Linewidth',1); hold on;
plot(y,(tfd12),'Linewidth',1); hold on;
plot(y,(tfd2),'Linewidth',1); hold on;
plot(y,(tff2),'Linewidth',1); hold on;
legend('F1','F12','P1','P12','P2','D1','D12','D2','F2')
legend('location','eastoutside')
xlim([1961 2010])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('fmort (g m^-^2 y^-^1)')
%title('Hist')
stamp('')
print('-dpng',[ppath 'Hist_onedeg_fmort_annual_tempSc_alts_time_mean.png'])


