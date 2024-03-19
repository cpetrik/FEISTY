% Plot mean biomasses of zoop, det, temp

clear 
close all

%%
fpath='/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/HADdown/';

% Figures
ppath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CC/';

%% 
load([fpath 'feisty_hadley_gridspec.mat'],'TMO','LAT','LON')
load([fpath 'Data_grid_nemuro_hadley.mat'],'GRD');

load([fpath 'nemuro_hadley_means.mat']);

%%
time = (TMO/365)+1900;
yrs= 1980:2100;
thist = find(yrs>=1980 & yrs<2001);
tssp = find(yrs>=2080 & yrs<2101);

%% empirical loss calc
D_dZm = 10 .^ (-2.617 + 1.989.*log10(D_Zm+eps) + 1.732e-2.*D_Tp);
D_dZl = 10 .^ (-2.954 + 2.228.*log10(D_Zl+eps) + 2.556e-2.*D_Tp);

dmz_tmean_hadley = double(mean(D_dZm,1,"omitnan"));
dlz_tmean_hadley = double(mean(D_dZl,1,"omitnan"));

dmz_hist_mean_hadley = double(mean(D_dZm,2,"omitnan"));
dlz_hist_mean_hadley = double(mean(D_dZl,2,"omitnan"));

dmz_ssp_mean_hadley = double(mean(D_dZm,2,"omitnan"));
dlz_ssp_mean_hadley = double(mean(D_dZl,2,"omitnan"));

%%
[ni,nj]=size(LON);

plotminlat=30; %Set these bounds for your data
plotmaxlat=48;
plotminlon=-134;
plotmaxlon=-115;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%% Plots in time
y = yrs;

% All size classes of all
figure(10)
subplot(3,1,1)
plot(y,ptemp_tmean_hadley,'r','Linewidth',1); hold on;
plot(y,btemp_tmean_hadley,'b','Linewidth',1); hold on;
legend('TP','TB')
xlim([y(1) y(end)])
ylabel('Temperature (^oC)')
title('HAD')

subplot(3,1,2)
plot(y,log10(mz_tmean_hadley),'color',[0 0.5 0.75],'Linewidth',1); hold on;
plot(y,log10(lz_tmean_hadley),'color',[0 0 0.75],'Linewidth',1); hold on;
plot(y,log10(det_tmean_hadley),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
legend('MZ','LZ','Det')
xlim([y(1) y(end)])
ylabel('log_1_0 Biomass (g m^-^2)')

subplot(3,1,3)
plot(y,(dmz_tmean_hadley),'color',[0 0.5 0.75],'Linewidth',1); hold on;
plot(y,(dlz_tmean_hadley),'color',[0 0 0.75],'Linewidth',1); hold on;
legend('dMZ','dLZ')
xlim([y(1) y(end)])
%ylim([-3 2])
xlabel('Time (mo)')
ylabel('Loss rate (g m^-^2 d^-^1)')
stamp('')
print('-dpng',[ppath 'HAD_mean_forcing_ts.png'])

%% Plots in space
H_MZ = NaN*ones(ni,nj);
H_LZ = NaN*ones(ni,nj);
H_MZloss=NaN*ones(ni,nj);
H_LZloss=NaN*ones(ni,nj);
H_Tp = NaN*ones(ni,nj);
H_Tb = NaN*ones(ni,nj);
H_Det = NaN*ones(ni,nj);

S_MZ = NaN*ones(ni,nj);
S_LZ = NaN*ones(ni,nj);
S_MZloss=NaN*ones(ni,nj);
S_LZloss=NaN*ones(ni,nj);
S_Tp = NaN*ones(ni,nj);
S_Tb = NaN*ones(ni,nj);
S_Det = NaN*ones(ni,nj);

H_MZ(GRD.ID) = mz_hist_mean_hadley;
H_LZ(GRD.ID) = lz_hist_mean_hadley;
H_MZloss(GRD.ID)=dmz_hist_mean_hadley;
H_LZloss(GRD.ID)=dlz_hist_mean_hadley;
H_Tp(GRD.ID) = ptemp_hist_mean_hadley;
H_Tb(GRD.ID) = btemp_hist_mean_hadley;
H_Det(GRD.ID) = det_hist_mean_hadley;

S_MZ(GRD.ID) = mz_ssp_mean_hadley;
S_LZ(GRD.ID) = lz_ssp_mean_hadley;
S_MZloss(GRD.ID)=dmz_ssp_mean_hadley;
S_LZloss(GRD.ID)=dlz_ssp_mean_hadley;
S_Tp(GRD.ID) = ptemp_ssp_mean_hadley;
S_Tb(GRD.ID) = btemp_ssp_mean_hadley;
S_Det(GRD.ID) = det_ssp_mean_hadley;

%% Temp
figure(1)
% Hist TP
subplot('Position',[0 0.51 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(H_Tp))
cmocean('thermal')
clim([5 20]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean Pelagic Temp (^oC) 1980-2000','HorizontalAlignment','center')

% Hist TB
subplot('Position',[0.5 0.51 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(H_Tb))
cmocean('thermal')
clim([0 15]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean Bottom Temp (^oC) 1980-2000','HorizontalAlignment','center')

% SSP TP
subplot('Position',[0 0 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(S_Tp))
cmocean('thermal')
clim([5 20]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean Pelagic Temp (^oC) 2080-2100','HorizontalAlignment','center')

% SSP TB
subplot('Position',[0.5 0 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(S_Tb))
cmocean('thermal')
clim([0 15]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean Bottom Temp (^oC) 2080-2100','HorizontalAlignment','center')
print('-dpng',[ppath 'temp_nemuro_hadley_hist_spp.png'])

%% Hist zoop
figure(2)
% Hist TP
subplot('Position',[0 0.51 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(H_MZ))
cmocean('tempo')
clim([0 20]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean MZ (g m^-^2) 1980-2000','HorizontalAlignment','center')

% Hist TB
subplot('Position',[0.5 0.51 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(H_LZ))
cmocean('tempo')
clim([0 15]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean LZ (g m^-^2) 1980-2000','HorizontalAlignment','center')

% SSP TP
subplot('Position',[0 0 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(H_MZloss))
cmocean('tempo')
clim([0 1]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean MZ loss (g m^-^2 d^-^1) 1980-2000','HorizontalAlignment','center')

% SSP TB
subplot('Position',[0.5 0 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(H_LZloss))
cmocean('tempo')
clim([0 0.5]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean LZ loss (g m^-^2 d^-^1) 1980-2000','HorizontalAlignment','center')
print('-dpng',[ppath 'zoop_nemuro_hadley_hist.png'])

%% SSP zoop
figure(3)
% Hist TP
subplot('Position',[0 0.51 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(S_MZ))
cmocean('tempo')
clim([0 20]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean MZ (g m^-^2) 2080-2100','HorizontalAlignment','center')

% Hist TB
subplot('Position',[0.5 0.51 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(S_LZ))
cmocean('tempo')
clim([0 15]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean LZ (g m^-^2) 2080-2100','HorizontalAlignment','center')

% SSP TP
subplot('Position',[0 0 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(S_MZloss))
cmocean('tempo')
clim([0 1]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean MZ loss (g m^-^2 d^-^1) 2080-2100','HorizontalAlignment','center')

% SSP TB
subplot('Position',[0.5 0 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(S_LZloss))
cmocean('tempo')
clim([0 0.5]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean LZ loss (g m^-^2 d^-^1) 2080-2100','HorizontalAlignment','center')
print('-dpng',[ppath 'zoop_nemuro_hadley_ssp.png'])

%% Det
figure(4)
% Hist TP
subplot('Position',[0 0.5 0.5 0.5])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(H_Det))
cmocean('tempo')
clim([0 3]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean Det (g m^-^2) 1980-2000','HorizontalAlignment','center')

% Hist TB
subplot('Position',[0.5 0.5 0.5 0.5])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(S_Det))
cmocean('tempo')
clim([0 3]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean Det (g m^-^2) 2080-2100','HorizontalAlignment','center')


subplot('Position',[0 0 0.5 0.5])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(H_Det))
cmocean('tempo')
clim([-2 2]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean Det (log_1_0 g m^-^2) 1980-2000','HorizontalAlignment','center')

% Hist TB
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(S_Det))
cmocean('tempo')
clim([-2 2]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean Det (log_1_0 g m^-^2) 2080-2100','HorizontalAlignment','center')
print('-dpng',[ppath 'det_nemuro_hadley_hist_spp.png'])

%% Hist zoop log10
figure(6)
% Hist TP
subplot('Position',[0 0.51 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(H_MZ))
cmocean('tempo')
clim([0.5 1.5]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean MZ (log_1_0 g m^-^2) 1980-2000','HorizontalAlignment','center')

% Hist TB
subplot('Position',[0.5 0.51 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(H_LZ))
cmocean('tempo')
clim([0.5 1.5]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean LZ (log_1_0 g m^-^2) 1980-2000','HorizontalAlignment','center')

% SSP TP
subplot('Position',[0 0 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(H_MZloss))
cmocean('tempo')
clim([-2 1]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean MZ loss (log_1_0 g m^-^2 d^-^1) 1980-2000','HorizontalAlignment','center')

% SSP TB
subplot('Position',[0.5 0 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(H_LZloss))
cmocean('tempo')
clim([-2 1]);
colorbar('h');
set(gcf,'renderer','painters')
text(0, 0.95,'mean LZ loss (log_1_0 g m^-^2 d^-^1) 1980-2000','HorizontalAlignment','center')
print('-dpng',[ppath 'zoop_log10_nemuro_hadley_hist.png'])
