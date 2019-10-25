% Plot mean biomasses of zoop to compare to POEM results
% COMPARE BY LME, COLOR-CODE BY TEMP, Plot GAM fit
% Uses 3 Knots and Probit link fn

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

gpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
cpath = [Pdrpbx 'Princeton/POEM_other/cobalt_data/'];
pp = [Pdrpbx 'Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([gpath 'esm26_area_1deg.mat']);
load([gpath 'LME_clim_temp_zoop_det_npp.mat']);
load([gpath 'LME_depth_area.mat'],'lme_depth','lme_shal_frac');

%% POEM results
frate = 0.3;
tfish = num2str(100+int64(10*frate));
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = ['All_fish',tfish(2:end)];
tharv = ['Harvest all fish ' num2str(frate)];
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'LME_ZBratios_clim_fished_',harv,'_' cfile '.mat']);
load([fpath 'All_gam_fits_DvD.mat']);

%% Assign a color to each LME based on temp
%tmap=colormap(jet(66));
tmap=cmocean('thermal',66);
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);
close all

%% Scatter plot with temp color
%
figure(1)
subplot('Position',[0.1 0.5 0.38 0.4])
plot(ZlDet(:,1),ZlDet(:,2),'k','LineWidth',1.5); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
scatter(log10(RatZlDet),FracPD,30,lme_ptemp(:,1),'filled'); hold on;
cmocean('thermal');
%colorbar
xlabel('log_1_0 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])
text(-0.95, 0.95,'A','FontSize',14)

%npp
subplot('Position',[0.52 0.5 0.38 0.4])
plot(npp(:,1),npp(:,2),'k','LineWidth',1.5); hold on;
plot(npp(:,1),npp(:,2)+2*npp(:,3),'--k'); hold on;
plot(npp(:,1),npp(:,2)-2*npp(:,3),'--k'); hold on;
scatter(log10(lme_npp/365),FracPD,30,lme_ptemp(:,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.925 0.5 0.025 0.395])
%title('Color = LME mean temp')
%ylabel('P / (P+D)')
xlabel('log_1_0 NPP')
axis([-0.5 1 0 1])
text(-0.45, 0.95, 'B','FontSize',14)

% set(gca,'FontSize',30,'XTick',[-1,-0.5,0,0.5],'XTickLabel',...
%     {'-1.0','-0.5','0.0','0.5'},...
%     'YTick',0:0.5:1,'YTickLabel',0:0.5:1)
print('-dpng',[ppath 'lme_scatter_GAMfit_PD_ZlDet_NPP.png'])


