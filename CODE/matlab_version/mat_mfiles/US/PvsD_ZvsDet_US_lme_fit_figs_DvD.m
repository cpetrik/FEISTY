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
ppath = [pp cfile '/US/'];
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

% US LMEs only
keep = [1:7,10,54,55,65];

%% Scatter plot black
%
figure(1)
subplot(3,3,1)
scatter(log10(RatZlDet(keep)),FracPD(keep),10,'k','filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
axis([-1 0.75 0 1])
ylabel('P / (P+D)')

subplot(3,3,4)
scatter(log10(RatZlDet(keep)),FracPF(keep),10,'k','filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,4),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)+2*ZlDet(:,5),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)-2*ZlDet(:,5),'--k'); hold on;
axis([-1 0.75 0 1])
ylabel('P / (P+F)')

subplot(3,3,7)
scatter(log10(RatZlDet(keep)),FracLM(keep),10,'k','filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,6),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)+2*ZlDet(:,7),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)-2*ZlDet(:,7),'--k'); hold on;
axis([-1 0.75 0 1])
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
print('-dpng',[ppath 'US_lme_scatter_ZlDet_GAMfit_black.png'])

%% Scatter plot with temp color
%
figure(2)
subplot(3,3,1)
scatter(log10(RatZlDet(keep)),FracPD(keep),10,lme_ptemp((keep),1),'filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
%title('Color = LME mean temp')
cmocean('thermal');
%colorbar
% xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(3,3,4)
scatter(log10(RatZlDet(keep)),FracPF(keep),10,lme_ptemp((keep),1),'filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,4),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)+2*ZlDet(:,5),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)-2*ZlDet(:,5),'--k'); hold on;
cmocean('thermal');
colorbar('Position',[0.375 0.4 0.02 0.25])
% xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
axis([-1 0.75 0 1])

subplot(3,3,7)
scatter(log10(RatZlDet(keep)),FracLM(keep),10,lme_ptemp((keep),1),'filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,6),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)+2*ZlDet(:,7),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)-2*ZlDet(:,7),'--k'); hold on;
cmocean('thermal');
%colorbar
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
axis([-1 0.75 0 1])
print('-dpng',[ppath 'US_lme_scatter_ZlDet_GAMfit_colorT_colorbar.png'])


%% Just P:D
figure(3)
plot(ZlDet(:,1),ZlDet(:,2),'LineWidth',2,'color',[0 0.5 0.75]); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
plot(log10(RatZlDet(keep)),FracPD(keep),'k.','MarkerSize',25); hold on;
axis([-1 0.75 0 1])
set(gca,'FontSize',18,'XTick',-0.9:0.3:0.6,'XTickLabel',-0.9:0.3:0.6,...
    'YTick',0:0.5:1,'YTickLabel',0:0.5:1)
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
print('-dpng',[ppath 'US_lme_scatter_ZlDet_GAMfit_PD.png'])


%% WITH OTHER GAMS
figure(4)
subplot(3,3,1)
scatter(log10(RatZlDet(keep)),FracPD(keep),10,'k','filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
%title('Color = LME mean temp')
%xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(3,3,4)
scatter(log10(RatZlDet(keep)),FracPF(keep),10,'k','filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,4),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)+2*ZlDet(:,5),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)-2*ZlDet(:,5),'--k'); hold on;
%xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
axis([-1 0.75 0 1])

subplot(3,3,7)
scatter(log10(RatZlDet(keep)),FracLM(keep),10,'k','filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,6),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)+2*ZlDet(:,7),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)-2*ZlDet(:,7),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
axis([-1 0.75 0 1])

%ptemp
subplot(3,3,2)
scatter(lme_ptemp((keep),1),FracPD(keep),10,'k','filled'); hold on;
plot(ptemp(:,1),ptemp(:,2),'k'); hold on;
plot(ptemp(:,1),ptemp(:,2)+2*ptemp(:,3),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,2)-2*ptemp(:,3),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+D)')
axis([-2 30 0 1])

subplot(3,3,5)
scatter(lme_ptemp((keep),1),FracPF(keep),10,'k','filled'); hold on;
plot(ptemp(:,1),ptemp(:,4),'k'); hold on;
plot(ptemp(:,1),ptemp(:,4)+2*ptemp(:,5),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,4)-2*ptemp(:,5),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+F)')
axis([-2 30 0 1])

subplot(3,3,8)
scatter(lme_ptemp((keep),1),FracLM(keep),10,'k','filled'); hold on;
plot(ptemp(:,1),ptemp(:,6),'k'); hold on;
plot(ptemp(:,1),ptemp(:,6)+2*ptemp(:,7),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,6)-2*ptemp(:,7),'--k'); hold on;
xlabel('Tpel (^oC)')
%ylabel('L / (L+M)')
axis([-2 30 0 1])


%npp
subplot(3,3,3)
scatter(log10(lme_npp(keep)/365),FracPD(keep),10,'k','filled'); hold on;
plot(npp(:,1),npp(:,2),'k'); hold on;
plot(npp(:,1),npp(:,2)+2*npp(:,3),'--k'); hold on;
plot(npp(:,1),npp(:,2)-2*npp(:,3),'--k'); hold on;
%title('Color = LME mean temp')
%ylabel('P / (P+D)')
axis([-0.5 1 0 1])

subplot(3,3,6)
scatter(log10(lme_npp(keep)/365),FracPF(keep),10,'k','filled'); hold on;
plot(npp(:,1),npp(:,4),'k'); hold on;
plot(npp(:,1),npp(:,4)+2*npp(:,5),'--k'); hold on;
plot(npp(:,1),npp(:,4)-2*npp(:,5),'--k'); hold on;
%ylabel('P / (P+F)')
axis([-0.5 1 0 1])

subplot(3,3,9)
scatter(log10(lme_npp(keep)/365),FracLM(keep),10,'k','filled'); hold on;
plot(npp(:,1),npp(:,6),'k'); hold on;
plot(npp(:,1),npp(:,6)+2*npp(:,7),'--k'); hold on;
plot(npp(:,1),npp(:,6)-2*npp(:,7),'--k'); hold on;
xlabel('log10 NPP')
%ylabel('L / (L+M)')
axis([-0.5 1 0 1])
print('-dpng',[ppath 'US_lme_scatter_ZlDet_pelT_NPP_BW_points.png'])

%% Just P:D, vary font size
figure(5)
plot(ZlDet(:,1),ZlDet(:,2),'LineWidth',3,'color',[0 0.5 0.75]); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k','LineWidth',2); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k','LineWidth',2); hold on;
%plot(log10(RatZlDet(keep)),FracPD,'k.','MarkerSize',25); hold on;
scatter(log10(RatZlDet(keep)),FracPD(keep),50,'k','filled'); hold on;
axis([-1 0.75 0 1])
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
set(gca,'FontSize',30,'XTick',[-1,-0.5,0,0.5],'XTickLabel',...
    {'-1.0','-0.5','0.0','0.5'},...
    'YTick',0:0.5:1,'YTickLabel',0:0.5:1)
print('-dpng',[ppath 'US_lme_scatter_ZlDet_GAMfit_PD_font30.png'])
%%
figure(6)
subplot(2,2,1)
scatter(log10(RatZlDet(keep)),FracPD(keep),20,'k','filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
axis([-1 0.75 0 1])
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
set(gca,'FontSize',14)

subplot(2,2,2)
scatter(log10(RatZlDet(keep)),FracPD(keep),20,'k','filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
axis([-1 0.75 0 1])
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
set(gca,'FontSize',14)

subplot(2,2,3)
scatter(log10(RatZlDet(keep)),FracPD(keep),20,'k','filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
axis([-1 0.75 0 1])
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
set(gca,'FontSize',14)

subplot(2,2,4)
scatter(log10(RatZlDet(keep)),FracPD(keep),20,'k','filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
axis([-1 0.75 0 1])
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
set(gca,'FontSize',14)
print('-dpng',[ppath 'US_lme_scatter_ZlDet_GAMfit_PD_sub4.png'])

