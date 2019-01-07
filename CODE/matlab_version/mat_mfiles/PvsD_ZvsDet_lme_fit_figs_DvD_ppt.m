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
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];

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

%% Scatter plot
%
figure(1)
subplot(3,3,1)
scatter(log10(RatZlDet),FracPD,10,lme_ptemp(:,1),'filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
%title('Color = LME mean temp')
cmocean('thermal');
%colorbar
% xlabel('log10 ZoopLoss:Det')
% ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(3,3,4)
scatter(log10(RatZlDet),FracPF,10,lme_ptemp(:,1),'filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,4),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)+2*ZlDet(:,5),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)-2*ZlDet(:,5),'--k'); hold on;
cmocean('thermal');
colorbar('Position',[0.375 0.4 0.02 0.25])
% xlabel('log10 ZoopLoss:Det')
% ylabel('P / (P+F)')
axis([-1 0.75 0 1])

subplot(3,3,7)
scatter(log10(RatZlDet),FracLM,10,lme_ptemp(:,1),'filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,6),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)+2*ZlDet(:,7),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)-2*ZlDet(:,7),'--k'); hold on;
cmocean('thermal');
%colorbar
% xlabel('log10 ZoopLoss:Det')
% ylabel('L / (L+M)')
axis([-1 0.75 0 1])
print('-dpng',[ppath 'lme_scatter_ZlDet_GAMfit_colorT_colorbar.png'])


%% Just P:D
figure(2)
subplot(2,2,3)
plot(log10(RatZlDet),FracPD,'k.','MarkerSize',10); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
%title('Color = LME mean temp')
% xlabel('log10 ZoopLoss:Det')
% ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(2,2,1)
scatter(log10(RatZlDet),FracPD,15,lme_ptemp(:,1),'filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
cmocean('thermal');
% title('Color = LME mean temp')
% xlabel('log10 ZoopLoss:Det')
% ylabel('P / (P+D)')
axis([-1 0.75 0 1])
print('-dpng',[ppath 'lme_scatter_ZlDet_GAMfit_PD.png'])

%% Just fits
figure(3)
subplot(3,3,1)
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
% xlabel('log10 ZoopLoss:Det')
% ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(3,3,4)
plot(ZlDet(:,1),ZlDet(:,4),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)+2*ZlDet(:,5),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)-2*ZlDet(:,5),'--k'); hold on;
% xlabel('log10 ZoopLoss:Det')
% ylabel('P / (P+F)')
axis([-1 0.75 0 1])

subplot(3,3,7)
plot(ZlDet(:,1),ZlDet(:,6),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)+2*ZlDet(:,7),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)-2*ZlDet(:,7),'--k'); hold on;
% xlabel('log10 ZoopLoss:Det')
% ylabel('L / (L+M)')
axis([-1 0.75 0 1])
print('-dpng',[ppath 'lme_scatter_ZlDet_GAMfit_noData.png'])


%% WITH OTHER GAMS
% Plot GAM fit color
figure(4)
subplot(3,3,1)
scatter(log10(RatZlDet),FracPD,10,lme_ptemp(:,1),'filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
cmocean('thermal');
%title('Color = LME mean temp')
%xlabel('log10 ZoopLoss:Det')
% ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(3,3,4)
scatter(log10(RatZlDet),FracPF,10,lme_ptemp(:,1),'filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,4),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)+2*ZlDet(:,5),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)-2*ZlDet(:,5),'--k'); hold on;
cmocean('thermal');
%xlabel('log10 ZoopLoss:Det')
% ylabel('P / (P+F)')
axis([-1 0.75 0 1])

subplot(3,3,7)
scatter(log10(RatZlDet),FracLM,10,lme_ptemp(:,1),'filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,6),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)+2*ZlDet(:,7),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)-2*ZlDet(:,7),'--k'); hold on;
cmocean('thermal');
% xlabel('log10 ZoopLoss:Det')
% ylabel('L / (L+M)')
axis([-1 0.75 0 1])

%ptemp
subplot(3,3,2)
plot(ptemp(:,1),ptemp(:,2),'k'); hold on;
plot(ptemp(:,1),ptemp(:,2)+2*ptemp(:,3),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,2)-2*ptemp(:,3),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+D)')
axis([-2 30 0 1])

subplot(3,3,5)
plot(ptemp(:,1),ptemp(:,4),'k'); hold on;
plot(ptemp(:,1),ptemp(:,4)+2*ptemp(:,5),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,4)-2*ptemp(:,5),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+F)')
axis([-2 30 0 1])

subplot(3,3,8)
plot(ptemp(:,1),ptemp(:,6),'k'); hold on;
plot(ptemp(:,1),ptemp(:,6)+2*ptemp(:,7),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,6)-2*ptemp(:,7),'--k'); hold on;
% xlabel('Tpel (^oC)')
%ylabel('L / (L+M)')
axis([-2 30 0 1])


%npp
subplot(3,3,3)
plot(npp(:,1),npp(:,2),'k'); hold on;
plot(npp(:,1),npp(:,2)+2*npp(:,3),'--k'); hold on;
plot(npp(:,1),npp(:,2)-2*npp(:,3),'--k'); hold on;
%title('Color = LME mean temp')
%ylabel('P / (P+D)')
axis([-0.5 1 0 1])

subplot(3,3,6)
plot(npp(:,1),npp(:,4),'k'); hold on;
plot(npp(:,1),npp(:,4)+2*npp(:,5),'--k'); hold on;
plot(npp(:,1),npp(:,4)-2*npp(:,5),'--k'); hold on;
%ylabel('P / (P+F)')
axis([-0.5 1 0 1])

subplot(3,3,9)
plot(npp(:,1),npp(:,6),'k'); hold on;
plot(npp(:,1),npp(:,6)+2*npp(:,7),'--k'); hold on;
plot(npp(:,1),npp(:,6)-2*npp(:,7),'--k'); hold on;
% xlabel('log10 NPP')
%ylabel('L / (L+M)')
axis([-0.5 1 0 1])
print('-dpng',[ppath 'lme_scatter_ZlDet_pelT_NPP_colorT.png'])

%%
figure(5)
subplot(3,3,1)
plot(log10(RatZlDet),FracPD,'.k','MarkerSize',10); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
%title('Color = LME mean temp')
%xlabel('log10 ZoopLoss:Det')
% ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(3,3,4)
plot(log10(RatZlDet),FracPF,'.k','MarkerSize',10); hold on;
plot(ZlDet(:,1),ZlDet(:,4),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)+2*ZlDet(:,5),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)-2*ZlDet(:,5),'--k'); hold on;
%xlabel('log10 ZoopLoss:Det')
% ylabel('P / (P+F)')
axis([-1 0.75 0 1])

subplot(3,3,7)
plot(log10(RatZlDet),FracLM,'.k','MarkerSize',10); hold on;
plot(ZlDet(:,1),ZlDet(:,6),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)+2*ZlDet(:,7),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)-2*ZlDet(:,7),'--k'); hold on;
% xlabel('log10 ZoopLoss:Det')
% ylabel('L / (L+M)')
axis([-1 0.75 0 1])

%ptemp
subplot(3,3,2)
plot(lme_ptemp(:,1),FracPD,'.k','MarkerSize',10); hold on;
plot(ptemp(:,1),ptemp(:,2),'k'); hold on;
plot(ptemp(:,1),ptemp(:,2)+2*ptemp(:,3),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,2)-2*ptemp(:,3),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+D)')
axis([-2 30 0 1])

subplot(3,3,5)
plot(lme_ptemp(:,1),FracPF,'.k','MarkerSize',10); hold on;
plot(ptemp(:,1),ptemp(:,4),'k'); hold on;
plot(ptemp(:,1),ptemp(:,4)+2*ptemp(:,5),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,4)-2*ptemp(:,5),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+F)')
axis([-2 30 0 1])

subplot(3,3,8)
plot(lme_ptemp(:,1),FracLM,'.k','MarkerSize',10); hold on;
plot(ptemp(:,1),ptemp(:,6),'k'); hold on;
plot(ptemp(:,1),ptemp(:,6)+2*ptemp(:,7),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,6)-2*ptemp(:,7),'--k'); hold on;
% xlabel('Tpel (^oC)')
%ylabel('L / (L+M)')
axis([-2 30 0 1])


%npp
subplot(3,3,3)
plot(log10(lme_npp/365),FracPD,'.k','MarkerSize',10); hold on;
plot(npp(:,1),npp(:,2),'k'); hold on;
plot(npp(:,1),npp(:,2)+2*npp(:,3),'--k'); hold on;
plot(npp(:,1),npp(:,2)-2*npp(:,3),'--k'); hold on;
%title('Color = LME mean temp')
%ylabel('P / (P+D)')
axis([-0.5 1 0 1])

subplot(3,3,6)
plot(log10(lme_npp/365),FracPF,'.k','MarkerSize',10); hold on;
plot(npp(:,1),npp(:,4),'k'); hold on;
plot(npp(:,1),npp(:,4)+2*npp(:,5),'--k'); hold on;
plot(npp(:,1),npp(:,4)-2*npp(:,5),'--k'); hold on;
%ylabel('P / (P+F)')
axis([-0.5 1 0 1])

subplot(3,3,9)
plot(log10(lme_npp/365),FracLM,'.k','MarkerSize',10); hold on;
plot(npp(:,1),npp(:,6),'k'); hold on;
plot(npp(:,1),npp(:,6)+2*npp(:,7),'--k'); hold on;
plot(npp(:,1),npp(:,6)-2*npp(:,7),'--k'); hold on;
% xlabel('log10 NPP')
%ylabel('L / (L+M)')
axis([-0.5 1 0 1])
print('-dpng',[ppath 'lme_scatter_ZlDet_pelT_NPP_BW_points.png'])

%% B&W just fits
figure(6)
subplot(3,3,1)
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
%title('Color = LME mean temp')
%xlabel('log10 ZoopLoss:Det')
% ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(3,3,4)
plot(ZlDet(:,1),ZlDet(:,4),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)+2*ZlDet(:,5),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)-2*ZlDet(:,5),'--k'); hold on;
%xlabel('log10 ZoopLoss:Det')
% ylabel('P / (P+F)')
axis([-1 0.75 0 1])

subplot(3,3,7)
plot(ZlDet(:,1),ZlDet(:,6),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)+2*ZlDet(:,7),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)-2*ZlDet(:,7),'--k'); hold on;
% xlabel('log10 ZoopLoss:Det')
% ylabel('L / (L+M)')
axis([-1 0.75 0 1])

%ptemp
subplot(3,3,2)
plot(ptemp(:,1),ptemp(:,2),'k'); hold on;
plot(ptemp(:,1),ptemp(:,2)+2*ptemp(:,3),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,2)-2*ptemp(:,3),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+D)')
axis([-2 30 0 1])

subplot(3,3,5)
plot(ptemp(:,1),ptemp(:,4),'k'); hold on;
plot(ptemp(:,1),ptemp(:,4)+2*ptemp(:,5),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,4)-2*ptemp(:,5),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+F)')
axis([-2 30 0 1])

subplot(3,3,8)
plot(ptemp(:,1),ptemp(:,6),'k'); hold on;
plot(ptemp(:,1),ptemp(:,6)+2*ptemp(:,7),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,6)-2*ptemp(:,7),'--k'); hold on;
% xlabel('Tpel (^oC)')
%ylabel('L / (L+M)')
axis([-2 30 0 1])

%npp
subplot(3,3,3)
plot(npp(:,1),npp(:,2),'k'); hold on;
plot(npp(:,1),npp(:,2)+2*npp(:,3),'--k'); hold on;
plot(npp(:,1),npp(:,2)-2*npp(:,3),'--k'); hold on;
%title('Color = LME mean temp')
%ylabel('P / (P+D)')
axis([-0.5 1 0 1])

subplot(3,3,6)
plot(npp(:,1),npp(:,4),'k'); hold on;
plot(npp(:,1),npp(:,4)+2*npp(:,5),'--k'); hold on;
plot(npp(:,1),npp(:,4)-2*npp(:,5),'--k'); hold on;
%ylabel('P / (P+F)')
axis([-0.5 1 0 1])

subplot(3,3,9)
plot(npp(:,1),npp(:,6),'k'); hold on;
plot(npp(:,1),npp(:,6)+2*npp(:,7),'--k'); hold on;
plot(npp(:,1),npp(:,6)-2*npp(:,7),'--k'); hold on;
% xlabel('log10 NPP')
%ylabel('L / (L+M)')
axis([-0.5 1 0 1])
%text(-0.4,0.9,'L')
print('-dpng',[ppath 'lme_scatter_ZlDet_pelT_NPP_BW_fits.png'])



