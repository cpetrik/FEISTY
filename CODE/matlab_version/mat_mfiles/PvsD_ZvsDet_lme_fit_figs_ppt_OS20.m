% Plot mean biomasses of zoop to compare to POEM results
% COMPARE BY LME, COLOR-CODE BY TEMP, Plot GAM fit
% Uses 3 Knots and Probit link fn

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';

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
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
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

%% Scatter plot black
%
figure(1)
subplot(2,2,1)
scatter(log10(RatZlDet),FracPD,20,'k','filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
axis([-1 0.75 0 1])

subplot(2,2,2)
scatter(lme_ptemp(:,1),FracPD,20,'k','filled'); hold on;
plot(ptemp(:,1),ptemp(:,2),'k'); hold on;
plot(ptemp(:,1),ptemp(:,2)+2*ptemp(:,3),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,2)-2*ptemp(:,3),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+D)')
axis([-2 30 0 1])
print('-dpng',[ppath 'lme_scatter_ZlDet_GAMfit_black_FracPD_temp.png'])

figure(2)
subplot(2,2,1)
scatter(log10(RatZlDet),FracPD,20,lme_ptemp(:,1),'filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
%title('Color = LME mean temp')
cmocean('thermal');
%colorbar
% xlabel('log10 ZoopLoss:Det')
% ylabel('P / (P+D)')
axis([-1 0.75 0 1])

%ptemp
subplot(2,2,2)
scatter(lme_ptemp(:,1),FracPD,20,'k','filled'); hold on;
plot(ptemp(:,1),ptemp(:,2),'k'); hold on;
plot(ptemp(:,1),ptemp(:,2)+2*ptemp(:,3),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,2)-2*ptemp(:,3),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+D)')
axis([-2 30 0 1])
print('-dpng',[ppath 'lme_scatter_ZlDet_GAMfit_ptemp_FracPD_temp.png'])

