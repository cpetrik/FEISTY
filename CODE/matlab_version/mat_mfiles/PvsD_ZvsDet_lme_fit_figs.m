% Plot mean biomasses of zoop to compare to POEM results
% COMPARE BY LME, COLOR-CODE BY TEMP, Plot GAM fit

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
cfile = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
harv = ['All_fish',tfish(2:end)];
tharv = ['Harvest all fish ' num2str(frate)];
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'LME_ZBratios_clim_fished_',harv,'_' cfile '.mat']);
load([fpath 'ZlDet_gam_fits.mat']);

%% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);
close all

%% Scatter plot
% Plot GAM fit
figure(1)
subplot(3,2,1)
for i=1:66
    plot(log10(RatZlDet(i)),FracPD(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(PDgam(:,3),PDgam(:,1),'k'); hold on;
plot(PDgam(:,3),PDgam(:,1)+2*PDgam(:,2),'--k'); hold on;
plot(PDgam(:,3),PDgam(:,1)-2*PDgam(:,2),'--k'); hold on;
%title('Color = LME mean temp')
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(3,2,3)
for i=1:66
    plot(log10(RatZlDet(i)),FracPF(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(PFgam(:,3),PFgam(:,1),'k'); hold on;
plot(PFgam(:,3),PFgam(:,1)+2*PFgam(:,2),'--k'); hold on;
plot(PFgam(:,3),PFgam(:,1)-2*PFgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
axis([-1 0.75 0 1])

subplot(3,2,5)
for i=1:66
    plot(log10(RatZlDet(i)),FracLM(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(LMgam(:,3),LMgam(:,1),'k'); hold on;
plot(LMgam(:,3),LMgam(:,1)+2*LMgam(:,2),'--k'); hold on;
plot(LMgam(:,3),LMgam(:,1)-2*LMgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
axis([-1 0.75 0 1])
print('-dpng',[ppath 'lme_scatter_ZlDet_GAMfit_colorT.png'])

%%
figure(2)
subplot(3,2,1)
plot(log10(RatZlDet),FracPD,'.k','MarkerSize',10); hold on;
plot(PDgam(:,3),PDgam(:,1),'k'); hold on;
plot(PDgam(:,3),PDgam(:,1)+2*PDgam(:,2),'--k'); hold on;
plot(PDgam(:,3),PDgam(:,1)-2*PDgam(:,2),'--k'); hold on;
%title('Color = LME mean temp')
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(3,2,3)
plot(log10(RatZlDet),FracPF,'.k','MarkerSize',10); hold on;
plot(PFgam(:,3),PFgam(:,1),'k'); hold on;
plot(PFgam(:,3),PFgam(:,1)+2*PFgam(:,2),'--k'); hold on;
plot(PFgam(:,3),PFgam(:,1)-2*PFgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
axis([-1 0.75 0 1])

subplot(3,2,5)
plot(log10(RatZlDet),FracLM,'.k','MarkerSize',10); hold on;
plot(LMgam(:,3),LMgam(:,1),'k'); hold on;
plot(LMgam(:,3),LMgam(:,1)+2*LMgam(:,2),'--k'); hold on;
plot(LMgam(:,3),LMgam(:,1)-2*LMgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
axis([-1 0.75 0 1])
print('-dpng',[ppath 'lme_scatter_ZlDet_GAMfit_BW.png'])


%% Just P:D
figure(5)
subplot(2,2,3)
plot(log10(RatZlDet),FracPD,'.k','MarkerSize',10); hold on;
plot(PDgam(:,3),PDgam(:,1),'k'); hold on;
plot(PDgam(:,3),PDgam(:,1)+2*PDgam(:,2),'--k'); hold on;
plot(PDgam(:,3),PDgam(:,1)-2*PDgam(:,2),'--k'); hold on;
%title('Color = LME mean temp')
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(2,2,1)
for i=1:66
    plot(log10(RatZlDet(i)),FracPD(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(PDgam(:,3),PDgam(:,1),'k'); hold on;
plot(PDgam(:,3),PDgam(:,1)+2*PDgam(:,2),'--k'); hold on;
plot(PDgam(:,3),PDgam(:,1)-2*PDgam(:,2),'--k'); hold on;
title('Color = LME mean temp')
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])
print('-dpng',[ppath 'lme_scatter_ZlDet_GAMfit_PD.png'])

%% Just fits
figure(6)
subplot(3,2,1)
plot(PDgam(:,3),PDgam(:,1),'k'); hold on;
plot(PDgam(:,3),PDgam(:,1)+2*PDgam(:,2),'--k'); hold on;
plot(PDgam(:,3),PDgam(:,1)-2*PDgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(3,2,3)
plot(PFgam(:,3),PFgam(:,1),'k'); hold on;
plot(PFgam(:,3),PFgam(:,1)+2*PFgam(:,2),'--k'); hold on;
plot(PFgam(:,3),PFgam(:,1)-2*PFgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
axis([-1 0.75 0 1])

subplot(3,2,5)
plot(LMgam(:,3),LMgam(:,1),'k'); hold on;
plot(LMgam(:,3),LMgam(:,1)+2*LMgam(:,2),'--k'); hold on;
plot(LMgam(:,3),LMgam(:,1)-2*LMgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
axis([-1 0.75 0 1])
print('-dpng',[ppath 'lme_scatter_ZlDet_GAMfit_noData.png'])


%% WITH OTHER GAMS
% Plot GAM fit color
figure(7)
subplot(3,4,1)
for i=1:66
    plot(log10(RatZlDet(i)),FracPD(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(PDgam(:,3),PDgam(:,1),'k'); hold on;
plot(PDgam(:,3),PDgam(:,1)+2*PDgam(:,2),'--k'); hold on;
plot(PDgam(:,3),PDgam(:,1)-2*PDgam(:,2),'--k'); hold on;
%title('Color = LME mean temp')
%xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(3,4,5)
for i=1:66
    plot(log10(RatZlDet(i)),FracPF(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(PFgam(:,3),PFgam(:,1),'k'); hold on;
plot(PFgam(:,3),PFgam(:,1)+2*PFgam(:,2),'--k'); hold on;
plot(PFgam(:,3),PFgam(:,1)-2*PFgam(:,2),'--k'); hold on;
%xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
axis([-1 0.75 0 1])

subplot(3,4,9)
for i=1:66
    plot(log10(RatZlDet(i)),FracLM(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(LMgam(:,3),LMgam(:,1),'k'); hold on;
plot(LMgam(:,3),LMgam(:,1)+2*LMgam(:,2),'--k'); hold on;
plot(LMgam(:,3),LMgam(:,1)-2*LMgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
axis([-1 0.75 0 1])

subplot(3,4,2)
plot(PDtemp(:,3),PDtemp(:,1),'k'); hold on;
plot(PDtemp(:,3),PDtemp(:,1)+2*PDtemp(:,2),'--k'); hold on;
plot(PDtemp(:,3),PDtemp(:,1)-2*PDtemp(:,2),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+D)')
axis([-2 30 0 1])

subplot(3,4,6)
plot(PFtemp(:,3),PFtemp(:,1),'k'); hold on;
plot(PFtemp(:,3),PFtemp(:,1)+2*PFtemp(:,2),'--k'); hold on;
plot(PFtemp(:,3),PFtemp(:,1)-2*PFtemp(:,2),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+F)')
axis([-2 30 0 1])

subplot(3,4,10)
plot(LMtemp(:,3),LMtemp(:,1),'k'); hold on;
plot(LMtemp(:,3),LMtemp(:,1)+2*LMtemp(:,2),'--k'); hold on;
plot(LMtemp(:,3),LMtemp(:,1)-2*LMtemp(:,2),'--k'); hold on;
xlabel('Tpel (^oC)')
%ylabel('L / (L+M)')
axis([-2 30 0 1])

subplot(3,4,3)
plot(PDfrac(:,3),PDfrac(:,1),'k'); hold on;
plot(PDfrac(:,3),PDfrac(:,1)+2*PDfrac(:,2),'--k'); hold on;
plot(PDfrac(:,3),PDfrac(:,1)-2*PDfrac(:,2),'--k'); hold on;
%title('Color = LME mean temp')
%xlabel('Frac<200m')
%ylabel('P / (P+D)')
axis([0 1 0 1])

subplot(3,4,7)
plot(PFfrac(:,3),PFfrac(:,1),'k'); hold on;
plot(PFfrac(:,3),PFfrac(:,1)+2*PFfrac(:,2),'--k'); hold on;
plot(PFfrac(:,3),PFfrac(:,1)-2*PFfrac(:,2),'--k'); hold on;
%xlabel('Frac<200m')
%ylabel('P / (P+F)')
axis([0 1 0 1])

subplot(3,4,11)
plot(LMfrac(:,3),LMfrac(:,1),'k'); hold on;
plot(LMfrac(:,3),LMfrac(:,1)+2*LMfrac(:,2),'--k'); hold on;
plot(LMfrac(:,3),LMfrac(:,1)-2*LMfrac(:,2),'--k'); hold on;
xlabel('Frac<200m')
%ylabel('L / (L+M)')
axis([0 1 0 1])

subplot(3,4,4)
plot(PDnpp(:,3),PDnpp(:,1),'k'); hold on;
plot(PDnpp(:,3),PDnpp(:,1)+2*PDnpp(:,2),'--k'); hold on;
plot(PDnpp(:,3),PDnpp(:,1)-2*PDnpp(:,2),'--k'); hold on;
%title('Color = LME mean temp')
%ylabel('P / (P+D)')
axis([-0.5 1 0 1])

subplot(3,4,8)
plot(PFnpp(:,3),PFnpp(:,1),'k'); hold on;
plot(PFnpp(:,3),PFnpp(:,1)+2*PFnpp(:,2),'--k'); hold on;
plot(PFnpp(:,3),PFnpp(:,1)-2*PFnpp(:,2),'--k'); hold on;
%ylabel('P / (P+F)')
axis([-0.5 1 0 1])

subplot(3,4,12)
plot(LMnpp(:,3),LMnpp(:,1),'k'); hold on;
plot(LMnpp(:,3),LMnpp(:,1)+2*LMnpp(:,2),'--k'); hold on;
plot(LMnpp(:,3),LMnpp(:,1)-2*LMnpp(:,2),'--k'); hold on;
xlabel('log10 NPP')
%ylabel('L / (L+M)')
axis([-0.5 1 0 1])
print('-dpng',[ppath 'lme_scatter_ZlDet_allGAMfit_colorT.png'])

%%
figure(10)
subplot(3,4,1)
plot(log10(RatZlDet),FracPD,'.k','MarkerSize',10); hold on;
plot(PDgam(:,3),PDgam(:,1),'k'); hold on;
plot(PDgam(:,3),PDgam(:,1)+2*PDgam(:,2),'--k'); hold on;
plot(PDgam(:,3),PDgam(:,1)-2*PDgam(:,2),'--k'); hold on;
%title('Color = LME mean temp')
%xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(3,4,5)
plot(log10(RatZlDet),FracPF,'.k','MarkerSize',10); hold on;
plot(PFgam(:,3),PFgam(:,1),'k'); hold on;
plot(PFgam(:,3),PFgam(:,1)+2*PFgam(:,2),'--k'); hold on;
plot(PFgam(:,3),PFgam(:,1)-2*PFgam(:,2),'--k'); hold on;
%xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
axis([-1 0.75 0 1])

subplot(3,4,9)
plot(log10(RatZlDet),FracLM,'.k','MarkerSize',10); hold on;
plot(LMgam(:,3),LMgam(:,1),'k'); hold on;
plot(LMgam(:,3),LMgam(:,1)+2*LMgam(:,2),'--k'); hold on;
plot(LMgam(:,3),LMgam(:,1)-2*LMgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
axis([-1 0.75 0 1])

subplot(3,4,2)
plot(lme_ptemp(:,1),FracPD,'.k','MarkerSize',10); hold on;
plot(PDtemp(:,3),PDtemp(:,1),'k'); hold on;
plot(PDtemp(:,3),PDtemp(:,1)+2*PDtemp(:,2),'--k'); hold on;
plot(PDtemp(:,3),PDtemp(:,1)-2*PDtemp(:,2),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+D)')
axis([-2 30 0 1])

subplot(3,4,6)
plot(lme_ptemp(:,1),FracPF,'.k','MarkerSize',10); hold on;
plot(PFtemp(:,3),PFtemp(:,1),'k'); hold on;
plot(PFtemp(:,3),PFtemp(:,1)+2*PFtemp(:,2),'--k'); hold on;
plot(PFtemp(:,3),PFtemp(:,1)-2*PFtemp(:,2),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+F)')
axis([-2 30 0 1])

subplot(3,4,10)
plot(lme_ptemp(:,1),FracLM,'.k','MarkerSize',10); hold on;
plot(LMtemp(:,3),LMtemp(:,1),'k'); hold on;
plot(LMtemp(:,3),LMtemp(:,1)+2*LMtemp(:,2),'--k'); hold on;
plot(LMtemp(:,3),LMtemp(:,1)-2*LMtemp(:,2),'--k'); hold on;
xlabel('Tpel (^oC)')
%ylabel('L / (L+M)')
axis([-2 30 0 1])

subplot(3,4,3)
plot(lme_shal_frac(:,1),FracPD,'.k','MarkerSize',10); hold on;
plot(PDfrac(:,3),PDfrac(:,1),'k'); hold on;
plot(PDfrac(:,3),PDfrac(:,1)+2*PDfrac(:,2),'--k'); hold on;
plot(PDfrac(:,3),PDfrac(:,1)-2*PDfrac(:,2),'--k'); hold on;
%title('Color = LME mean temp')
%xlabel('Frac<200m')
%ylabel('P / (P+D)')
axis([0 1 0 1])

subplot(3,4,7)
plot(lme_shal_frac(:,1),FracPF,'.k','MarkerSize',10); hold on;
plot(PFfrac(:,3),PFfrac(:,1),'k'); hold on;
plot(PFfrac(:,3),PFfrac(:,1)+2*PFfrac(:,2),'--k'); hold on;
plot(PFfrac(:,3),PFfrac(:,1)-2*PFfrac(:,2),'--k'); hold on;
%xlabel('Frac<200m')
%ylabel('P / (P+F)')
axis([0 1 0 1])

subplot(3,4,11)
plot(lme_shal_frac(:,1),FracLM,'.k','MarkerSize',10); hold on;
plot(LMfrac(:,3),LMfrac(:,1),'k'); hold on;
plot(LMfrac(:,3),LMfrac(:,1)+2*LMfrac(:,2),'--k'); hold on;
plot(LMfrac(:,3),LMfrac(:,1)-2*LMfrac(:,2),'--k'); hold on;
xlabel('Frac<200m')
%ylabel('L / (L+M)')
axis([0 1 0 1])

subplot(3,4,4)
plot(log10(lme_npp/365),FracPD,'.k','MarkerSize',10); hold on;
plot(PDnpp(:,3),PDnpp(:,1),'k'); hold on;
plot(PDnpp(:,3),PDnpp(:,1)+2*PDnpp(:,2),'--k'); hold on;
plot(PDnpp(:,3),PDnpp(:,1)-2*PDnpp(:,2),'--k'); hold on;
%title('Color = LME mean temp')
%ylabel('P / (P+D)')
axis([-0.5 1 0 1])

subplot(3,4,8)
plot(log10(lme_npp/365),FracPF,'.k','MarkerSize',10); hold on;
plot(PFnpp(:,3),PFnpp(:,1),'k'); hold on;
plot(PFnpp(:,3),PFnpp(:,1)+2*PFnpp(:,2),'--k'); hold on;
plot(PFnpp(:,3),PFnpp(:,1)-2*PFnpp(:,2),'--k'); hold on;
%ylabel('P / (P+F)')
axis([-0.5 1 0 1])

subplot(3,4,12)
plot(log10(lme_npp/365),FracLM,'.k','MarkerSize',10); hold on;
plot(LMnpp(:,3),LMnpp(:,1),'k'); hold on;
plot(LMnpp(:,3),LMnpp(:,1)+2*LMnpp(:,2),'--k'); hold on;
plot(LMnpp(:,3),LMnpp(:,1)-2*LMnpp(:,2),'--k'); hold on;
xlabel('log10 NPP')
%ylabel('L / (L+M)')
axis([-0.5 1 0 1])
print('-dpng',[ppath 'lme_scatter_ZlDet_allGAMfit_BW_points.png'])

%% B&W
figure(9)
subplot(3,4,1)
%plot(log10(RatZlDet),FracPD,'.k','MarkerSize',10); hold on;
plot(PDgam(:,3),PDgam(:,1),'k'); hold on;
plot(PDgam(:,3),PDgam(:,1)+2*PDgam(:,2),'--k'); hold on;
plot(PDgam(:,3),PDgam(:,1)-2*PDgam(:,2),'--k'); hold on;
%title('Color = LME mean temp')
%xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])
set(gca','YTick',0:0.5:1)
title('A')

subplot(3,4,5)
%plot(log10(RatZlDet),FracPF,'.k','MarkerSize',10); hold on;
plot(PFgam(:,3),PFgam(:,1),'k'); hold on;
plot(PFgam(:,3),PFgam(:,1)+2*PFgam(:,2),'--k'); hold on;
plot(PFgam(:,3),PFgam(:,1)-2*PFgam(:,2),'--k'); hold on;
%xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
axis([-1 0.75 0 1])
set(gca,'YTick',0:0.5:1)
title('E')

subplot(3,4,9)
%plot(log10(RatZlDet),FracLM,'.k','MarkerSize',10); hold on;
plot(LMgam(:,3),LMgam(:,1),'k'); hold on;
plot(LMgam(:,3),LMgam(:,1)+2*LMgam(:,2),'--k'); hold on;
plot(LMgam(:,3),LMgam(:,1)-2*LMgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
axis([-1 0.75 0 1])
set(gca,'YTick',0:0.5:1)
title('I')

subplot(3,4,2)
plot(PDtemp(:,3),PDtemp(:,1),'k'); hold on;
plot(PDtemp(:,3),PDtemp(:,1)+2*PDtemp(:,2),'--k'); hold on;
plot(PDtemp(:,3),PDtemp(:,1)-2*PDtemp(:,2),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+D)')
axis([-2 30 0 1])
title('B')

subplot(3,4,6)
plot(PFtemp(:,3),PFtemp(:,1),'k'); hold on;
plot(PFtemp(:,3),PFtemp(:,1)+2*PFtemp(:,2),'--k'); hold on;
plot(PFtemp(:,3),PFtemp(:,1)-2*PFtemp(:,2),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+F)')
axis([-2 30 0 1])
title('F')

subplot(3,4,10)
plot(LMtemp(:,3),LMtemp(:,1),'k'); hold on;
plot(LMtemp(:,3),LMtemp(:,1)+2*LMtemp(:,2),'--k'); hold on;
plot(LMtemp(:,3),LMtemp(:,1)-2*LMtemp(:,2),'--k'); hold on;
xlabel('Tpel (^oC)')
%ylabel('L / (L+M)')
axis([-2 30 0 1])
title('J')

subplot(3,4,3)
plot(PDfrac(:,3),PDfrac(:,1),'k'); hold on;
plot(PDfrac(:,3),PDfrac(:,1)+2*PDfrac(:,2),'--k'); hold on;
plot(PDfrac(:,3),PDfrac(:,1)-2*PDfrac(:,2),'--k'); hold on;
%title('Color = LME mean temp')
%xlabel('Frac<200m')
%ylabel('P / (P+D)')
axis([0 1 0 1])
title('C')

subplot(3,4,7)
plot(PFfrac(:,3),PFfrac(:,1),'k'); hold on;
plot(PFfrac(:,3),PFfrac(:,1)+2*PFfrac(:,2),'--k'); hold on;
plot(PFfrac(:,3),PFfrac(:,1)-2*PFfrac(:,2),'--k'); hold on;
%xlabel('Frac<200m')
%ylabel('P / (P+F)')
axis([0 1 0 1])
title('G')

subplot(3,4,11)
plot(LMfrac(:,3),LMfrac(:,1),'k'); hold on;
plot(LMfrac(:,3),LMfrac(:,1)+2*LMfrac(:,2),'--k'); hold on;
plot(LMfrac(:,3),LMfrac(:,1)-2*LMfrac(:,2),'--k'); hold on;
xlabel('Frac<200m')
%ylabel('L / (L+M)')
axis([0 1 0 1])
title('K')

subplot(3,4,4)
plot(PDnpp(:,3),PDnpp(:,1),'k'); hold on;
plot(PDnpp(:,3),PDnpp(:,1)+2*PDnpp(:,2),'--k'); hold on;
plot(PDnpp(:,3),PDnpp(:,1)-2*PDnpp(:,2),'--k'); hold on;
%title('Color = LME mean temp')
%ylabel('P / (P+D)')
axis([-0.5 1 0 1])
title('D')

subplot(3,4,8)
plot(PFnpp(:,3),PFnpp(:,1),'k'); hold on;
plot(PFnpp(:,3),PFnpp(:,1)+2*PFnpp(:,2),'--k'); hold on;
plot(PFnpp(:,3),PFnpp(:,1)-2*PFnpp(:,2),'--k'); hold on;
%ylabel('P / (P+F)')
axis([-0.5 1 0 1])
title('H')

subplot(3,4,12)
plot(LMnpp(:,3),LMnpp(:,1),'k'); hold on;
plot(LMnpp(:,3),LMnpp(:,1)+2*LMnpp(:,2),'--k'); hold on;
plot(LMnpp(:,3),LMnpp(:,1)-2*LMnpp(:,2),'--k'); hold on;
xlabel('log10 NPP')
%ylabel('L / (L+M)')
axis([-0.5 1 0 1])
title('L')
print('-dpng',[ppath 'lme_scatter_ZlDet_allGAMfit_BW_fits.png'])

%% Color; ZB and T
figure(8)
subplot('position',[0.15 0.7 0.35 0.28])
for i=1:66
    plot(log10(RatZlDet(i)),FracPD(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(PDgam(:,3),PDgam(:,1),'k'); hold on;
plot(PDgam(:,3),PDgam(:,1)+2*PDgam(:,2),'--k'); hold on;
plot(PDgam(:,3),PDgam(:,1)-2*PDgam(:,2),'--k'); hold on;
%title('Color = LME mean temp')
%xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])
set(gca,'XTickLabel','','YTick',0:0.5:1)

subplot('position',[0.15 0.39 0.35 0.28])
for i=1:66
    plot(log10(RatZlDet(i)),FracPF(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(PFgam(:,3),PFgam(:,1),'k'); hold on;
plot(PFgam(:,3),PFgam(:,1)+2*PFgam(:,2),'--k'); hold on;
plot(PFgam(:,3),PFgam(:,1)-2*PFgam(:,2),'--k'); hold on;
%xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
axis([-1 0.75 0 1])
set(gca,'XTickLabel','','YTick',0:0.5:1)

subplot('position',[0.15 0.08 0.35 0.28])
for i=1:66
    plot(log10(RatZlDet(i)),FracLM(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(LMgam(:,3),LMgam(:,1),'k'); hold on;
plot(LMgam(:,3),LMgam(:,1)+2*LMgam(:,2),'--k'); hold on;
plot(LMgam(:,3),LMgam(:,1)-2*LMgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
axis([-1 0.75 0 1])
set(gca,'YTick',0:0.5:1)

subplot('position',[0.55 0.7 0.35 0.28])
plot(PDtemp(:,3),PDtemp(:,1),'k'); hold on;
plot(PDtemp(:,3),PDtemp(:,1)+2*PDtemp(:,2),'--k'); hold on;
plot(PDtemp(:,3),PDtemp(:,1)-2*PDtemp(:,2),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+D)')
axis([-2 30 0 1])
set(gca,'XTickLabel','','YTick',0:0.5:1)

subplot('position',[0.55 0.39 0.35 0.28])
plot(PFtemp(:,3),PFtemp(:,1),'k'); hold on;
plot(PFtemp(:,3),PFtemp(:,1)+2*PFtemp(:,2),'--k'); hold on;
plot(PFtemp(:,3),PFtemp(:,1)-2*PFtemp(:,2),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+F)')
axis([-2 30 0 1])
set(gca,'XTickLabel','','YTick',0:0.5:1)

subplot('position',[0.55 0.08 0.35 0.28])
plot(LMtemp(:,3),LMtemp(:,1),'k'); hold on;
plot(LMtemp(:,3),LMtemp(:,1)+2*LMtemp(:,2),'--k'); hold on;
plot(LMtemp(:,3),LMtemp(:,1)-2*LMtemp(:,2),'--k'); hold on;
xlabel('Tpel (^oC)')
%ylabel('L / (L+M)')
axis([-2 30 0 1])
set(gca,'YTick',0:0.5:1)
print('-dpng',[ppath 'lme_scatter_ZlDet_GAMfit_colorT_v3.png'])

