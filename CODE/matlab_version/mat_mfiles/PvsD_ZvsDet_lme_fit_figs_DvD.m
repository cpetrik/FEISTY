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
% Plot GAM fit
figure(1)
subplot(3,2,1)
for i=1:66
    plot(log10(RatZlDet(i)),FracPD(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
%title('Color = LME mean temp')
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])
text(-0.9,0.9,'A')

subplot(3,2,3)
for i=1:66
    plot(log10(RatZlDet(i)),FracPF(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(ZlDet(:,1),ZlDet(:,4),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)+2*ZlDet(:,5),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)-2*ZlDet(:,5),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
axis([-1 0.75 0 1])
text(-0.9,0.9,'B')

subplot(3,2,5)
for i=1:66
    plot(log10(RatZlDet(i)),FracLM(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(ZlDet(:,1),ZlDet(:,6),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)+2*ZlDet(:,7),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)-2*ZlDet(:,7),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
axis([-1 0.75 0 1])
text(-0.9,0.9,'C')
print('-dpng',[ppath 'lme_scatter_ZlDet_GAMfit_colorT.png'])

%%
figure(11)
subplot(3,2,1)
scatter(log10(RatZlDet),FracPD,10,lme_ptemp(:,1),'filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
%title('Color = LME mean temp')
cmocean('thermal');
%colorbar
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])
text(-0.9,0.9,'A')

subplot(3,2,3)
scatter(log10(RatZlDet),FracPF,10,lme_ptemp(:,1),'filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,4),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)+2*ZlDet(:,5),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)-2*ZlDet(:,5),'--k'); hold on;
cmocean('thermal');
colorbar('Position',[0.475 0.325 0.03 0.39])
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
axis([-1 0.75 0 1])
text(-0.9,0.9,'B')

subplot(3,2,5)
scatter(log10(RatZlDet),FracLM,10,lme_ptemp(:,1),'filled'); hold on;
plot(ZlDet(:,1),ZlDet(:,6),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)+2*ZlDet(:,7),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)-2*ZlDet(:,7),'--k'); hold on;
cmocean('thermal');
%colorbar
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
axis([-1 0.75 0 1])
text(-0.9,0.9,'C')
print('-dpng',[ppath 'lme_scatter_ZlDet_GAMfit_colorT_v2.png'])

%%
figure(2)
subplot(3,2,1)
plot(log10(RatZlDet),FracPD,'k.','MarkerSize',10); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
%title('Color = LME mean temp')
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])
text(-0.9,0.9,'A')

subplot(3,2,3)
plot(log10(RatZlDet),FracPF,'k.','MarkerSize',10); hold on;
plot(ZlDet(:,1),ZlDet(:,4),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)+2*ZlDet(:,5),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)-2*ZlDet(:,5),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
axis([-1 0.75 0 1])
text(-0.9,0.9,'B')

subplot(3,2,5)
plot(log10(RatZlDet),FracLM,'k.','MarkerSize',10); hold on;
plot(ZlDet(:,1),ZlDet(:,6),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)+2*ZlDet(:,7),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)-2*ZlDet(:,7),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
axis([-1 0.75 0 1])
text(-0.9,0.9,'C')
print('-dpng',[ppath 'lme_scatter_ZlDet_GAMfit_BW.png'])


%% Just P:D
figure(3)
subplot(2,2,3)
plot(log10(RatZlDet),FracPD,'k.','MarkerSize',10); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
%title('Color = LME mean temp')
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(2,2,1)
for i=1:66
    plot(log10(RatZlDet(i)),FracPD(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
title('Color = LME mean temp')
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])
print('-dpng',[ppath 'lme_scatter_ZlDet_GAMfit_PD.png'])

%% Just fits
figure(4)
subplot(3,2,1)
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(3,2,3)
plot(ZlDet(:,1),ZlDet(:,4),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)+2*ZlDet(:,5),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)-2*ZlDet(:,5),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
axis([-1 0.75 0 1])

subplot(3,2,5)
plot(ZlDet(:,1),ZlDet(:,6),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)+2*ZlDet(:,7),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)-2*ZlDet(:,7),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
axis([-1 0.75 0 1])
print('-dpng',[ppath 'lme_scatter_ZlDet_GAMfit_noData.png'])


%% WITH OTHER GAMS
% Plot GAM fit color
figure(5)
subplot(3,4,1)
for i=1:66
    plot(log10(RatZlDet(i)),FracPD(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
%title('Color = LME mean temp')
%xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(3,4,5)
for i=1:66
    plot(log10(RatZlDet(i)),FracPF(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(ZlDet(:,1),ZlDet(:,4),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)+2*ZlDet(:,5),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)-2*ZlDet(:,5),'--k'); hold on;
%xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
axis([-1 0.75 0 1])

subplot(3,4,9)
for i=1:66
    plot(log10(RatZlDet(i)),FracLM(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(ZlDet(:,1),ZlDet(:,6),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)+2*ZlDet(:,7),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)-2*ZlDet(:,7),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
axis([-1 0.75 0 1])

%ptemp
subplot(3,4,2)
plot(ptemp(:,1),ptemp(:,2),'k'); hold on;
plot(ptemp(:,1),ptemp(:,2)+2*ptemp(:,3),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,2)-2*ptemp(:,3),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+D)')
axis([-2 30 0 1])

subplot(3,4,6)
plot(ptemp(:,1),ptemp(:,4),'k'); hold on;
plot(ptemp(:,1),ptemp(:,4)+2*ptemp(:,5),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,4)-2*ptemp(:,5),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+F)')
axis([-2 30 0 1])

subplot(3,4,10)
plot(ptemp(:,1),ptemp(:,6),'k'); hold on;
plot(ptemp(:,1),ptemp(:,6)+2*ptemp(:,7),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,6)-2*ptemp(:,7),'--k'); hold on;
xlabel('Tpel (^oC)')
%ylabel('L / (L+M)')
axis([-2 30 0 1])

%Frac200
subplot(3,4,3)
plot(Frac200(:,1),Frac200(:,2),'k'); hold on;
plot(Frac200(:,1),Frac200(:,2)+2*Frac200(:,3),'--k'); hold on;
plot(Frac200(:,1),Frac200(:,2)-2*Frac200(:,3),'--k'); hold on;
%title('Color = LME mean temp')
%xlabel('Frac<200m')
%ylabel('P / (P+D)')
axis([0 1 0 1])

subplot(3,4,7)
plot(Frac200(:,1),Frac200(:,4),'k'); hold on;
plot(Frac200(:,1),Frac200(:,4)+2*Frac200(:,5),'--k'); hold on;
plot(Frac200(:,1),Frac200(:,4)-2*Frac200(:,5),'--k'); hold on;
%xlabel('Frac<200m')
%ylabel('P / (P+F)')
axis([0 1 0 1])

subplot(3,4,11)
plot(Frac200(:,1),Frac200(:,6),'k'); hold on;
plot(Frac200(:,1),Frac200(:,6)+2*Frac200(:,7),'--k'); hold on;
plot(Frac200(:,1),Frac200(:,6)-2*Frac200(:,7),'--k'); hold on;
xlabel('Frac<200m')
%ylabel('L / (L+M)')
axis([0 1 0 1])

%npp
subplot(3,4,4)
plot(npp(:,1),npp(:,2),'k'); hold on;
plot(npp(:,1),npp(:,2)+2*npp(:,3),'--k'); hold on;
plot(npp(:,1),npp(:,2)-2*npp(:,3),'--k'); hold on;
%title('Color = LME mean temp')
%ylabel('P / (P+D)')
axis([-0.5 1 0 1])

subplot(3,4,8)
plot(npp(:,1),npp(:,4),'k'); hold on;
plot(npp(:,1),npp(:,4)+2*npp(:,5),'--k'); hold on;
plot(npp(:,1),npp(:,4)-2*npp(:,5),'--k'); hold on;
%ylabel('P / (P+F)')
axis([-0.5 1 0 1])

subplot(3,4,12)
plot(npp(:,1),npp(:,6),'k'); hold on;
plot(npp(:,1),npp(:,6)+2*npp(:,7),'--k'); hold on;
plot(npp(:,1),npp(:,6)-2*npp(:,7),'--k'); hold on;
xlabel('log10 NPP')
%ylabel('L / (L+M)')
axis([-0.5 1 0 1])
print('-dpng',[ppath 'lme_scatter_ZlDet_allGAMfit_colorT.png'])

%%
figure(6)
subplot(3,4,1)
plot(log10(RatZlDet),FracPD,'.k','MarkerSize',10); hold on;
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
%title('Color = LME mean temp')
%xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])
text(-0.9,0.9,'A')

subplot(3,4,5)
plot(log10(RatZlDet),FracPF,'.k','MarkerSize',10); hold on;
plot(ZlDet(:,1),ZlDet(:,4),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)+2*ZlDet(:,5),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)-2*ZlDet(:,5),'--k'); hold on;
%xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
axis([-1 0.75 0 1])
text(-0.9,0.9,'E')

subplot(3,4,9)
plot(log10(RatZlDet),FracLM,'.k','MarkerSize',10); hold on;
plot(ZlDet(:,1),ZlDet(:,6),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)+2*ZlDet(:,7),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)-2*ZlDet(:,7),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
axis([-1 0.75 0 1])
text(-0.9,0.9,'I')

%ptemp
subplot(3,4,2)
plot(lme_ptemp(:,1),FracPD,'.k','MarkerSize',10); hold on;
plot(ptemp(:,1),ptemp(:,2),'k'); hold on;
plot(ptemp(:,1),ptemp(:,2)+2*ptemp(:,3),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,2)-2*ptemp(:,3),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+D)')
axis([-2 30 0 1])
text(-1,0.9,'B')

subplot(3,4,6)
plot(lme_ptemp(:,1),FracPF,'.k','MarkerSize',10); hold on;
plot(ptemp(:,1),ptemp(:,4),'k'); hold on;
plot(ptemp(:,1),ptemp(:,4)+2*ptemp(:,5),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,4)-2*ptemp(:,5),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+F)')
axis([-2 30 0 1])
text(-1,0.9,'F')

subplot(3,4,10)
plot(lme_ptemp(:,1),FracLM,'.k','MarkerSize',10); hold on;
plot(ptemp(:,1),ptemp(:,6),'k'); hold on;
plot(ptemp(:,1),ptemp(:,6)+2*ptemp(:,7),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,6)-2*ptemp(:,7),'--k'); hold on;
xlabel('Tpel (^oC)')
%ylabel('L / (L+M)')
axis([-2 30 0 1])
text(-1,0.9,'J')

%Frac200
subplot(3,4,3)
plot(lme_shal_frac(:,1),FracPD,'.k','MarkerSize',10); hold on;
plot(Frac200(:,1),Frac200(:,2),'k'); hold on;
plot(Frac200(:,1),Frac200(:,2)+2*Frac200(:,3),'--k'); hold on;
plot(Frac200(:,1),Frac200(:,2)-2*Frac200(:,3),'--k'); hold on;
%title('Color = LME mean temp')
%xlabel('Frac<200m')
%ylabel('P / (P+D)')
axis([0 1 0 1])
text(0.1,0.9,'C')

subplot(3,4,7)
plot(lme_shal_frac(:,1),FracPF,'.k','MarkerSize',10); hold on;
plot(Frac200(:,1),Frac200(:,4),'k'); hold on;
plot(Frac200(:,1),Frac200(:,4)+2*Frac200(:,5),'--k'); hold on;
plot(Frac200(:,1),Frac200(:,4)-2*Frac200(:,5),'--k'); hold on;
%xlabel('Frac<200m')
%ylabel('P / (P+F)')
axis([0 1 0 1])
text(0.1,0.9,'G')

subplot(3,4,11)
plot(lme_shal_frac(:,1),FracLM,'.k','MarkerSize',10); hold on;
plot(Frac200(:,1),Frac200(:,6),'k'); hold on;
plot(Frac200(:,1),Frac200(:,6)+2*Frac200(:,7),'--k'); hold on;
plot(Frac200(:,1),Frac200(:,6)-2*Frac200(:,7),'--k'); hold on;
xlabel('Frac<200m')
%ylabel('L / (L+M)')
axis([0 1 0 1])
text(0.1,0.9,'K')

%npp
subplot(3,4,4)
plot(log10(lme_npp/365),FracPD,'.k','MarkerSize',10); hold on;
plot(npp(:,1),npp(:,2),'k'); hold on;
plot(npp(:,1),npp(:,2)+2*npp(:,3),'--k'); hold on;
plot(npp(:,1),npp(:,2)-2*npp(:,3),'--k'); hold on;
%title('Color = LME mean temp')
%ylabel('P / (P+D)')
axis([-0.5 1 0 1])
text(-0.4,0.9,'D')

subplot(3,4,8)
plot(log10(lme_npp/365),FracPF,'.k','MarkerSize',10); hold on;
plot(npp(:,1),npp(:,4),'k'); hold on;
plot(npp(:,1),npp(:,4)+2*npp(:,5),'--k'); hold on;
plot(npp(:,1),npp(:,4)-2*npp(:,5),'--k'); hold on;
%ylabel('P / (P+F)')
axis([-0.5 1 0 1])
text(-0.4,0.9,'H')

subplot(3,4,12)
plot(log10(lme_npp/365),FracLM,'.k','MarkerSize',10); hold on;
plot(npp(:,1),npp(:,6),'k'); hold on;
plot(npp(:,1),npp(:,6)+2*npp(:,7),'--k'); hold on;
plot(npp(:,1),npp(:,6)-2*npp(:,7),'--k'); hold on;
xlabel('log10 NPP')
%ylabel('L / (L+M)')
axis([-0.5 1 0 1])
text(-0.4,0.9,'L')
print('-dpng',[ppath 'lme_scatter_ZlDet_allGAMfit_BW_points.png'])

%% B&W just fits
figure(7)
subplot(3,4,1)
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
%title('Color = LME mean temp')
%xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])

subplot(3,4,5)
plot(ZlDet(:,1),ZlDet(:,4),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)+2*ZlDet(:,5),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,4)-2*ZlDet(:,5),'--k'); hold on;
%xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
axis([-1 0.75 0 1])

subplot(3,4,9)
plot(ZlDet(:,1),ZlDet(:,6),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)+2*ZlDet(:,7),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,6)-2*ZlDet(:,7),'--k'); hold on;
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
axis([-1 0.75 0 1])

%ptemp
subplot(3,4,2)
plot(ptemp(:,1),ptemp(:,2),'k'); hold on;
plot(ptemp(:,1),ptemp(:,2)+2*ptemp(:,3),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,2)-2*ptemp(:,3),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+D)')
axis([-2 30 0 1])

subplot(3,4,6)
plot(ptemp(:,1),ptemp(:,4),'k'); hold on;
plot(ptemp(:,1),ptemp(:,4)+2*ptemp(:,5),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,4)-2*ptemp(:,5),'--k'); hold on;
%xlabel('Tpel (^oC)')
%ylabel('P / (P+F)')
axis([-2 30 0 1])

subplot(3,4,10)
plot(ptemp(:,1),ptemp(:,6),'k'); hold on;
plot(ptemp(:,1),ptemp(:,6)+2*ptemp(:,7),'--k'); hold on;
plot(ptemp(:,1),ptemp(:,6)-2*ptemp(:,7),'--k'); hold on;
xlabel('Tpel (^oC)')
%ylabel('L / (L+M)')
axis([-2 30 0 1])

%Frac200
subplot(3,4,3)
plot(Frac200(:,1),Frac200(:,2),'k'); hold on;
plot(Frac200(:,1),Frac200(:,2)+2*Frac200(:,3),'--k'); hold on;
plot(Frac200(:,1),Frac200(:,2)-2*Frac200(:,3),'--k'); hold on;
%title('Color = LME mean temp')
%xlabel('Frac<200m')
%ylabel('P / (P+D)')
axis([0 1 0 1])

subplot(3,4,7)
plot(Frac200(:,1),Frac200(:,4),'k'); hold on;
plot(Frac200(:,1),Frac200(:,4)+2*Frac200(:,5),'--k'); hold on;
plot(Frac200(:,1),Frac200(:,4)-2*Frac200(:,5),'--k'); hold on;
%xlabel('Frac<200m')
%ylabel('P / (P+F)')
axis([0 1 0 1])

subplot(3,4,11)
plot(Frac200(:,1),Frac200(:,6),'k'); hold on;
plot(Frac200(:,1),Frac200(:,6)+2*Frac200(:,7),'--k'); hold on;
plot(Frac200(:,1),Frac200(:,6)-2*Frac200(:,7),'--k'); hold on;
xlabel('Frac<200m')
%ylabel('L / (L+M)')
axis([0 1 0 1])

%npp
subplot(3,4,4)
plot(npp(:,1),npp(:,2),'k'); hold on;
plot(npp(:,1),npp(:,2)+2*npp(:,3),'--k'); hold on;
plot(npp(:,1),npp(:,2)-2*npp(:,3),'--k'); hold on;
%title('Color = LME mean temp')
%ylabel('P / (P+D)')
axis([-0.5 1 0 1])

subplot(3,4,8)
plot(npp(:,1),npp(:,4),'k'); hold on;
plot(npp(:,1),npp(:,4)+2*npp(:,5),'--k'); hold on;
plot(npp(:,1),npp(:,4)-2*npp(:,5),'--k'); hold on;
%ylabel('P / (P+F)')
axis([-0.5 1 0 1])

subplot(3,4,12)
plot(npp(:,1),npp(:,6),'k'); hold on;
plot(npp(:,1),npp(:,6)+2*npp(:,7),'--k'); hold on;
plot(npp(:,1),npp(:,6)-2*npp(:,7),'--k'); hold on;
xlabel('log10 NPP')
%ylabel('L / (L+M)')
axis([-0.5 1 0 1])
%text(-0.4,0.9,'L')
print('-dpng',[ppath 'lme_scatter_ZlDet_allGAMfit_BW_fits.png'])


%% WITH NPP GAM
% Plot GAM fit color
figure(8)
subplot(2,2,1)
for i=1:66
    plot(log10(RatZlDet(i)),FracPD(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(ZlDet(:,1),ZlDet(:,2),'k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)+2*ZlDet(:,3),'--k'); hold on;
plot(ZlDet(:,1),ZlDet(:,2)-2*ZlDet(:,3),'--k'); hold on;
%title('Color = LME mean temp')
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
axis([-1 0.75 0 1])

%npp
subplot(2,2,2)
for i=1:66
    plot(log10(lme_npp(i)/365),FracPD(i),'.','MarkerSize',10,'color',tmap(tid(i,2),:)); hold on;
end
plot(npp(:,1),npp(:,2),'k'); hold on;
plot(npp(:,1),npp(:,2)+2*npp(:,3),'--k'); hold on;
plot(npp(:,1),npp(:,2)-2*npp(:,3),'--k'); hold on;
%title('Color = LME mean temp')
%ylabel('P / (P+D)')
xlabel('log10 NPP')
axis([-0.5 1 0 1])
print('-dpng',[ppath 'lme_scatter_PDgam_ZlDet_NPP_colorT.png'])


