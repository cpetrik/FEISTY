% Plot mean biomasses of zoop to compare to POEM results

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

%% grid info
[ni,nj]=size(lon);

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

AREA_OCN = max(area,1);

geolat_t=lat;
geolon_t=lon;

%% Assign a color to each LME based on temp
%tmap=colormap(jet(66));
tmap=cmocean('thermal',66);
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);
close all

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
%load([fpath 'Means_Climatol_' harv '_' cfile '.mat']);
load([fpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);
load([fpath 'LME_clim_fished_' harv '_' cfile '.mat']);
lme_area_km2 = lme_area * 1e-6;

%%
% now units in g m-2 or g m-2 d-1

lme_zl = lme_zl*365;
lme_det = lme_det*365;
lme_npp = lme_npp*365;

FracZDet = lme_z ./ (lme_z+lme_det);
FracZB = lme_z ./ (lme_z+lme_mbio(:,9));
FracZlDet = lme_zl ./ (lme_zl+lme_det);
FracZlB = lme_zl ./ (lme_zl+lme_mbio(:,9));

RatZDet = lme_z ./ (lme_det);
RatZB = lme_z ./ (lme_mbio(:,9));
RatZlDet = lme_zl ./ (lme_det);
RatZlB = lme_zl ./ (lme_mbio(:,9));

F = lme_mbio(:,1)+lme_mbio(:,4);
P = lme_mbio(:,2)+lme_mbio(:,5)+lme_mbio(:,7);
D = lme_mbio(:,3)+lme_mbio(:,6)+lme_mbio(:,8);
M = lme_mbio(:,4)+lme_mbio(:,5)+lme_mbio(:,6);
L = lme_mbio(:,7)+lme_mbio(:,8);
FracPD = P ./ (P+D);
FracPF = P ./ (P+F);
FracLM = L ./ (L+M);

lme=1:66;

T=table(lme',lme_ptemp(:,1),RatZDet,RatZB,RatZlDet,RatZlB,...
    FracPD,FracPF,FracLM,lme_depth,lme_shal_frac,lme_npp,'VariableNames',...
    {'LME','LME_ptemp','RatZDet','RatZB','RatZlDet','RatZlB','FracPD',...
    'FracPF','FracLM','LME_depth','LME_Frac200','NPP'});
writetable(T,[fpath 'LME_ZBratios_clim_fished_',harv,'.csv'],'Delimiter',',');
save([fpath 'LME_ZBratios_clim_fished_',harv,'_' cfile '.mat'],...
    'RatZDet','RatZB','RatZlDet','RatZlB','FracPD','FracPF','FracLM',...
    'lme_ptemp','lme_depth','lme_shal_frac','lme_npp');

%% Scatter plot
% COMPARE BY LME?
% COLOR-CODE BY TEMP?
figure(1)
subplot(2,2,1)
for i=1:66
    plot(FracZDet(i),FracPD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('Zoop / (Zoop+Det)')
ylabel('P / (P+D)')
ylim([-0.1 1.1])
title('Color = LME mean temp')

subplot(2,2,2)
for i=1:66
    plot(FracZB(i),FracPD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('Zoop / (Zoop+Bent)')
ylabel('P / (P+D)')
ylim([-0.1 1.1])

subplot(2,2,3)
for i=1:66
    plot(FracZlDet(i),FracPD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('ZoopLoss / (ZoopLoss+Det)')
ylabel('P / (P+D)')
ylim([-0.1 1.1])

subplot(2,2,4)
for i=1:66
    plot(FracZlB(i),FracPD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('ZoopLoss / (ZoopLoss+Bent)')
ylabel('P / (P+D)')
ylim([-0.1 1.1])
%print('-dpng',[ppath 'lme_scatter_frac_PD.png'])

%%
figure(2)
subplot(2,2,1)
for i=1:66
    plot(FracZDet(i),FracPF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('Zoop / (Zoop+Det)')
ylabel('P / (P+F)')
ylim([-0.1 1.1])
title('Color = LME mean temp')

subplot(2,2,2)
for i=1:66
    plot(FracZB(i),FracPF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('Zoop / (Zoop+Bent)')
ylabel('P / (P+F)')
ylim([-0.1 1.1])

subplot(2,2,3)
for i=1:66
    plot(FracZlDet(i),FracPF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('ZoopLoss / (ZoopLoss+Det)')
ylabel('P / (P+F)')
ylim([-0.1 1.1])

subplot(2,2,4)
for i=1:66
    plot(FracZlB(i),FracPF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('ZoopLoss / (ZoopLoss+Bent)')
ylabel('P / (P+F)')
%print('-dpng',[ppath 'lme_scatter_frac_PF.png'])
ylim([-0.1 1.1])


figure(3)
subplot(2,2,1)
for i=1:66
    plot(FracZDet(i),FracLM(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('Zoop / (Zoop+Det)')
ylabel('L / (L+M)')
ylim([-0.1 1.1])
title('Color = LME mean temp')

subplot(2,2,2)
for i=1:66
    plot(FracZB(i),FracLM(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('Zoop / (Zoop+Bent)')
ylabel('L / (L+M)')
ylim([-0.1 1.1])

subplot(2,2,3)
for i=1:66
    plot(FracZlDet(i),FracLM(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('ZoopLoss / (ZoopLoss+Det)')
ylabel('L / (L+M)')
ylim([-0.1 1.1])

subplot(2,2,4)
for i=1:66
    plot(FracZlB(i),FracLM(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('ZoopLoss / (ZoopLoss+Bent)')
ylabel('L / (L+M)')
%print('-dpng',[ppath 'lme_scatter_frac_LM.png'])
ylim([-0.1 1.1])

%
figure(4)
subplot(2,2,1)
for i=1:66
    plot(log10(RatZDet(i)),FracPD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 Zoop:Det')
ylabel('P / (P+D)')
ylim([-0.1 1.1])
title('Color = LME mean temp')

subplot(2,2,2)
for i=1:66
    plot(log10(RatZB(i)),FracPD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 Zoop:Bent')
ylabel('P / (P+D)')
ylim([-0.1 1.1])

subplot(2,2,3)
for i=1:66
    plot(log10(RatZlDet(i)),FracPD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
ylim([-0.1 1.1])

subplot(2,2,4)
for i=1:66
    plot(log10(RatZlB(i)),FracPD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 ZoopLoss:Bent')
ylabel('P / (P+D)')
ylim([-0.1 1.1])
%print('-dpng',[ppath 'lme_scatter_ratio_PD.png'])


figure(5)
subplot(2,2,1)
for i=1:66
    plot(log10(RatZDet(i)),FracPF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 Zoop:Det')
ylabel('P / (P+F)')
ylim([-0.1 1.1])
title('Color = LME mean temp')

subplot(2,2,2)
for i=1:66
    plot(log10(RatZB(i)),FracPF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 Zoop:Bent')
ylabel('P / (P+F)')
ylim([-0.1 1.1])

subplot(2,2,3)
for i=1:66
    plot(log10(RatZlDet(i)),FracPF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
ylim([-0.1 1.1])

subplot(2,2,4)
for i=1:66
    plot(log10(RatZlB(i)),FracPF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 ZoopLoss:Bent')
ylabel('P / (P+F)')
ylim([-0.1 1.1])
%print('-dpng',[ppath 'lme_scatter_ratio_PF.png'])


figure(6)
subplot(2,2,1)
for i=1:66
    plot(log10(RatZDet(i)),FracLM(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 Zoop:Det')
ylabel('L / (L+M)')
ylim([-0.1 1.1])
title('Color = LME mean temp')

subplot(2,2,2)
for i=1:66
    plot(log10(RatZB(i)),FracLM(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 Zoop:Bent')
ylabel('L / (L+M)')
ylim([-0.1 1.1])

subplot(2,2,3)
for i=1:66
    plot(log10(RatZlDet(i)),FracLM(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
ylim([-0.1 1.1])

subplot(2,2,4)
for i=1:66
    plot(log10(RatZlB(i)),FracLM(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 ZoopLoss:Bent')
ylabel('L / (L+M)')
ylim([-0.1 1.1])
%print('-dpng',[ppath 'lme_scatter_ratio_LM.png'])

%% Best ones
figure(7)
subplot(3,2,1)
for i=1:66
    plot(FracZB(i),FracPD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('Zoop / (Zoop+Bent)')
ylabel('P / (P+D)')
ylim([-0.1 1.1])
title('Color = LME mean temp')
subplot(3,2,2)
for i=1:66
    plot(FracZlB(i),FracPD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('ZoopLoss / (ZoopLoss+Bent)')
ylabel('P / (P+D)')
ylim([-0.1 1.1])

subplot(3,2,3)
for i=1:66
    plot(FracZB(i),FracPF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('Zoop / (Zoop+Bent)')
ylabel('P / (P+F)')
ylim([-0.1 1.1])
subplot(3,2,4)
for i=1:66
    plot(FracZlB(i),FracPF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('ZoopLoss / (ZoopLoss+Bent)')
ylabel('P / (P+F)')
ylim([-0.1 1.1])

subplot(3,2,5)
for i=1:66
    plot(FracZB(i),FracLM(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('Zoop / (Zoop+Bent)')
ylabel('L / (L+M)')
ylim([-0.1 1.1])
subplot(3,2,6)
for i=1:66
    plot(FracZlB(i),FracLM(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('ZoopLoss / (ZoopLoss+Bent)')
ylabel('L / (L+M)')
ylim([-0.1 1.1])
%print('-dpng',[ppath 'lme_scatter_frac_bent.png'])


%%
figure(8)
subplot(3,2,1)
for i=1:66
    plot(log10(RatZB(i)),FracPD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 Zoop:Bent')
ylabel('P / (P+D)')
ylim([-0.1 1.1])
title('Color = LME mean temp')
subplot(3,2,2)
for i=1:66
    plot(log10(RatZlB(i)),FracPD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 ZoopLoss:Bent')
ylabel('P / (P+D)')
ylim([-0.1 1.1])

subplot(3,2,3)
for i=1:66
    plot(log10(RatZB(i)),FracPF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 Zoop:Bent')
ylabel('P / (P+F)')
ylim([-0.1 1.1])
subplot(3,2,4)
for i=1:66
    plot(log10(RatZlB(i)),FracPF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 ZoopLoss:Bent')
ylabel('P / (P+F)')
ylim([-0.1 1.1])

subplot(3,2,5)
for i=1:66
    plot(log10(RatZB(i)),FracLM(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 Zoop:Bent')
ylabel('L / (L+M)')
ylim([-0.1 1.1])
subplot(3,2,6)
for i=1:66
    plot(log10(RatZlB(i)),FracLM(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 ZoopLoss:Bent')
ylabel('L / (L+M)')
ylim([-0.1 1.1])
%print('-dpng',[ppath 'lme_scatter_ratio_bent.png'])

%%
figure(9)
subplot(3,2,1)
for i=1:66
    plot(log10(RatZDet(i)),FracPD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 Zoop:Det')
ylabel('P / (P+D)')
ylim([-0.1 1.1])
title('Color = LME mean temp')
subplot(3,2,2)
for i=1:66
    plot(log10(RatZlDet(i)),FracPD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
ylim([-0.1 1.1])

subplot(3,2,3)
for i=1:66
    plot(log10(RatZDet(i)),FracPF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 Zoop:Det')
ylabel('P / (P+F)')
ylim([-0.1 1.1])
subplot(3,2,4)
for i=1:66
    plot(log10(RatZlDet(i)),FracPF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
ylim([-0.1 1.1])

subplot(3,2,5)
for i=1:66
    plot(log10(RatZDet(i)),FracLM(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 Zoop:Det')
ylabel('L / (L+M)')
ylim([-0.1 1.1])
subplot(3,2,6)
for i=1:66
    plot(log10(RatZlDet(i)),FracLM(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
xlabel('log10 ZoopLoss:Det')
ylabel('L / (L+M)')
ylim([-0.1 1.1])
print('-dpng',[ppath 'lme_scatter_ratio_det.png'])



