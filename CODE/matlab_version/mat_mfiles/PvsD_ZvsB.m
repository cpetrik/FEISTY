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
% plot info
[ni,nj]=size(lon);

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

geolat_t=lat;
geolon_t=lon;

load([cpath 'cobalt_zoop_biom_means.mat'],'mz_mean_clim','lz_mean_clim','mzloss_mean_clim','lzloss_mean_clim')
load([cpath 'cobalt_det_biom_means.mat'],'det_mean_clim')

frate = 0.3;
tfish = num2str(100+int64(10*frate));
cfile = 'Dc_enc70_cmax-metab20_b175_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
harv = ['All_fish',tfish(2:end)];
tharv = ['Harvest all fish ' num2str(frate)];
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'Means_Climatol_' harv '_' cfile '.mat']);

%% plot info
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
geolon_c = double(lon);
geolat_c = double(lat);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%%
% Convert units
%ESM2.6 already in mg C m-2 or mg C m-2 d-1

z_mean = mz_mean_clim + lz_mean_clim;
z_loss = mzloss_mean_clim+lzloss_mean_clim;

z_mean_grid = z_mean(ID);
z_loss_grid = z_loss(ID);
det_grid = det_mean_clim(ID) * 365;
FracZDet = z_mean_grid ./ (z_mean_grid+det_grid);
FracZB = z_mean_grid ./ (z_mean_grid+b_mean);
FracZlDet = z_loss_grid ./ (z_loss_grid+det_grid);
FracZlB = z_loss_grid ./ (z_loss_grid+b_mean);

RatZDet = z_mean_grid ./ (det_grid);
RatZB = z_mean_grid ./ (b_mean);
RatZlDet = z_loss_grid ./ (det_grid);
RatZlB = z_loss_grid ./ (b_mean);

F = sf_mean+mf_mean;
P = sp_mean+mp_mean+lp_mean;
D = sd_mean+md_mean+ld_mean;
M = mf_mean+mp_mean+md_mean;
L = lp_mean+ld_mean;
FracPD = P ./ (P+D);
FracPF = P ./ (P+F);
FracML = M ./ (M+L);

%% Scatter plot
% COMPARE BY LME?
% COLOR-CODE BY TEMP?
figure(1)
subplot(2,2,1)
plot((FracZDet),FracPD,'.k')
xlabel('Zoop / (Zoop+Det)')
ylabel('P / (P+D)')
ylim([-0.1 1.1])

subplot(2,2,2)
plot((FracZB),FracPD,'.k')
xlabel('Zoop / (Zoop+Bent)')
ylabel('P / (P+D)')
ylim([-0.1 1.1])

subplot(2,2,3)
plot((FracZlDet),FracPD,'.k')
xlabel('ZoopLoss / (ZoopLoss+Det)')
ylabel('P / (P+D)')
ylim([-0.1 1.1])

subplot(2,2,4)
plot((FracZlB),FracPD,'.k')
xlabel('ZoopLoss / (ZoopLoss+Bent)')
ylabel('P / (P+D)')
ylim([-0.1 1.1])
print('-dpng',[ppath cfile '_scatter_frac_PD.png'])


figure(2)
subplot(2,2,1)
plot((FracZDet),FracPF,'.k')
xlabel('Zoop / (Zoop+Det)')
ylabel('P / (P+F)')
ylim([-0.1 1.1])

subplot(2,2,2)
plot((FracZB),FracPF,'.k')
xlabel('Zoop / (Zoop+Bent)')
ylabel('P / (P+F)')
ylim([-0.1 1.1])

subplot(2,2,3)
plot((FracZlDet),FracPF,'.k')
xlabel('ZoopLoss / (ZoopLoss+Det)')
ylabel('P / (P+F)')
ylim([-0.1 1.1])

subplot(2,2,4)
plot((FracZlB),FracPF,'.k')
xlabel('ZoopLoss / (ZoopLoss+Bent)')
ylabel('P / (P+F)')
print('-dpng',[ppath cfile '_scatter_frac_PF.png'])
ylim([-0.1 1.1])


figure(3)
subplot(2,2,1)
plot((FracZDet),FracML,'.k')
xlabel('Zoop / (Zoop+Det)')
ylabel('M / (M+L)')
ylim([-0.1 1.1])

subplot(2,2,2)
plot((FracZB),FracML,'.k')
xlabel('Zoop / (Zoop+Bent)')
ylabel('M / (M+L)')
ylim([-0.1 1.1])

subplot(2,2,3)
plot((FracZlDet),FracML,'.k')
xlabel('ZoopLoss / (ZoopLoss+Det)')
ylabel('M / (M+L)')
ylim([-0.1 1.1])

subplot(2,2,4)
plot((FracZlB),FracML,'.k')
xlabel('ZoopLoss / (ZoopLoss+Bent)')
ylabel('M / (M+L)')
print('-dpng',[ppath cfile '_scatter_frac_ML.png'])
ylim([-0.1 1.1])

%
figure(4)
subplot(2,2,1)
plot(log10(RatZDet),FracPD,'.k')
xlabel('log10 Zoop:Det')
ylabel('P / (P+D)')
ylim([-0.1 1.1])

subplot(2,2,2)
plot(log10(RatZB),FracPD,'.k')
xlabel('log10 Zoop:Bent')
ylabel('P / (P+D)')
ylim([-0.1 1.1])

subplot(2,2,3)
plot(log10(RatZlDet),FracPD,'.k')
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+D)')
ylim([-0.1 1.1])

subplot(2,2,4)
plot(log10(RatZlB),FracPD,'.k')
xlabel('log10 ZoopLoss:Bent')
ylabel('P / (P+D)')
ylim([-0.1 1.1])
print('-dpng',[ppath cfile '_scatter_ratio_PD.png'])


figure(5)
subplot(2,2,1)
plot(log10(RatZDet),FracPF,'.k')
xlabel('log10 Zoop:Det')
ylabel('P / (P+F)')
ylim([-0.1 1.1])

subplot(2,2,2)
plot(log10(RatZB),FracPF,'.k')
xlabel('log10 Zoop:Bent')
ylabel('P / (P+F)')
ylim([-0.1 1.1])

subplot(2,2,3)
plot(log10(RatZlDet),FracPF,'.k')
xlabel('log10 ZoopLoss:Det')
ylabel('P / (P+F)')
ylim([-0.1 1.1])

subplot(2,2,4)
plot(log10(RatZlB),FracPF,'.k')
xlabel('log10 ZoopLoss:Bent')
ylabel('P / (P+F)')
ylim([-0.1 1.1])
print('-dpng',[ppath cfile '_scatter_ratio_PF.png'])


figure(6)
subplot(2,2,1)
plot(log10(RatZDet),FracML,'.k')
xlabel('log10 Zoop:Det')
ylabel('M / (M+L)')
ylim([-0.1 1.1])

subplot(2,2,2)
plot(log10(RatZB),FracML,'.k')
xlabel('log10 Zoop:Bent')
ylabel('M / (M+L)')
ylim([-0.1 1.1])

subplot(2,2,3)
plot(log10(RatZlDet),FracML,'.k')
xlabel('log10 ZoopLoss:Det')
ylabel('M / (M+L)')
ylim([-0.1 1.1])

subplot(2,2,4)
plot(log10(RatZlB),FracML,'.k')
xlabel('log10 ZoopLoss:Bent')
ylabel('M / (M+L)')
ylim([-0.1 1.1])
print('-dpng',[ppath cfile '_scatter_ratio_ML.png'])

%% Pcolor
zees = linspace(min(z_mean_grid),max(z_mean_grid),10);
zls = linspace(min(z_loss_grid),max(z_loss_grid),10);
dees = linspace(min(det_grid),max(det_grid),10);
bees = linspace(min(b_mean),max(b_mean),10);

[zgrid1,dgrid] = meshgrid(zees,dees);
[zgrid2,bgrid] = meshgrid(zees,bees);
[zlgrid1,dlgrid] = meshgrid(zls,dees);
[zlgrid2,blgrid] = meshgrid(zls,bees);

[Nz,zees,zbin] = histcounts(z_mean_grid,zees);
[Nzl,zls,zlbin] = histcounts(z_loss_grid,zls);
[Nd,dees,dbin] = histcounts(det_grid,dees);
[Nb,bees,bbin] = histcounts(b_mean,bees);

%%
ZDfpd = NaN*ones(10);
ZDfpf = NaN*ones(10);
ZDfml = NaN*ones(10);
for z=1:10
    for d=1:10
        zid = find(zbin==z);
        did = find(dbin==d);
        id = intersect(zid,did);
        ZDfpd(d,z) = nanmean(FracPD(id));
        ZDfpf(d,z) = nanmean(FracPF(id));
        ZDfml(d,z) = nanmean(FracML(id));
    end
end

ZBfpd = NaN*ones(10);
ZBfpf = NaN*ones(10);
ZBfml = NaN*ones(10);
for z=1:10
    for b=1:10
        zid = find(zbin==z);
        bid = find(bbin==b);
        id = intersect(zid,bid);
        ZBfpd(b,z) = nanmean(FracPD(id));
        ZBfpf(b,z) = nanmean(FracPF(id));
        ZBfml(b,z) = nanmean(FracML(id));
    end
end

ZlDfpd = NaN*ones(10);
ZlDfpf = NaN*ones(10);
ZlDfml = NaN*ones(10);
for z=1:10
    for d=1:10
        zid = find(zlbin==z);
        did = find(dbin==d);
        id = intersect(zid,did);
        ZlDfpd(d,z) = nanmean(FracPD(id));
        ZlDfpf(d,z) = nanmean(FracPF(id));
        ZlDfml(d,z) = nanmean(FracML(id));
    end
end

ZlBfpd = NaN*ones(10);
ZlBfpf = NaN*ones(10);
ZlBfml = NaN*ones(10);
for z=1:10
    for b=1:10
        zid = find(zlbin==z);
        bid = find(bbin==b);
        id = intersect(zid,bid);
        ZlBfpd(b,z) = nanmean(FracPD(id));
        ZlBfpf(b,z) = nanmean(FracPF(id));
        ZlBfml(b,z) = nanmean(FracML(id));
    end
end

%%
% figure(10)
% pcolor(zgrid1,dgrid,ZDfpd)
% colorbar
% colormap('jet')
% caxis([0 1])
% xlabel('Mean zoop biom')
% ylabel('Mean det biom')
% title('P / (P+D)')
%
% figure(11)
% pcolor(zgrid1,bgrid,ZBfpd)
% colorbar
% colormap('jet')
% caxis([0 1])
% xlabel('Mean zoop biom')
% ylabel('Mean bent biom')
% title('P / (P+D)')
% 
% figure(12)
% pcolor(zlgrid1,dgrid,ZlDfpd)
% colorbar
% colormap('jet')
% caxis([0 1])
% xlabel('Mean zoop loss')
% ylabel('Mean det biom')
% title('P / (P+D)')
% 
% figure(13)
% pcolor(zlgrid1,bgrid,ZlBfpd)
% colorbar
% colormap('jet')
% caxis([0 1])
% xlabel('Mean zoop loss')
% ylabel('Mean bent biom')
% title('P / (P+D)')


figure(14)
subplot(1,2,1)
pcolor(zgrid1,dgrid,ZDfpd)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop biom')
ylabel('Mean det biom')
title('P / (P+D)')

subplot(1,2,2)
pcolor(zgrid1,bgrid,ZBfpd)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop biom')
ylabel('Mean bent biom')
title('P / (P+D)')
print('-dpng',[ppath cfile '_pcolor_zbiom_PD.png'])


figure(15)
subplot(1,2,1)
pcolor(zlgrid1,dgrid,ZlDfpd)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop loss')
ylabel('Mean det biom')
title('P / (P+D)')

subplot(1,2,2)
pcolor(zlgrid1,bgrid,ZlBfpd)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop loss')
ylabel('Mean bent biom')
title('P / (P+D)')
print('-dpng',[ppath cfile '_pcolor_zloss_PD.png'])


figure(16)
subplot(2,2,1)
pcolor(zgrid1,dgrid,ZDfpd)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop biom')
ylabel('Mean det biom')
title('P / (P+D)')

subplot(2,2,2)
pcolor(zgrid1,bgrid,ZBfpd)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop biom')
ylabel('Mean bent biom')
title('P / (P+D)')

subplot(2,2,3)
pcolor(zlgrid1,dgrid,ZlDfpd)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop loss')
ylabel('Mean det biom')
title('P / (P+D)')

subplot(2,2,4)
pcolor(zlgrid1,bgrid,ZlBfpd)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop loss')
ylabel('Mean bent biom')
title('P / (P+D)')
print('-dpng',[ppath cfile '_pcolor_PD.png'])

%%
figure(17)
subplot(2,2,1)
pcolor(zgrid1,dgrid,ZDfpf)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop biom')
ylabel('Mean det biom')
title('P / (P+F)')

subplot(2,2,2)
pcolor(zgrid1,bgrid,ZBfpf)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop biom')
ylabel('Mean bent biom')
title('P / (P+F)')

subplot(2,2,3)
pcolor(zlgrid1,dgrid,ZlDfpf)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop loss')
ylabel('Mean det biom')
title('P / (P+F)')

subplot(2,2,4)
pcolor(zlgrid1,bgrid,ZlBfpf)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop loss')
ylabel('Mean bent biom')
title('P / (P+F)')
print('-dpng',[ppath cfile '_pcolor_PF.png'])


figure(18)
subplot(2,2,1)
pcolor(zgrid1,dgrid,ZDfml)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop biom')
ylabel('Mean det biom')
title('M / (M+L)')

subplot(2,2,2)
pcolor(zgrid1,bgrid,ZBfml)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop biom')
ylabel('Mean bent biom')
title('M / (M+L)')

subplot(2,2,3)
pcolor(zlgrid1,dgrid,ZlDfml)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop loss')
ylabel('Mean det biom')
title('M / (M+L)')

subplot(2,2,4)
pcolor(zlgrid1,bgrid,ZlBfml)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop loss')
ylabel('Mean bent biom')
title('M / (M+L)')
print('-dpng',[ppath cfile '_pcolor_ML.png'])

%% CONTOURS
%P:D
figure(19)
subplot(2,2,1)
contourf(zgrid1,dgrid,ZDfpd)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop biom')
ylabel('Mean det biom')
title('P / (P+D)')

subplot(2,2,2)
contourf(zgrid1,bgrid,ZBfpd)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop biom')
ylabel('Mean bent biom')
title('P / (P+D)')

subplot(2,2,3)
contourf(zlgrid1,dgrid,ZlDfpd)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop loss')
ylabel('Mean det biom')
title('P / (P+D)')

subplot(2,2,4)
contourf(zlgrid1,bgrid,ZlBfpd)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop loss')
ylabel('Mean bent biom')
title('P / (P+D)')
print('-dpng',[ppath cfile '_contourf_PD.png'])

% P:F
figure(20)
subplot(2,2,1)
contourf(zgrid1,dgrid,ZDfpf)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop biom')
ylabel('Mean det biom')
title('P / (P+F)')

subplot(2,2,2)
contourf(zgrid1,bgrid,ZBfpf)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop biom')
ylabel('Mean bent biom')
title('P / (P+F)')

subplot(2,2,3)
contourf(zlgrid1,dgrid,ZlDfpf)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop loss')
ylabel('Mean det biom')
title('P / (P+F)')

subplot(2,2,4)
contourf(zlgrid1,bgrid,ZlBfpf)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop loss')
ylabel('Mean bent biom')
title('P / (P+F)')
print('-dpng',[ppath cfile '_contourf_PF.png'])


figure(21)
subplot(2,2,1)
contourf(zgrid1,dgrid,ZDfml)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop biom')
ylabel('Mean det biom')
title('M / (M+L)')

subplot(2,2,2)
contourf(zgrid1,bgrid,ZBfml)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop biom')
ylabel('Mean bent biom')
title('M / (M+L)')

subplot(2,2,3)
contourf(zlgrid1,dgrid,ZlDfml)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop loss')
ylabel('Mean det biom')
title('M / (M+L)')

subplot(2,2,4)
contourf(zlgrid1,bgrid,ZlBfml)
colorbar
cmocean('balance')
caxis([0 1])
xlabel('Mean zoop loss')
ylabel('Mean bent biom')
title('M / (M+L)')
print('-dpng',[ppath cfile '_contourf_ML.png'])



