%POEM catch vs. SAUP catch by LME
%Use same methods as Stock et al. 2017 to reduce SAUP dataset

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'LME_clim_temp.mat']);

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')

AREA_OCN = max(area,1);

% POEM file info
frate = 0.3;
tfish = num2str(100+int64(10*frate));
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

%%
bees = 0.05:0.025:0.1;
cmaxs = 20:5:30;
r = NaN*ones(length(bees),length(cmaxs),4);
rmse = NaN*ones(length(bees),length(cmaxs),4);
F = NaN*ones(length(bees),length(cmaxs),4);
for e=1:length(bees)
    for c=1:length(cmaxs)
        bent_eff = bees(e);
        h = cmaxs(c);
        tcfn = num2str(h);
        tbe = num2str(100+int64(100*bent_eff));
        
        cfile = ['Dc_enc70-b200_cm',tcfn,...
            '_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE',...
            tbe(2:end),'_noCC_RE00100'];
        ppath = [pp cfile '/'];
        dpath = [dp cfile '/'];
        
        load([dpath 'LME_SAUP_stats_' cfile '.mat'],'fish_stat')
        r(e,c,:)    = fish_stat(:,1);
        rmse(e,c,:) = fish_stat(:,2);
        F(e,c,:)    = fish_stat(:,3);
        
    end
end


%% Plots
nc = length(cmaxs);
ne = length(bees);
encs2 = [bees 0.125];
cmaxs2 = [cmaxs 35];
[cgrid,egrid]=meshgrid(cmaxs2,encs2);
r2  = NaN*ones(ne+1,nc+1,4);
rmse2  = NaN*ones(ne+1,nc+1,4);
F2  = NaN*ones(ne+1,nc+1,4);
r2(1:ne,1:nc,:) = r;
rmse2(1:ne,1:nc,:) = rmse;
F2(1:ne,1:nc,:) = F;

cmap_ther = colormap(cmocean('thermal')); close all;
cmap_revt = flipud(cmap_ther);

cfile2 = ['Dc_enc70-b200_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_noCC_RE00100'];

% r
figure(1)
subplot(2,2,1)
pcolor(egrid,cgrid,squeeze(r2(:,:,1)))
cmocean('thermal')
colorbar
%caxis([0 1])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('corr All')

subplot(2,2,2)
pcolor(egrid,cgrid,squeeze(r2(:,:,2)))
cmocean('thermal')
colorbar
%caxis([0 1])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('corr F')

subplot(2,2,3)
pcolor(egrid,cgrid,squeeze(r2(:,:,3)))
cmocean('thermal')
colorbar
%caxis([0 1])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('corr P')

subplot(2,2,4)
pcolor(egrid,cgrid,squeeze(r2(:,:,4)))
cmocean('thermal')
colorbar
%caxis([0 1])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('corr D')
print('-dpng',[pp cfile2 '_CmaxBE_SAUP_r.png'])

%% rmse
figure(2)
subplot(2,2,1)
pcolor(egrid,cgrid,squeeze(rmse2(:,:,1)))
colormap(cmap_revt)
colorbar
%caxis([0 1])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('RMSE All')

subplot(2,2,2)
pcolor(egrid,cgrid,squeeze(rmse2(:,:,2)))
colormap(cmap_revt)
colorbar
%caxis([0 1])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('RMSE F')

subplot(2,2,3)
pcolor(egrid,cgrid,squeeze(rmse2(:,:,3)))
colormap(cmap_revt)
colorbar
%caxis([0 1])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('RMSE P')

subplot(2,2,4)
pcolor(egrid,cgrid,squeeze(rmse2(:,:,4)))
colormap(cmap_revt)
colorbar
%caxis([0 1])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('RMSE D')
print('-dpng',[pp cfile2 '_CmaxBE_SAUP_RMSE.png'])

%% Fmed
figure(3)
subplot(2,2,1)
pcolor(egrid,cgrid,squeeze(F2(:,:,1)))
cmocean('balance')
colorbar
caxis([0 2])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('Fmed All')

subplot(2,2,2)
pcolor(egrid,cgrid,squeeze(F2(:,:,2)))
cmocean('balance')
colorbar
caxis([0 2])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('Fmed F')

subplot(2,2,3)
pcolor(egrid,cgrid,squeeze(F2(:,:,3)))
cmocean('balance')
colorbar
caxis([0 2])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('Fmed P')

subplot(2,2,4)
pcolor(egrid,cgrid,squeeze(F2(:,:,4)))
cmocean('balance')
colorbar
caxis([0 2])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('Fmed D')
print('-dpng',[pp cfile2 '_CmaxBE_SAUP_Fmed.png'])

%% Combined
corr_all  = NaN*ones(ne+1,nc+1);
rmse_all  = NaN*ones(ne+1,nc+1);
F_all  = NaN*ones(ne+1,nc+1);
corr_PD  = NaN*ones(ne+1,nc+1);
rmse_PD  = NaN*ones(ne+1,nc+1);
F_PD  = NaN*ones(ne+1,nc+1);

corr_all(1:ne,1:nc) = normalize(r(:,:,1))+normalize(r(:,:,2))+normalize(r(:,:,3))...
    +normalize(r(:,:,4));
max(corr_all(:)) %e=90,c=30

rmse_all(1:ne,1:nc) = normalize(rmse(:,:,1))+normalize(rmse(:,:,2))+...
    normalize(rmse(:,:,3))+normalize(rmse(:,:,4));
min(rmse_all(:)) %e=100,c=30

F_all(1:ne,1:nc) = normalize(F(:,:,1))+normalize(F(:,:,2))+normalize(F(:,:,3))...
    +normalize(F(:,:,4));
min(abs(1-F_all(:))) %e=10,c=20
aF = abs(1-F_all);

corr_PD(1:ne,1:nc) = normalize(r(:,:,3))+normalize(r(:,:,4));
max(corr_PD(:)) %e=60,c=20

rmse_PD(1:ne,1:nc) = normalize(rmse(:,:,3))+normalize(rmse(:,:,4));
min(rmse_PD(:)) %e=20,c=10

F_PD(1:ne,1:nc) = normalize(F(:,:,3))+normalize(F(:,:,4));
min(abs(1-F_PD(:))) %e=20,c=30
aFPD = abs(1-F_PD);

%% Figs
figure(4)
subplot(2,3,1)
pcolor(egrid,cgrid,corr_all)
cmocean('thermal')
colorbar
%caxis([0 2])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('Corr all combined')

subplot(2,3,4)
pcolor(egrid,cgrid,corr_PD)
cmocean('thermal')
colorbar
%caxis([0 2])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('Corr P&D')

subplot(2,3,2)
pcolor(egrid,cgrid,rmse_all)
colormap(cmap_revt)
colorbar
%caxis([0 2])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('RMSE all combined')

subplot(2,3,5)
pcolor(egrid,cgrid,rmse_PD)
colormap(cmap_revt)
colorbar
%caxis([0 2])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('RMSE P&D')

subplot(2,3,3)
pcolor(egrid,cgrid,F_all)
cmocean('balance')
colorbar
%caxis([0 2])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('Fmed all combined')

subplot(2,3,6)
pcolor(egrid,cgrid,F_PD)
cmocean('balance')
colorbar
%caxis([0 2])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('Fmed P&D')
print('-dpng',[pp cfile2 '_CmaxBE_SAUP_combos.png'])

%% P&D combine metrics
nr = normalize(corr_PD);
nrmse = normalize(max(rmse_PD(:))-rmse_PD);
nF = normalize(abs(1-F_PD));

Pall = NaN*ones(ne+1,nc+1);
Dall = NaN*ones(ne+1,nc+1);
p1=r(:,:,3);
p2=rmse(:,:,3);
p3=F(:,:,3);
d1=r(:,:,4);
d2=rmse(:,:,4);
d3=F(:,:,4);

mPD = nr+nrmse+nF;
%P
Pall(1:ne,1:nc) = normalize(p1) + normalize(max(p2(:))-p2) + normalize(abs(1-p3));
%D
Dall(1:ne,1:nc) = normalize(d1) + normalize(max(d2(:))-d2) + normalize(abs(1-d3));

figure
subplot(1,3,2)
pcolor(egrid,cgrid,mPD)
cmocean('thermal')
colorbar
%caxis([0 2])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('All metrics P&D')

subplot(1,3,1)
pcolor(egrid,cgrid,Pall)
cmocean('thermal')
colorbar
%caxis([0 2])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('All metrics P')

subplot(1,3,3)
pcolor(egrid,cgrid,Dall)
cmocean('thermal')
colorbar
%caxis([0 2])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('All metrics D')
print('-dpng',[pp cfile2 '_CmaxBE_SAUP_comboPD.png'])

cp = '/Volumes/GFDL/CSV/Matlab_new_size/';
save([cp 'Bio_rates/' cfile2 '_clim_',harv,'_SAUPcomp_CmaxBE.mat']);


