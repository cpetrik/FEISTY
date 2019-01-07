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
r = NaN*ones(length(bees),length(cmaxs),2);
rmse = NaN*ones(length(bees),length(cmaxs),2);
F = NaN*ones(length(bees),length(cmaxs),2);
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
        
        load([dpath 'LME_DvD_stats_' cfile '.mat'],'fish_stat')
        r(e,c,:)    = fish_stat(1,:);
        rmse(e,c,:) = fish_stat(2,:);
        F(e,c,:)    = fish_stat(3,:);
        
    end
end


%% Plots
nc = length(cmaxs);
ne = length(bees);
encs2 = [bees 0.125];
cmaxs2 = [cmaxs 35];
[cgrid,egrid]=meshgrid(cmaxs2,encs2);
r2  = NaN*ones(ne+1,nc+1,2);
rmse2  = NaN*ones(ne+1,nc+1,2);
F2  = NaN*ones(ne+1,nc+1,2);
r2(1:ne,1:nc,:) = r;
rmse2(1:ne,1:nc,:) = rmse;
F2(1:ne,1:nc,:) = F;

cmap_ther = colormap(cmocean('thermal')); close all;
cmap_revt = flipud(cmap_ther);

cfile2 = ['Dc_enc70-b200_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_noCC_RE00100'];

%% r
figure(1)
subplot(3,2,1)
pcolor(egrid,cgrid,squeeze(r2(:,:,1)))
cmocean('thermal')
colorbar
caxis([0.22 0.34])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('corr All')

subplot(3,2,2)
pcolor(egrid,cgrid,squeeze(r2(:,:,2)))
cmocean('thermal')
colorbar
caxis([0.22 0.34])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('corr noLELC')

% rmse
subplot(3,2,3)
pcolor(egrid,cgrid,squeeze(rmse2(:,:,1)))
colormap(cmap_revt)
colorbar
caxis([0.28 0.36])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('RMSE All')

subplot(3,2,4)
pcolor(egrid,cgrid,squeeze(rmse2(:,:,2)))
colormap(cmap_revt)
colorbar
caxis([0.28 0.36])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('RMSE noLELC')

% Fmed
subplot(3,2,5)
pcolor(egrid,cgrid,squeeze(F2(:,:,1)))
cmocean('balance')
colorbar
caxis([0.5 1.5])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('Fmed All')

subplot(3,2,6)
pcolor(egrid,cgrid,squeeze(F2(:,:,2)))
cmocean('balance')
colorbar
caxis([0.5 1.5])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('Fmed noLELC')
stamp([harv '_' cfile])
print('-dpng',[pp cfile2 '_CmaxBE_DvD_stats.png'])

%% Combined
Pall = NaN*ones(ne+1,nc+1);
Dall = NaN*ones(ne+1,nc+1);
p1=r(:,:,1);
p2=rmse(:,:,1);
p3=F(:,:,1);
d1=r(:,:,2);
d2=rmse(:,:,2);
d3=F(:,:,2);

%All LMEs
Pall(1:ne,1:nc) = normalize(p1) + normalize(max(p2(:))-p2) + normalize(abs(1-p3));
%No LELC
Dall(1:ne,1:nc) = normalize(d1) + normalize(max(d2(:))-d2) + normalize(abs(1-d3));

figure
subplot(1,2,1)
pcolor(egrid,cgrid,Pall)
cmocean('thermal')
colorbar
caxis([0.6 2.2])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('All metrics all LMEs')

subplot(1,2,2)
pcolor(egrid,cgrid,Dall)
cmocean('thermal')
colorbar
caxis([0.6 2.2])
set(gca,'XTick',bees,'XTickLabel',bees,...
    'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('BE')
ylabel('Cmax coeff')
title('All metrics no LELC')
print('-dpng',[pp cfile2 '_CmaxBE_DvD_combo.png'])

cp = '/Volumes/GFDL/CSV/Matlab_new_size/';
save([cp 'Bio_rates/' cfile2 '_clim_',harv,'_DvDcomp_CmaxBE.mat']);


