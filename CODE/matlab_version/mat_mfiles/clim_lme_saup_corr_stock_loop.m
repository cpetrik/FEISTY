%POEM catch vs. SAUP catch by LME
%Use same methods as Stock et al. 2017 to reduce SAUP dataset

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/CSV/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'LME_clim_temp.mat']);

%use weighted catches
load([spath 'SAUP_LME_Catch_annual.mat'],'yr','totcatch','lme_catch',...
    'Flme_wcatch','Dlme_wcatch','Plme_wcatch');

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')

AREA_OCN = max(area,1);

% "lfile" never changes, has lme areas
lfile = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100';
lpath = ['/Volumes/GFDL/NC/Matlab_new_size/' lfile '/'];
load([lpath 'LME_clim_fished03_' lfile '.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;

% POEM file info
frate = 0.3;
tfish = num2str(100+int64(10*frate));
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';
cfile2 = ['Dc_enc-b200_cm_m-b175-k09_fcrit20',...
            '_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100'];


%% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);
close all

load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

x=-6:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% SAUP
Flme_catch_all = nansum(Flme_wcatch,3);
Plme_catch_all = nansum(Plme_wcatch,3);
Dlme_catch_all = nansum(Dlme_wcatch,3);

%1950-2006 SAUP average
id = find(yr>1950 & yr<=2006);

slme_mcatch = nanmean(lme_catch(id,:));
slme_mcatch = slme_mcatch';
Fslme_mcatch = nanmean(Flme_catch_all(id,:));
Fslme_mcatch = Fslme_mcatch';
Pslme_mcatch = nanmean(Plme_catch_all(id,:));
Pslme_mcatch = Pslme_mcatch';
Dslme_mcatch = nanmean(Dlme_catch_all(id,:));
Dslme_mcatch = Dslme_mcatch';

slme_mcatch10 = NaN*ones(size(slme_mcatch));
Flme_mcatch10 = NaN*ones(size(slme_mcatch));
Plme_mcatch10 = NaN*ones(size(slme_mcatch));
Dlme_mcatch10 = NaN*ones(size(slme_mcatch));

%Top 10 yrs by LME SAUP 
for i=1:66
    [sort_lme_catch,ix] = sort(lme_catch(:,i),'descend');
    sort_Flme_catch = Flme_catch_all(ix,i);
    sort_Plme_catch = Plme_catch_all(ix,i);
    sort_Dlme_catch = Dlme_catch_all(ix,i);
    slme_mcatch10(i) = nanmean(sort_lme_catch(1:10));
    Flme_mcatch10(i) = nanmean(sort_Flme_catch(1:10));
    Plme_mcatch10(i) = nanmean(sort_Plme_catch(1:10));
    Dlme_mcatch10(i) = nanmean(sort_Dlme_catch(1:10));
end

% MT/km2
slme_mcatch10 = slme_mcatch10 ./ lme_area_km2;
Flme_mcatch10 = Flme_mcatch10 ./ lme_area_km2;
Plme_mcatch10 = Plme_mcatch10 ./ lme_area_km2;
Dlme_mcatch10 = Dlme_mcatch10 ./ lme_area_km2;

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

%% Loop over params
encs = 10:10:100;
cmaxs = 10:10:50;
frate = 0.3;

%%
r = NaN*ones(length(encs),length(cmaxs),4);
rmse = NaN*ones(length(encs),length(cmaxs),4);
F = NaN*ones(length(encs),length(cmaxs),4);
for c=1:length(cmaxs)
    for e=1:length(encs)
        gam = encs(e);
        h = cmaxs(c);
        tcfn = num2str(h);
        tefn = num2str(round(gam));
        
        cfile = ['Dc_enc',tefn,'-b200_cm',tcfn,'_m-b175-k09_fcrit20',...
            '_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100'];
        fpath=['/Volumes/GFDL/CSV/Matlab_new_size/' cfile '/'];
        
        ppath = [pp cfile '/'];
        dpath = [dp cfile '/'];
        
        load([dpath 'LME_clim_',harv,'_loop.mat']);
        
        %% POEM LME biomass in MT
        plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
        plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
        plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
        plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
        % MT/km2
        plme_mcatch = plme_mcatch ./ lme_area_km2;
        plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;
        plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
        plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;
        
        l10p=log10(plme_mcatch);
        l10pF=log10(plme_Fmcatch);
        l10pP=log10(plme_Pmcatch);
        l10pD=log10(plme_Dmcatch);
        
        %% Stats
        % Drop Arctic, Antarctic, Hawaii, Australia -------------------------
        %r
        r(e,c,1)=corr(l10s(keep),l10p(keep));
        r(e,c,2)=corr(l10sF(keep),l10pF(keep));
        r(e,c,3)=corr(l10sP(keep),l10pP(keep));
        r(e,c,4)=corr(l10sD(keep),l10pD(keep));
        
        %root mean square error
        o=l10s(keep);
        p=l10p(keep);
        n = length(o);
        num=nansum((p-o).^2);
        rmse(e,c,1) = sqrt(num/n);
        
        o=l10sF(keep);
        p=l10pF(keep);
        n = length(o);
        num=nansum((p-o).^2);
        rmse(e,c,2) = sqrt(num/n);
        
        o=l10sP(keep);
        p=l10pP(keep);
        n = length(o);
        num=nansum((p-o).^2);
        rmse(e,c,3) = sqrt(num/n);
        
        o=l10sD(keep);
        p=l10pD(keep);
        n = length(o);
        num=nansum((p-o).^2);
        rmse(e,c,4) = sqrt(num/n);
        
        %Fmed
        F(e,c,1)=10^(median(l10s(keep)-l10p(keep)));
        F(e,c,2)=10^(median(l10sF(keep)-l10pF(keep)));
        F(e,c,3)=10^(median(l10sP(keep)-l10pP(keep)));
        F(e,c,4)=10^(median(l10sD(keep)-l10pD(keep)));

    end
end


%% Plots
nc = length(cmaxs);
ne = length(encs);
encs2 = [encs 110];
cmaxs2 = [cmaxs 60];
[cgrid,egrid]=meshgrid(cmaxs2,encs2);
r2  = NaN*ones(ne+1,nc+1,4);
rmse2  = NaN*ones(ne+1,nc+1,4);
F2  = NaN*ones(ne+1,nc+1,4);
r2(1:ne,1:nc,:) = r;
rmse2(1:ne,1:nc,:) = rmse;
F2(1:ne,1:nc,:) = F;

cmap_ther = colormap(cmocean('thermal')); close all;
cmap_revt = flipud(cmap_ther);

% r
figure(1)
subplot(2,2,1)
pcolor(egrid,cgrid,squeeze(r2(:,:,1)))
cmocean('thermal')
colorbar
%caxis([0 1])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('corr All')

subplot(2,2,2)
pcolor(egrid,cgrid,squeeze(r2(:,:,2)))
cmocean('thermal')
colorbar
%caxis([0 1])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('corr F')

subplot(2,2,3)
pcolor(egrid,cgrid,squeeze(r2(:,:,3)))
cmocean('thermal')
colorbar
%caxis([0 1])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('corr P')

subplot(2,2,4)
pcolor(egrid,cgrid,squeeze(r2(:,:,4)))
cmocean('thermal')
colorbar
%caxis([0 1])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('corr D')
print('-dpng',[pp cfile2 '_EncCmax_SAUP_r.png'])

%% rmse
figure(2)
subplot(2,2,1)
pcolor(egrid,cgrid,squeeze(rmse2(:,:,1)))
colormap(cmap_revt)
colorbar
%caxis([0 1])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('RMSE All')

subplot(2,2,2)
pcolor(egrid,cgrid,squeeze(rmse2(:,:,2)))
colormap(cmap_revt)
colorbar
%caxis([0 1])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('RMSE F')

subplot(2,2,3)
pcolor(egrid,cgrid,squeeze(rmse2(:,:,3)))
colormap(cmap_revt)
colorbar
%caxis([0 1])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('RMSE P')

subplot(2,2,4)
pcolor(egrid,cgrid,squeeze(rmse2(:,:,4)))
colormap(cmap_revt)
colorbar
%caxis([0 1])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('RMSE D')
print('-dpng',[pp cfile2 '_EncCmax_SAUP_RMSE.png'])

%% Fmed
figure(3)
subplot(2,2,1)
pcolor(egrid,cgrid,squeeze(F2(:,:,1)))
cmocean('balance')
colorbar
caxis([0 2])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('Fmed All')

subplot(2,2,2)
pcolor(egrid,cgrid,squeeze(F2(:,:,2)))
cmocean('balance')
colorbar
caxis([0 2])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('Fmed F')

subplot(2,2,3)
pcolor(egrid,cgrid,squeeze(F2(:,:,3)))
cmocean('balance')
colorbar
caxis([0 2])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('Fmed P')

subplot(2,2,4)
pcolor(egrid,cgrid,squeeze(F2(:,:,4)))
cmocean('balance')
colorbar
caxis([0 2])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('Fmed D')
print('-dpng',[pp cfile2 '_EncCmax_SAUP_Fmed.png'])

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

rmse_PD(1:ne,1:nc) = normalize(F(:,:,3))+normalize(F(:,:,4));
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
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('Corr all combined')

subplot(2,3,4)
pcolor(egrid,cgrid,corr_PD)
cmocean('thermal')
colorbar
%caxis([0 2])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('Corr P&D')

subplot(2,3,2)
pcolor(egrid,cgrid,rmse_all)
colormap(cmap_revt)
colorbar
%caxis([0 2])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('RMSE all combined')

subplot(2,3,5)
pcolor(egrid,cgrid,rmse_PD)
colormap(cmap_revt)
colorbar
%caxis([0 2])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('RMSE P&D')

subplot(2,3,3)
pcolor(egrid,cgrid,F_all)
cmocean('balance')
colorbar
%caxis([0 2])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('Fmed all combined')

subplot(2,3,6)
pcolor(egrid,cgrid,F_PD)
cmocean('balance')
colorbar
%caxis([0 2])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('Fmed P&D')
print('-dpng',[pp cfile2 '_EncCmax_SAUP_combos.png'])

%% P&D combine metrics
nr = normalize(corr_PD);
nrmse = normalize(max(rmse_PD(:))-rmse_PD);
nF = normalize(abs(1-F_PD));

mPD = nr+nrmse+nF;

figure
pcolor(egrid,cgrid,mPD)
cmocean('thermal')
colorbar
%caxis([0 2])
set(gca,'XTick',encs(2:2:end),'XTickLabel',encs(1:2:end),...
        'YTick',cmaxs,'YTickLabel',cmaxs)
xlabel('Enc coeff')
ylabel('Cmax coeff')
title('All metrics P&D')
print('-dpng',[pp cfile2 '_EncCmax_SAUP_comboPD.png'])

save([dp 'Bio_rates/' cfile2 '_clim_',harv,'_loop_SAUPcomp_EncCmax.mat']);


