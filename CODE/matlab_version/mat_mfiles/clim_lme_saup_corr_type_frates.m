%POEM catch vs. SAUP catch by LME
%Find fishing selectivity that gives highest corr for each type

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

%use weighted catches
load([spath 'SAUP_LME_Catch_annual.mat'],'yr','totcatch','lme_catch',...
    'Flme_wcatch','Dlme_wcatch','Plme_wcatch');

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')
cmap1(1,:)=[1 1 1];
cmap1(2,:)=cmap_ppt(1,:);
cmap1(3,:)=cmap_ppt(3,:);
cmap1(4,:)=cmap_ppt(5,:);

cmap2(1,:)=cmap_ppt(1,:);
cmap2(2,:)=cmap_ppt(3,:);
cmap2(3,:)=cmap_ppt(5,:);

cmap3(1,:)=[0 0 0];
cmap3(2,:)=cmap_ppt(1,:);
cmap3(3,:)=cmap_ppt(3,:);
cmap3(4,:)=cmap_ppt(5,:);

cmap4(1,:)=cmap_ppt(3,:);
cmap4(2,:)=cmap_ppt(1,:);
cmap4(3,:)=cmap_ppt(5,:);

AREA_OCN = max(area,1);

lfile = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100';
lpath = [dp lfile '/'];
load([lpath 'LME_clim_fished03_' lfile '.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;

% plot info
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

%%
Flme_catch_all = nansum(Flme_wcatch,3);
Plme_catch_all = nansum(Plme_wcatch,3);
Dlme_catch_all = nansum(Dlme_wcatch,3);

%1956-2005 SAUP average
id = find(yr>1955 & yr<=2005);

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
%Top 10 yrs SAUP
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

%% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);

keep = [1:9,11:34,36:38,46:54,59:60,62,65];

x=-6:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% Loop over POEM params
fqs = 0.5:0.5:3.5;
pqs = 0.25:0.25:1.5;
dqs = 1:7;
frate = 0.1;

%rall = NaN*ones(length(fqs),length(pqs),length(dqs));
rall = NaN*ones(length(fqs),1,1);
rF = rall;
rP = rall;
rD = rall;
rmse = rall;
rmseF = rall;
rmseP = rall;
rmseD = rall;

sname = 'Clim_';
dp = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_CC100_lgRE00100_mdRE00100';
datap = '/Volumes/GFDL/CSV/Matlab_new_size/';
dpath = [datap char(dp) '/'];
ppath = [pp char(dp) '/'];

for fq=1:length(fqs)
    for pq=4;%1:length(pqs)
        for dq=7;%1:length(dqs)
            MFsel = fqs(fq);
            LPsel = pqs(pq);
            LDsel = dqs(dq);
            tF = num2str(1000+int64(100*frate*MFsel));
            tP = num2str(1000+int64(100*frate*LPsel));
            tD = num2str(1000+int64(100*frate*LDsel));
            
            charv = ['fish_lF',tF(2:end),'_qP',tP(2:end),'_qD',tD(2:end)];
            harv = ['lF',tF(2:end),'_qP',tP(2:end),'_qD',tD(2:end)];
            cfile = [harv,'_cm20_m-b175-k09_D075_J100_A050_Sm025_nmort1_BE05_CC100_RE00100'];
            
            if exist([dpath 'LME_clim_fished_',harv,'_loop_' dp '.mat'],'file')
                load([dpath 'LME_clim_fished_',harv,'_loop_' dp '.mat'],'lme_mcatch');
                
                close all
                
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
                
                %log10 Difference
                diff_catch = (log10(plme_mcatch) - log10(slme_mcatch10)) ./ log10((slme_mcatch10+eps));
                Fdiff_catch = (log10(plme_Fmcatch) - log10(Flme_mcatch10)) ./ log10((Flme_mcatch10+eps));
                Pdiff_catch = (log10(plme_Pmcatch) - log10(Plme_mcatch10)) ./ log10((Plme_mcatch10+eps));
                Ddiff_catch = (log10(plme_Dmcatch) - log10(Dlme_mcatch10)) ./ log10((Dlme_mcatch10+eps));
                
                code = [1:66]';
                T = table(code,slme_mcatch10,plme_mcatch,Flme_mcatch10,plme_Fmcatch,...
                    Plme_mcatch10,plme_Pmcatch,Dlme_mcatch10,plme_Dmcatch,'VariableNames',...
                    {'lme','saupAll','poemAll','saupF','poemF','saupP','poemP','saupD','poemD'});
                writetable(T,[dpath 'LME_saup_catch_clim_fished',harv,'.csv']);
                
                %% Plots
                l10s=log10(slme_mcatch10+eps);
                l10p=log10(plme_mcatch);
                l10sF=log10(Flme_mcatch10+eps);
                l10pF=log10(plme_Fmcatch);
                l10sP=log10(Plme_mcatch10+eps);
                l10pP=log10(plme_Pmcatch);
                l10sD=log10(Dlme_mcatch10+eps);
                l10pD=log10(plme_Dmcatch);
                
                
                %% Drop Arctic, Antarctic, Hawaii, Australia -------------------------
                
                % Stats
                %r
                rall(fq,1,1)=corr(l10s(keep),l10p(keep));
                rF(fq,1,1)=corr(l10sF(keep),l10pF(keep));
                rP(fq,1,1)=corr(l10sP(keep),l10pP(keep));
                rD(fq,1,1)=corr(l10sD(keep),l10pD(keep));
                
                %root mean square error
                o=l10s(keep);
                p=l10p(keep);
                n = length(o);
                num=nansum((p-o).^2);
                rmse(fq,1,1) = sqrt(num/n);
                
                o=l10sF(keep);
                p=l10pF(keep);
                n = length(o);
                num=nansum((p-o).^2);
                rmseF(fq,1,1) = sqrt(num/n);
                
                o=l10sP(keep);
                p=l10pP(keep);
                n = length(o);
                num=nansum((p-o).^2);
                rmseP(fq,1,1) = sqrt(num/n);
                
                o=l10sD(keep);
                p=l10pD(keep);
                n = length(o);
                num=nansum((p-o).^2);
                rmseD(fq,1,1) = sqrt(num/n);
                
            end %if
            
        end
    end
end

%% Plots of r and RMSE
cfile2 = ['_cm20_m-b175-k09_D075_J100_A050_Sm025_nmort1_BE05_CC100_RE00100'];

nf = length(fqs);
np = length(pqs);
nd = length(dqs);

%%
figure(10)
subplot(2,2,1)
bar(squeeze(rall))
xlabel('F sel')
title('Corr all fish')

subplot(2,2,2)
bar(squeeze(rF))
xlabel('F sel')
title('Corr F')

subplot(2,2,3)
bar(squeeze(rP))
xlabel('F sel')
title('Corr P')

subplot(2,2,4)
bar(squeeze(rD))
xlabel('F sel')
title('Corr D')
stamp(cfile2)
%print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_corrD_frates_PD'])

%%
figure(22)
subplot(2,2,1)
bar(squeeze(rmse))
xlabel('F sel')
title('RMSE all fish')

subplot(2,2,2)
bar(squeeze(rmseF))
xlabel('F sel')
title('RMSE F')

subplot(2,2,3)
bar(squeeze(rmseP))
xlabel('F sel')
title('RMSE P')

subplot(2,2,4)
bar(squeeze(rmseD))
xlabel('F sel')
title('RMSE D')
stamp(cfile2)
%print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_rmseD_frates_PD'])



%%
vals = NaN*ones(2,6);
qr = NaN*ones(2,6);

corr_all = 0.25*(normalize(rall)+normalize(rF)+normalize(rP)+normalize(rD));
rmse_all = normalize(rmse)+normalize(rmseF)+normalize(rmseP)+normalize(rmseD);
corr_PD = 0.5*(normalize(rP)+normalize(rD));
rmse_PD = normalize(rmseP)+normalize(rmseD);

ax1 = find(rall(:)==max(rall(:)));
fx1 = find(rF(:)==max(rF(:)));
px1 = find(rP(:)==max(rP(:)));
dx1 = find(rD(:)==max(rD(:)));
lx1 = find(corr_all(:)==max(corr_all(:)));
pdx1 = find(corr_PD(:)==max(corr_PD(:)));

ax2 = find(rmse(:)==min(rmse(:)));
fx2 = find(rmseF(:)==min(rmseF(:)));
px2 = find(rmseP(:)==min(rmseP(:)));
dx2 = find(rmseD(:)==min(rmseD(:)));
lx2 = find(rmse_all(:)==min(rmse_all(:)));
pdx2 = find(rmse_PD(:)==min(rmse_PD(:)));

vals(1,1) = rall(ax1);
vals(1,2) = rF(fx1);
vals(1,3) = rP(px1);
vals(1,4) = rD(dx1);
vals(1,5) = corr_all(lx1);
vals(1,6) = corr_PD(pdx1);
vals(2,1) = rmse(ax2);
vals(2,2) = rmseF(fx2);
vals(2,3) = rmseP(px2);
vals(2,4) = rmseD(dx2);
vals(2,5) = rmse_all(lx2);
vals(2,6) = rmse_PD(pdx2);

qr(1,1) = fqs(ax1);
qr(1,2) = fqs(fx1);
qr(1,3) = fqs(px1);
qr(1,4) = fqs(dx1);
qr(1,5) = fqs(lx1);
qr(1,6) = fqs(pdx1);

qr(2,1) = fqs(ax2);
qr(2,2) = fqs(fx2);
qr(2,3) = fqs(px2);
qr(2,4) = fqs(dx2);
qr(2,5) = fqs(lx2);
qr(2,6) = fqs(pdx2);

qr(:,7) = mean(qr,2);

% %%
% save([dpath 'Clim_fished_loop_LME_SAUP_catch_comp_qs.mat'],'rall','rF','rP',...
%     'rD','rmse','rmseF','rmseP','rmseD','fgrid3','pgrid3','dgrid3',...
%     'vals','qr','qrmse')
% csvwrite([dpath 'Clim_fished_loop_LME_SAUP_catch_comp_qs_corr.csv'],qr);
% csvwrite([dpath 'Clim_fished_loop_LME_SAUP_catch_comp_qs_rmse.csv'],qrmse);

