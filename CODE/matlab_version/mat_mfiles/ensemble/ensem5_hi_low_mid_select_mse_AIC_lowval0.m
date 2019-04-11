% FEISTY 5 full parameter sets of hi or low
% Only last 12 months of 150 years saved 

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
dp = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';

% param ensemble results
load([dp 'Climatol_ensemble_param5_hi_low_mid_lowval0.mat'],'rmse_all','mis_all','sim',...
    'lme_Fmcatch','lme_Pmcatch','lme_Dmcatch');

% climatol parameter set
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
dpath = ['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat']);

%%
nfile = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/';
load([nfile 'LHS_param5_hi_low.mat']);
fx5 = fx;
clear fx 

load([nfile 'LHS_param5_hi_low_mid.mat'],'fx');
fx3 = fx;
clear fx 

fx = [fx5;fx3];

% sfile = sim{M};
% sname = sfile(83:end);

%% SAU comparison
[r,rmse,ss,mis] = lme_saup_corr_stock_ensem_lowval0(lme_mcatch);

% SAUP data
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/SAUP_Stock_top10.mat');
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;
%Set F values within log10(5) = 0

l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

l10all = [l10sF(keep);l10sP(keep);l10sD(keep)];
%variance of catch observations
sig = var(l10all);
%num of observations
n = length(l10all);

%% Outputs
%mean square error
mse_all = rmse_all.^2;
mse = rmse.^2;

%put residuals of all fn types in one vector
mis_all_F = mis_all(:,:,2);
mis_all_PD = mis_all(:,:,3:4);
mis_PD = reshape(mis_all_PD,length(fx),45*2);
mis_combo = [mis_all_F,mis_PD];

mis_fn = mis(:,2:4);
mis_fn = reshape(mis_fn,45*3,1);
mis_fn = mis_fn(~isnan(mis_fn));

%% Classic AIC 
% AIC = -2*log(L) + 2*K
% log(L) = (-n/2) * log(2*pi*var) - (1/(2*var)) * sum(resid^2)

%logLike
LL = -n/2 * log(2*pi*sig) - (1/(2*sig)) * sum(mis_combo.^2,2);
LL0 = -n/2 * log(2*pi*sig) - (1/(2*sig)) * sum(mis_fn.^2);
LL_all = [LL0;LL];

caic_all = -2 * LL;
caic = -2 * LL0;

[caic_srt,idc] = sort(caic_all);
idc2 = find(caic_srt<caic);

caic_srt2 = nan*ones(length(caic_srt)+1,1);
caic_srt2(idc2) = caic_srt(idc2);
caic_srt2(idc2(end)+1) = caic;
caic_srt2((idc2(end)+2):end) = caic_srt((idc2(end)+1):end);

idc_srt2 = nan*ones(length(caic_srt)+1,1);
idc_srt2(idc2) = idc(idc2);
idc_srt2(idc2(end)+1) = 0;
idc_srt2((idc2(end)+2):end) = idc((idc2(end)+1):end);

cdel = caic_srt2 - caic_srt2(1);
cw = exp(-0.5*cdel) ./ sum( exp(-0.5*cdel) );

caicv(:,1) = idc_srt2;
caicv(:,2) = caic_srt2;
caicv(:,3) = cdel;
caicv(:,4) = cw;
cT = array2table(caicv,'VariableNames',{'ParamSet','AIC','delta','weight'});
writetable(cT,[nfile 'LHS_param5_hi_low_mid_AIC_classic_lowval0.csv'])

%% Built in Fn
%logLike LL_all 
k = 5*ones(length(LL),1);
[baic_all,bbic_all] = aicbic(LL,k,n);
[baic,bbic] = aicbic(LL0,5,n);

[baic_srt,idb] = sort(baic_all);
baic_srt2 = [baic;baic_srt];
bdel = baic_srt2 - baic_srt2(1);
bw = exp(-0.5*bdel) ./ sum( exp(-0.5*bdel) );

% GIVES SAME ANSWER AS CLASSIC

baicv(:,1) = [0;idb];
baicv(:,2) = baic_srt2;
baicv(:,3) = bdel;
baicv(:,4) = bw;
bT = array2table(baicv,'VariableNames',{'ParamSet','AIC','delta','weight'});
writetable(bT,[nfile 'LHS_param5_hi_low_mid_AIC_builtin_lowval0.csv'])


%% AICs <= AIC(orig) + 2
besti = find(baicv(:,2) <= baic+2);
pid = baicv(besti,1);
pid = pid(pid>0);
pset(:,1) = pid;
pset(:,2:6) = fx(pid,:);
pset(:,7) = baic_all(pid);

pT = array2table(pset,'VariableNames',{'ParamSet','Lambda','bMet','bEnc',...
    'aMet','aEnc','AIC'});
writetable(pT,[nfile 'LHS_param5_hi_low_mid_bestAIC_params_lowval0.csv'])

id1 = pid;

%% Plot RMSE in 3D space
figure(4)
subplot(2,2,1)
scatter3(rmse_all(3,:),rmse_all(4,:),rmse_all(2,:)); hold on;
scatter3(rmse_all(3,id1),rmse_all(4,id1),rmse_all(2,id1),'b','filled'); hold on;
scatter3(rmse(3),rmse(4),rmse(2),'k','filled');
xlabel('P RMSE')
ylabel('D RMSE')
zlabel('F RMSE')
% xlim([0.25 8])
% ylim([0.25 1])
% zlim([0 6])
subplot(2,2,2)
scatter(rmse_all(3,:),rmse_all(4,:)); hold on;
scatter(rmse_all(3,id1),rmse_all(4,id1),'b','filled'); hold on;
scatter(rmse(3),rmse(4),'k','filled');
xlabel('P RMSE')
ylabel('D RMSE')
% xlim([0.25 8])
% ylim([0.25 1])
subplot(2,2,3)
scatter(rmse_all(3,:),rmse_all(2,:)); hold on;
scatter(rmse_all(3,id1),rmse_all(2,id1),'b','filled'); hold on;
scatter(rmse(3),rmse(2),'k','filled');
xlabel('P RMSE')
ylabel('F RMSE')
% xlim([0.25 8])
% ylim([0 6])
subplot(2,2,4)
scatter(rmse_all(4,:),rmse_all(2,:)); hold on;
scatter(rmse_all(4,id1),rmse_all(2,id1),'b','filled'); hold on;
scatter(rmse(4),rmse(2),'k','filled');
xlabel('D RMSE')
ylabel('F RMSE')
% xlim([0.25 1])
% ylim([0 6])
print('-dpng',[pp 'RMSE_SAUP_type_best_AIC_param5_hi_low_mid_lowval0.png'])


%% realized param distr
figure(5)
for n=1:5
    subplot(3,2,n)
    hist(pset(:,n+1))
    xlim([plow(n) phi(n)])
    title(ptext{n})
end
print('-dpng',[pp 'param5_hi_lo_mid_distr_best_AIC_lowval0.png'])

% Assim = 0.75; amet = 3; gamma = 100 more often

%% vis best maps
% for j=1:length(id1)
%     M=id1(j);
%     sfile = sim{M};
%     sname = sfile(140:end);
%     load([sfile '_means.mat']);
%     
%     %% Maps
%     vis_climatol_ensem5(sname,sf_mean,sp_mean,sd_mean,...
%     mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean,pp);
%     
% end

% P & D ARE BAD IF A_MET = 6

%% vis best SAUP comp
% for j=1:length(id1)
%     M=id1(j);
%     sfile = sim{M};
%     sname = sfile(140:end);
%     
%     %% Comp
%     vis_clim_lme_saup_corr_stock(lme_Fmcatch(:,M),lme_Pmcatch(:,M),...
%         lme_Dmcatch(:,M),pp,sname);
%     
% end

%% vis good and bad F SAUP comp
% bad = [32,28,13];
% good = [2,46,42,38,43,39,40];
% Fbg = [bad,good];
% for j=1:length(Fbg)
%     M=Fbg(j);
%     sfile = sim{M};
%     sname = sfile(140:end);
%     
%     %% Comp
%     clim_catch_lme_saup_Fcorr_stock_LMEid(lme_Fmcatch(:,M),pp,sname);
%     
% end


%% Plot AIC in 3D space
% reduce param to bpow=0.15 and amet=3
idb = find(fx(:,2)==0.15);
ida = find(fx(:,4)==3);
idbot = intersect(ida,idb);
test = fx(idbot,:);
test(:,6) = baic_all(idbot);

x = test(:,1); 
y = test(:,3); 
z = test(:,5); 
v = test(:,6); 

lam = [0.6,0.675,0.75];
be = [0.15,0.2,0.25];
ae = [50,75,100];
[xq,yq,zq] = meshgrid(lam,be,ae);
% Interpolate the scattered data on the grid. Plot the results.
vq = griddata(x,y,z,v,xq,yq,zq);

%%
agrid = NaN(4,4,4);
bgrid = NaN(4,4,4);
cgrid = NaN(4,4,4);
egrid = NaN(4,4,4);

agrid(1:3,1:3,1:3) = xq;
bgrid(1:3,1:3,1:3) = yq;
egrid(1:3,1:3,1:3) = zq;
cgrid(1:3,1:3,1:3) = vq;

%%
figure(6)
subplot(2,3,1)
scatter3(test(:,1),test(:,3),test(:,5),100,test(:,6),'filled'); hold on;
xlabel('Assim')
ylabel('b_E')
zlabel('a_E')
colormap('jet')
caxis([360 385])
colorbar('northoutside')
%
subplot(2,3,4)
surf(squeeze(agrid(:,1,:)),squeeze(bgrid(:,1,:)),squeeze(egrid(:,1,:)),...
    squeeze(cgrid(:,1,:)),'FaceColor','interp');
hold on;
surf(squeeze(agrid(:,2,:)),squeeze(bgrid(:,2,:)),squeeze(egrid(:,2,:)),...
    squeeze(cgrid(:,2,:)),'FaceColor','interp');
hold on;
surf(squeeze(agrid(:,3,:)),squeeze(bgrid(:,3,:)),squeeze(egrid(:,3,:)),...
    squeeze(cgrid(:,3,:)),'FaceColor','interp');
colormap('jet')
caxis([360 385])
%
subplot(2,3,5)
surf(squeeze(agrid(1,:,:)),squeeze(bgrid(1,:,:)),squeeze(egrid(1,:,:)),...
    squeeze(cgrid(1,:,:)),'FaceColor','interp');
hold on;
surf(squeeze(agrid(2,:,:)),squeeze(bgrid(2,:,:)),squeeze(egrid(2,:,:)),...
    squeeze(cgrid(2,:,:)),'FaceColor','interp');
hold on;
surf(squeeze(agrid(3,:,:)),squeeze(bgrid(3,:,:)),squeeze(egrid(3,:,:)),...
    squeeze(cgrid(3,:,:)),'FaceColor','interp');
colormap('jet')
caxis([360 385])
%
subplot(2,3,6)
surf(squeeze(agrid(:,:,1)),squeeze(bgrid(:,:,1)),squeeze(egrid(:,:,1)),...
    squeeze(cgrid(:,:,1)),'FaceColor','interp');
hold on;
surf(squeeze(agrid(:,:,2)),squeeze(bgrid(:,:,2)),squeeze(egrid(:,:,2)),...
    squeeze(cgrid(:,:,2)),'FaceColor','interp');
hold on;
surf(squeeze(agrid(:,:,3)),squeeze(bgrid(:,:,3)),squeeze(egrid(:,:,3)),...
    squeeze(cgrid(:,:,3)),'FaceColor','interp');
colormap('jet')
caxis([360 385])

%
subplot(2,3,2)
surf(squeeze(agrid(:,1,:)),squeeze(bgrid(:,1,:)),squeeze(egrid(:,1,:)),...
    squeeze(cgrid(:,1,:)),'FaceColor','interp');
hold on;

surf(squeeze(agrid(1,:,:)),squeeze(bgrid(1,:,:)),squeeze(egrid(1,:,:)),...
    squeeze(cgrid(1,:,:)),'FaceColor','interp');
hold on;

surf(squeeze(agrid(:,:,3)),squeeze(bgrid(:,:,3)),squeeze(egrid(:,:,3)),...
    squeeze(cgrid(:,:,3)),'FaceColor','interp');
colormap('jet')
caxis([360 385])
title('AIC when b_M = 0.15 and a_M = 3')

%
subplot(2,3,3)
surf(squeeze(agrid(:,3,:)),squeeze(bgrid(:,3,:)),squeeze(egrid(:,3,:)),...
    squeeze(cgrid(:,3,:)),'FaceColor','interp');
hold on;

surf(squeeze(agrid(3,:,:)),squeeze(bgrid(3,:,:)),squeeze(egrid(3,:,:)),...
    squeeze(cgrid(3,:,:)),'FaceColor','interp');
hold on;

surf(squeeze(agrid(:,:,1)),squeeze(bgrid(:,:,1)),squeeze(egrid(:,:,1)),...
    squeeze(cgrid(:,:,1)),'FaceColor','interp');
hold on;
colormap('jet')
caxis([360 385])
print('-dpng',[pp 'param5_hi_lo_mid_AIC_best_met_lowval0.png'])






