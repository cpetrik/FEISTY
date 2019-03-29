% FEISTY 100 LHS parameter sets
% Only last 12 months of 150 years saved 

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
dp = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/Dc_D075_Sm025_nmort1_noCC_RE00100/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/param_ensemble/Dc_D075_Sm025_nmort1_noCC_RE00100/';

% param ensemble results
load([dp 'Climatol_ensemble_LHS100.mat']);

% climatol parameter set
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
dpath = ['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat']);

%%
nfile = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/';
load([nfile 'LHS_param15_100vals.mat']);

% sfile = sim{M};
% sname = sfile(83:end);

%% SAU comparison
[r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch);

%% Outputs
% Max like ratio test Sum of Squares (sum of squared residuals)?
%mean square error
mse_all = rmse_all.^2;
mse = rmse.^2;

%% Histograms of pure misfits
figure(10)
subplot(2,2,1)
hist(mis(:,2));
title('F')      % normal but shifted >0
xlim([-4 4])
subplot(2,2,2)
hist(mis(:,3));
title('P')      % slightly right-skewed
xlim([-3 3])
subplot(2,2,3)
hist(mis(:,4));
title('D')      % pretty normal
xlim([-2 2])
print('-dpng',[pp 'FEISTY_clim_SAUP_misfit.png'])

mis_combo = reshape(mis_all,100*45,5);

figure(11)
subplot(2,2,1)
hist(mis_combo(:,2));
title('F')      % model has many very low values
%xlim([-4 4])
subplot(2,2,2)
hist(mis_combo(:,3));
title('P')      % model has a few very low values
%xlim([-3 3])
subplot(2,2,3)
hist(mis_combo(:,4));
title('D')      % model has a few low values
%xlim([-2 2])
print('-dpng',[pp 'Ensemble_clim_SAUP_misfit.png'])

mis_combo_fn = reshape(mis_combo(:,2:4),4500*3,1);

figure(12)
subplot(2,2,1)
hist(mis_combo(:,2));
title('F')      % model has many very low values
%xlim([-4 4])
subplot(2,2,2)
hist(mis_combo(:,3));
title('P')      % model has a few very low values
%xlim([-3 3])
subplot(2,2,3)
hist(mis_combo(:,4));
title('D')      % model has a few low values
%xlim([-2 2])
subplot(2,2,4)
hist(mis_combo_fn);
title('F, P, & D')      % model has a few low values
%xlim([-2 2])
print('-dpng',[pp 'Ensemble_clim_SAUP_misfit2.png'])

%% Plot in 3d space
%x=P,y=D,z=F
%Sum of squares
figure(1)
subplot(2,2,1)
scatter3(ss_all(3,:),ss_all(4,:),ss_all(2,:)); hold on;
scatter3(ss(3),ss(4),ss(2),'k','filled'); 
xlabel('P SS')
ylabel('D SS')
zlabel('F SS')
subplot(2,2,2)
scatter(ss_all(3,:),ss_all(4,:)); hold on;
scatter(ss(3),ss(4),'k','filled');
xlabel('P SS')
ylabel('D SS')
subplot(2,2,3)
scatter(ss_all(3,:),ss_all(2,:)); hold on;
scatter(ss(3),ss(2),'k','filled');
xlabel('P SS')
ylabel('F SS')
subplot(2,2,4)
scatter(ss_all(4,:),ss_all(2,:)); hold on;
scatter(ss(4),ss(2),'k','filled');
xlabel('D SS')
ylabel('F SS')

%% MSE
figure(2)
subplot(2,2,1)
scatter3(mse_all(3,:),mse_all(4,:),mse_all(2,:)); hold on;
scatter3(mse(3),mse(4),mse(2),'k','filled');
xlabel('P MSE')
ylabel('D MSE')
zlabel('F MSE')
subplot(2,2,2)
scatter(mse_all(3,:),mse_all(4,:)); hold on;
scatter(mse(3),mse(4),'k','filled');
xlabel('P MSE')
ylabel('D MSE')
subplot(2,2,3)
scatter(mse_all(3,:),mse_all(2,:)); hold on;
scatter(mse(3),mse(2),'k','filled');
xlabel('P MSE')
ylabel('F MSE')
subplot(2,2,4)
scatter(mse_all(4,:),mse_all(2,:)); hold on;
scatter(mse(4),mse(2),'k','filled');
xlabel('D MSE')
ylabel('F MSE')

%% RMSE
figure(3)
subplot(2,2,1)
scatter3(rmse_all(3,:),rmse_all(4,:),rmse_all(2,:)); hold on;
scatter3(rmse(3),rmse(4),rmse(2),'k','filled');
xlabel('P RMSE')
ylabel('D RMSE')
zlabel('F RMSE')
subplot(2,2,2)
scatter(rmse_all(3,:),rmse_all(4,:)); hold on;
scatter(rmse(3),rmse(4),'k','filled');
xlabel('P RMSE')
ylabel('D RMSE')
xlim([0.5 2])
subplot(2,2,3)
scatter(rmse_all(3,:),rmse_all(2,:)); hold on;
scatter(rmse(3),rmse(2),'k','filled');
xlabel('P RMSE')
ylabel('F RMSE')
xlim([0.5 2])
subplot(2,2,4)
scatter(rmse_all(4,:),rmse_all(2,:)); hold on;
scatter(rmse(4),rmse(2),'k','filled');
xlabel('D RMSE')
ylabel('F RMSE')

%% low RMSEs
fid = find(rmse_all(2,:) < 5);
pid = find(rmse_all(3,:) < 0.9); %could go up to 1.1 to get lower Fs
did = find(rmse_all(4,:) < 0.6);
id1 = intersect(fid,did);
id2 = intersect(pid,id1);

pset = fx(id2,:);

%% RMSE
figure(4)
subplot(2,2,1)
scatter3(rmse_all(3,:),rmse_all(4,:),rmse_all(2,:)); hold on;
scatter3(rmse_all(3,id1),rmse_all(4,id1),rmse_all(2,id1),'b','filled'); hold on;
scatter3(rmse_all(3,id2),rmse_all(4,id2),rmse_all(2,id2),'r','filled'); hold on;
scatter3(rmse(3),rmse(4),rmse(2),'k','filled');
xlabel('P RMSE')
ylabel('D RMSE')
zlabel('F RMSE')
xlim([0.25 1.25])
ylim([0.25 1.25])
zlim([0 6])
subplot(2,2,2)
scatter(rmse_all(3,:),rmse_all(4,:)); hold on;
scatter(rmse_all(3,id1),rmse_all(4,id1),'b','filled'); hold on;
scatter(rmse_all(3,id2),rmse_all(4,id2),'r','filled'); hold on;
scatter(rmse(3),rmse(4),'k','filled');
xlabel('P RMSE')
ylabel('D RMSE')
xlim([0.25 1.25])
ylim([0.25 1.25])
subplot(2,2,3)
scatter(rmse_all(3,:),rmse_all(2,:)); hold on;
scatter(rmse_all(3,id1),rmse_all(2,id1),'b','filled'); hold on;
scatter(rmse_all(3,id2),rmse_all(2,id2),'r','filled'); hold on;
scatter(rmse(3),rmse(2),'k','filled');
xlabel('P RMSE')
ylabel('F RMSE')
xlim([0.25 1.25])
ylim([0 6])
subplot(2,2,4)
scatter(rmse_all(4,:),rmse_all(2,:)); hold on;
scatter(rmse_all(4,id1),rmse_all(2,id1),'b','filled'); hold on;
scatter(rmse_all(4,id2),rmse_all(2,id2),'r','filled'); hold on;
scatter(rmse(4),rmse(2),'k','filled');
xlabel('D RMSE')
ylabel('F RMSE')
xlim([0.25 1.25])
ylim([0 6])
print('-dpng',[pp 'RMSE_SAUP_type_best.png'])


%% realized param distr
figure(5)
for n=1:15
    subplot(4,4,n)
    hist(pset(:,n))
    xlim([plow(n) phi(n)])
    title(ptext{n})
end
print('-dpng',[pp 'param_distr_best.png'])

%% vis best maps
% for k=1:length(id2)
%     M=id2(k);
%     sfile = sim{M};
%     sname = sfile(83:end);
%     load([sfile '_means.mat']);
%     
%     %% Maps
%     vis_climatol_ensem(sname,sf_mean,sp_mean,sd_mean,...
%     mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean);
%     
% end


