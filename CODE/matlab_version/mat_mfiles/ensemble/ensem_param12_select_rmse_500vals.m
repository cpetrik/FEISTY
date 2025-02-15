% FEISTY 100 LHS parameter sets
% Only last 12 months of 150 years saved 
% Reduced parameter space based on results of 5 param hi/low runs

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
dp = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/Dc_D075_J100_Sm025_nmort1_noCC_RE00100/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/param_ensemble/Dc_D075_J100_Sm025_nmort1_noCC_RE00100/';

% param ensemble results
load([dp 'Climatol_ensemble_param12_LHS500.mat']);

% climatol parameter set
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
dpath = ['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat']);

%%
nfile = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/';
load([nfile 'LHS_param12_500vals.mat']);

% sfile = sim{M};
% sname = sfile(83:end);

%% SAU comparison
[r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch);

%% Outputs
% Max like ratio test Sum of Squares (sum of squared residuals)?
%mean square error
mse_all = rmse_all.^2;
mse = rmse.^2;

%% Plot in 3d space
%x=P,y=D,z=F

% RMSE
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
%xlim([0.5 2])
subplot(2,2,3)
scatter(rmse_all(3,:),rmse_all(2,:)); hold on;
scatter(rmse(3),rmse(2),'k','filled');
xlabel('P RMSE')
ylabel('F RMSE')
axis([0 4 0 4])
subplot(2,2,4)
scatter(rmse_all(4,:),rmse_all(2,:)); hold on;
scatter(rmse(4),rmse(2),'k','filled');
xlabel('D RMSE')
ylabel('F RMSE')

%% low RMSEs
fid = find(rmse_all(2,:) < 2);
pid = find(rmse_all(3,:) < 0.9); 
did = find(rmse_all(4,:) < 0.9);
id1 = intersect(fid,did);
id2 = intersect(pid,id1);

pset = fx2(id2,:);

pset2 = pset;
pset2(:,13) = rmse_all(2,id2)';
pset2(:,14) = rmse_all(3,id2)';
pset2(:,15) = rmse_all(4,id2)';

pT = array2table(pset2,'VariableNames',{'Lambda','K_a','amet','h','gam','kce',...
    'kt','bpow','benc','bcmx','bent_eff','A','rmseF','rmseP','rmseD'});
writetable(pT,[nfile 'LHS500_param12_bestRMSE_params.csv'])

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
xlim([0.25 8])
ylim([0.25 1])
zlim([0 6])
subplot(2,2,2)
scatter(rmse_all(3,:),rmse_all(4,:)); hold on;
scatter(rmse_all(3,id1),rmse_all(4,id1),'b','filled'); hold on;
scatter(rmse_all(3,id2),rmse_all(4,id2),'r','filled'); hold on;
scatter(rmse(3),rmse(4),'k','filled');
xlabel('P RMSE')
ylabel('D RMSE')
xlim([0.25 8])
ylim([0.25 1])
subplot(2,2,3)
scatter(rmse_all(3,:),rmse_all(2,:)); hold on;
scatter(rmse_all(3,id1),rmse_all(2,id1),'b','filled'); hold on;
scatter(rmse_all(3,id2),rmse_all(2,id2),'r','filled'); hold on;
scatter(rmse(3),rmse(2),'k','filled');
xlabel('P RMSE')
ylabel('F RMSE')
xlim([0.25 8])
ylim([0 6])
subplot(2,2,4)
scatter(rmse_all(4,:),rmse_all(2,:)); hold on;
scatter(rmse_all(4,id1),rmse_all(2,id1),'b','filled'); hold on;
scatter(rmse_all(4,id2),rmse_all(2,id2),'r','filled'); hold on;
scatter(rmse(4),rmse(2),'k','filled');
xlabel('D RMSE')
ylabel('F RMSE')
xlim([0.25 1])
ylim([0 6])
print('-dpng',[pp 'RMSE_SAUP_type_best_param12_500.png'])


%% realized param distr
figure(5)
for n=1:12
    subplot(4,3,n)
    hist(pset(:,n))
    xlim([plow(n) phi(n)])
    title(ptext{n})
end
print('-dpng',[pp 'param_distr_bestRMSE_param12_500.png'])

%% vis best maps
for k=1:length(id2)
    M=id2(k);
    sfile = sim{M};
    sname = sfile(88:end);
    load([sfile '_means.mat']);
    
    %% Maps
    vis_climatol_ensem13(sname,sf_mean,sp_mean,sd_mean,...
    mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean,pp);
    
end

%% vis maps w/ good F & D, but not P
id3 = setdiff(id1,id2);
pset3 = fx2(id3,:);

% realized param distr
figure(6)
for n=1:12
    subplot(4,4,n)
    hist(pset3(:,n))
    xlim([plow(n) phi(n)])
    title(ptext{n})
end
%print('-dpng',[pp 'param_distr_best_FD_param12_500.png'])

%%
% for k=1:length(id3)
%     M=id3(k);
%     sfile = sim{M};
%     sname = sfile(83:end);
%     load([sfile '_means.mat']);
%     
%     %% Maps
%     vis_climatol_ensem(sname,sf_mean,sp_mean,sd_mean,...
%     mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean);
%     
% end


