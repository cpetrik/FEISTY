% All correlations
clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
%keep = notLELC;
keep=1:66;

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath = [dp cfile '/'];
ppath = [pp cfile '/'];

%% Hist grid
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t','AREA_OCN');
grid = csvread([gpath 'grid_csv.csv']);
load([gpath 'lme_mask_esm2m.mat']);
load([cpath 'LME_hist_temp_zoop_det.mat'],'lme_ptemp','lme_area');

hlme_ptemp = lme_ptemp;
clear lme_ptemp

hlme_area_km2 = lme_area * 1e-6;
clear lme_area

hlme = lme_mask_esm2m';
hID = grid(:,1);

%% Fore grid
load([cpath 'LME_fore_temp.mat'],'lme_ptemp');

flme_ptemp = lme_ptemp;
clear lme_ptemp
 
%% FEISTY Hist 1951-2000
load([fpath 'LME_hist_',harv,'_' cfile '.mat'],'lme_mcatch',...
    'lme_mbio','lme_sbio');
load([fpath 'TEeffDet_Historic_All_fish03_' cfile '.mat'],'lme_te');

hlme_mcatch = nansum(lme_mcatch,2) * 1e-6;
hlme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
hlme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
hlme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
hlme_Bmbio = lme_mbio(:,9) * 1e-6;
hlme_Mmbio = (lme_mbio(:,4) + lme_mbio(:,5) + lme_mbio(:,6)) * 1e-6;
hlme_Lmbio = (lme_mbio(:,7) + lme_mbio(:,8)) * 1e-6;

% MT/km2
hlme_mcatch = hlme_mcatch ./ hlme_area_km2;
hlme_Fmcatch = hlme_Fmcatch ./ hlme_area_km2;
hlme_Pmcatch = hlme_Pmcatch ./ hlme_area_km2;
hlme_Dmcatch = hlme_Dmcatch ./ hlme_area_km2;
hlme_Bmbio = hlme_Bmbio ./ hlme_area_km2;
hlme_Mmbio = hlme_Mmbio ./ hlme_area_km2;
hlme_Lmbio = hlme_Lmbio ./ hlme_area_km2;

hFracPD = hlme_Pmcatch ./ (hlme_Pmcatch + hlme_Dmcatch);
hFracPF = hlme_Pmcatch ./ (hlme_Pmcatch + hlme_Fmcatch);
hFracLM = hlme_Lmbio ./ (hlme_Lmbio + hlme_Mmbio);

l10h=log10(hlme_mcatch);
l10hF=log10(hlme_Fmcatch);
l10hP=log10(hlme_Pmcatch);
l10hD=log10(hlme_Dmcatch);
l10hB=log10(hlme_Bmbio);

l10hL = log10(lme_te(:,2));
l10hHTL = log10(lme_te(:,4));
l10hLTL = log10(lme_te(:,6));

clear lme_mcatch lme_mbio lme_sbio lme_te

%% FEISTY Fore 2051-2100
load([fpath 'LME_fore_',harv,'_' cfile '.mat'],'lme_mcatch',...
    'lme_mbio','lme_sbio');
load([fpath 'TEeffDet_Forecast_All_fish03_' cfile '.mat'],'lme_te');

% POEM LME biomass in MT
plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
plme_Bmbio = lme_mbio(:,9) * 1e-6;
plme_Mmbio = (lme_mbio(:,4) + lme_mbio(:,5) + lme_mbio(:,6)) * 1e-6;
plme_Lmbio = (lme_mbio(:,7) + lme_mbio(:,8)) * 1e-6;
% MT/km2
plme_mcatch = plme_mcatch ./ hlme_area_km2;
plme_Fmcatch = plme_Fmcatch ./ hlme_area_km2;
plme_Pmcatch = plme_Pmcatch ./ hlme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ hlme_area_km2;
plme_Bmbio = plme_Bmbio ./ hlme_area_km2;
plme_Mmbio = plme_Mmbio ./ hlme_area_km2;
plme_Lmbio = plme_Lmbio ./ hlme_area_km2;

pFracPD = plme_Pmcatch ./ (plme_Pmcatch + plme_Dmcatch);
pFracPF = plme_Pmcatch ./ (plme_Pmcatch + plme_Fmcatch);
pFracLM = plme_Lmbio ./ (plme_Lmbio + plme_Mmbio);

l10p=log10(plme_mcatch);
l10pF=log10(plme_Fmcatch);
l10pP=log10(plme_Pmcatch);
l10pD=log10(plme_Dmcatch);
l10pB=log10(plme_Bmbio);

% FEISTY LME TEeffs
%     lme_te(L,2) = nanmean(TEeff_L
%     lme_te(L,4) = nanmean(TEeff_HTLd
%     lme_te(L,6) = nanmean(TEeff_LTLd
l10pL = log10(lme_te(:,2));
l10pHTL = log10(lme_te(:,4));
l10pLTL = log10(lme_te(:,6));

clear lme_mcatch lme_mbio lme_sbio lme_te


%% Percent difference?
% [rall,pall]=corr(l10h(keep),l10p(keep));
% [rF,pF]=corr(l10hF(keep),l10pF(keep));
% [rP,pP]=corr(l10hP(keep),l10pP(keep));
% [rD,pD]=corr(l10hD(keep),l10pD(keep));
% [rB,pB]=corr(l10hB(keep),l10pB(keep));
% [rPD,pPD]=corr(hFracPD(keep),pFracPD(keep));
% [rPF,pPF]=corr(hFracPF(keep),pFracPF(keep));
% [rLM,pLM]=corr(hFracLM(keep),pFracLM(keep));
% [rL,pL]=corr(l10hL(keep),l10pL(keep));
% [rHTL,pHTL]=corr(l10hHTL(keep),l10pHTL(keep));
% [rLTL,pLTL]=corr(l10hLTL(keep),l10pLTL(keep));
% 
% 
% % Table
% fish_stat(1,1) = rall;
% fish_stat(2,1) = rF;
% fish_stat(3,1) = rP;
% fish_stat(4,1) = rD;
% fish_stat(5,1) = rB;
% fish_stat(6,1) = rPD;
% fish_stat(7,1) = rPF;
% fish_stat(8,1) = rLM;
% fish_stat(9,1) = rL;
% fish_stat(10,1) = rHTL;
% fish_stat(11,1) = rLTL;
% 
% % save
% Fstat = array2table(fish_stat,'VariableNames',{'r'},...
%     'RowNames',{'All Fish','F','P','D','B','Frac Pel-Dem','Frac Pel-Forage',...
%     'Frac Large-Med','TEeffL','TEeffHTL','TEeffLTL'});
% writetable(Fstat,[fpath 'LME_hist_fore_stats_' cfile '.csv'],'Delimiter',',',...
%     'WriteRowNames',true)
% save([fpath 'LME_hist_fore_stats_' cfile '.mat'],'fish_stat')

%% Figures
x=-8:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

figure(1)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hF(keep),l10pF(keep),20,flme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.375 0.5 0.3 0.025],'orientation','horizontal')
%text(-2.75,0.75,['r = ' sprintf('%2.2f',rF) ' (p = ' sprintf('%2.2f',pF) ')'])
%text(-2.75,0.5,['RMSE = ' sprintf('%2.2f',rmseF)])
axis([-3 1 -3 1])
xlabel('Hist')
ylabel('RCP 8.5')
title('Forage Fishes')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hP(keep),l10pP(keep),20,flme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
%text(-5.5,1.0,['r = ' sprintf('%2.2f',rP) ' (p = ' sprintf('%2.2f',pP) ')'])
%text(-5.5,0.5,['RMSE = ' sprintf('%2.2f',rmseP)])
axis([-6 2 -6 2])
xlabel('Hist')
ylabel('RCP 8.5')
title('Large Pelagics')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hD(keep),l10pD(keep),20,flme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
%text(-1.75,1.4,['r = ' sprintf('%2.2f',rD) ' (p = ' sprintf('%2.2f',pD) ')'])
%text(-1.75,1.1,['RMSE = ' sprintf('%2.2f',rmseD)])
axis([-2 1 -2 1])
xlabel('Hist')
ylabel('RCP 8.5')
title('Demersals')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10h(keep),l10p(keep),20,flme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
%text(-1.75,1.4,['r = ' sprintf('%2.2f',rall) ' (p = ' sprintf('%2.2f',pall) ')'])
%text(-1.75,1.1,['RMSE = ' sprintf('%2.2f',rmse)])
axis([-2 1 -2 1])
xlabel('Hist')
ylabel('RCP 8.5')
title('All fishes')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Hist_RCP85_',harv,'_comp_types_temp.png'])

%% Fractions
figure(2)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(hFracPF(keep),pFracPF(keep),20,flme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.375 0.5 0.3 0.025],'orientation','horizontal')
%text(-0.1,1.1,['r = ' sprintf('%2.2f',rPF) ' (p = ' sprintf('%2.2f',pPF) ')'])
%text(-0.1,1.0,['RMSE = ' sprintf('%2.2f',rmsePF)])
axis([-0.2 1.2 -0.2 1.2])
xlabel('Hist')
ylabel('RCP 8.5')
title('P / (P+F)')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(hFracPD(keep),pFracPD(keep),20,flme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
%text(-0.1,1.1,['r = ' sprintf('%2.2f',rPD) ' (p = ' sprintf('%2.2f',pPD) ')'])
%text(-0.1,1.0,['RMSE = ' sprintf('%2.2f',rmsePD)])
axis([-0.2 1.2 -0.2 1.2])
xlabel('Hist')
ylabel('RCP 8.5')
title('P / (P+D)')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(hFracLM(keep),pFracLM(keep),20,flme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
%text(-0.1,1.1,['r = ' sprintf('%2.2f',rLM) ' (p = ' sprintf('%2.2f',pLM) ')'])
%text(-0.1,1.0,['RMSE = ' sprintf('%2.2f',rmseLM)])
axis([-0.2 1.2 -0.2 1.2])
xlabel('Hist')
ylabel('RCP 8.5')
title('L / (L+M)')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hB(keep),l10pB(keep),20,flme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
%text(-3.75,-0.3,['r = ' sprintf('%2.2f',rB) ' (p = ' sprintf('%2.2f',pB) ')'])
%text(-3.75,-0.7,['RMSE = ' sprintf('%2.2f',rmseB)])
axis([-4 0 -4 0])
xlabel('Hist')
ylabel('RCP 8.5')
title('Benthos')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Hist_RCP85_',harv,'_comp_fracs_temp.png'])

% benthos figs look the same scale, so mistake somewhere else

%% TEeffs
figure(3)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hL(keep),l10pL(keep),20,flme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.375 0.5 0.3 0.025],'orientation','horizontal')
%text(-4.75,-1.25,['r = ' sprintf('%2.2f',rL) ' (p = ' sprintf('%2.2f',pL) ')'])
%text(-4.75,-1.5,['RMSE = ' sprintf('%2.2f',rmseL)])
axis([-5 -2 -5 -2])
xlabel('Hist')
ylabel('RCP 8.5')
title('log_1_0 TEeff L')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hHTL(keep),l10pHTL(keep),20,flme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
%text(-3.25,-0.75,['r = ' sprintf('%2.2f',rHTL) ' (p = ' sprintf('%2.2f',pHTL) ')'])
%text(-3.25,-1.0,['RMSE = ' sprintf('%2.2f',rmseHTL)])
axis([-4 -0.5 -4 -0.5])
xlabel('Hist')
ylabel('RCP 8.5')
title('log_1_0 TEeff HTL')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hLTL(keep),l10pLTL(keep),20,flme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
%text(-2.4,-0.65,['r = ' sprintf('%2.2f',rLTL) ' (p = ' sprintf('%2.2f',pLTL) ')'])
%text(-2.4,-0.8,['RMSE = ' sprintf('%2.2f',rmseLTL)])
axis([-1.8 -0.8 -1.8 -0.8])
xlabel('Hist')
ylabel('RCP 8.5')
title('log_1_0 TEeff LTL')

stamp([harv '_' cfile])
print('-dpng',[ppath 'Hist_RCP85_',harv,'_comp_TEeffs_temp.png'])