% All correlations between Clim and NoNuUpdate Clim
clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

dp = '/Volumes/FEISTY/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
%keep = notLELC;
keep=1:66;

harv = 'All_fish03';

%% Climatol grid
Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
cdir='/Volumes/FEISTY/GCM_DATA/ESM26_hist/';
load([gpath 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([gpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([gpath 'LME_clim_temp_zoop_det.mat']);

clme = lme_mask_onedeg;

clme_ptemp = lme_ptemp;
clear lme_ptemp

%% FEISTY Climatol
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath = [dp cfile '/Climatology/'];
load([fpath 'LME_clim_fished_',harv,'_' cfile '.mat'],'lme_mcatch',...
    'lme_mbio','lme_area');
load([fpath 'TEeff_Climatol_All_fish03_' cfile '.mat'],'lme_te');
clme_area_km2 = lme_area * 1e-6;
clear lme_area

% POEM LME biomass in MT
plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
plme_Bmbio = lme_mbio(:,9) * 1e-6;
plme_Mmbio = (lme_mbio(:,4) + lme_mbio(:,5) + lme_mbio(:,6)) * 1e-6;
plme_Lmbio = (lme_mbio(:,7) + lme_mbio(:,8)) * 1e-6;
% MT/km2
plme_mcatch = plme_mcatch ./ clme_area_km2;
plme_Fmcatch = plme_Fmcatch ./ clme_area_km2;
plme_Pmcatch = plme_Pmcatch ./ clme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ clme_area_km2;
plme_Bmbio = plme_Bmbio ./ clme_area_km2;
plme_Mmbio = plme_Mmbio ./ clme_area_km2;
plme_Lmbio = plme_Lmbio ./ clme_area_km2;

pFracPD = plme_Pmcatch ./ (plme_Pmcatch + plme_Dmcatch);
pFracPF = plme_Pmcatch ./ (plme_Pmcatch + plme_Fmcatch);
pFracLM = plme_Lmbio ./ (plme_Lmbio + plme_Mmbio);

l10p=log10(plme_mcatch);
l10pF=log10(plme_Fmcatch);
l10pP=log10(plme_Pmcatch);
l10pD=log10(plme_Dmcatch);
l10pB=log10(plme_Bmbio);

% FEISTY LME TEeffs
l10pL = log10(lme_te(:,2));
l10pHTL = log10(lme_te(:,3));
l10pLTL = log10(lme_te(:,4));

clear lme_mcatch lme_mbio lme_te
 
%% FEISTY Clim NoNuUpdate
cfile2 = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A080_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/clim_complete/post_proc/pp_figs/',...
    cfile2,'/NoNuUpdate_'];
dpath=['/Volumes/FEISTY/NC/Clim_comp_tests/' cfile2 '/NoNuUpdate_'];
load([dpath 'LME_clim_fished_',harv,'_' cfile2 '.mat']);
load([dpath 'TEeffDet_Climatol_All_fish03_' cfile2 '.mat'],'lme_te');

hlme_mcatch = nansum(lme_mcatch,2) * 1e-6;
hlme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
hlme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
hlme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
hlme_Bmbio = lme_mbio(:,9) * 1e-6;
hlme_Mmbio = (lme_mbio(:,4) + lme_mbio(:,5) + lme_mbio(:,6)) * 1e-6;
hlme_Lmbio = (lme_mbio(:,7) + lme_mbio(:,8)) * 1e-6;

% MT/km2
hlme_area_km2 = clme_area_km2;
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
l10hHTL = log10(lme_te(:,3));
l10hLTL = log10(lme_te(:,4));

clear lme_mcatch lme_mbio lme_sbio

%% Stats
%r
[rall,pall]=corr(l10h(keep),l10p(keep));
[rF,pF]=corr(l10hF(keep),l10pF(keep));
[rP,pP]=corr(l10hP(keep),l10pP(keep));
[rD,pD]=corr(l10hD(keep),l10pD(keep));
[rB,pB]=corr(l10hB(keep),l10pB(keep));
[rPD,pPD]=corr(hFracPD(keep),pFracPD(keep));
[rPF,pPF]=corr(hFracPF(keep),pFracPF(keep));
[rLM,pLM]=corr(hFracLM(keep),pFracLM(keep));
[rL,pL]=corr(l10hL(keep),l10pL(keep));
[rHTL,pHTL]=corr(l10hHTL(keep),l10pHTL(keep));
[rLTL,pLTL]=corr(l10hLTL(keep),l10pLTL(keep));

%root mean square error
o=l10h(keep);
p=l10p(keep);
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);

o=l10hF(keep);
p=l10pF(keep);
n = length(o);
num=nansum((p-o).^2);
rmseF = sqrt(num/n);

o=l10hP(keep);
p=l10pP(keep);
n = length(o);
num=nansum((p-o).^2);
rmseP = sqrt(num/n);

o=l10hD(keep);
p=l10pD(keep);
n = length(o);
num=nansum((p-o).^2);
rmseD = sqrt(num/n);

o=l10hB(keep);
p=l10pB(keep);
n = length(o);
num=nansum((p-o).^2);
rmseB = sqrt(num/n);

o=hFracPD(keep);
p=pFracPD(keep);
n = length(o);
num=nansum((p-o).^2);
rmsePD = sqrt(num/n);

o=hFracPF(keep);
p=pFracPF(keep);
n = length(o);
num=nansum((p-o).^2);
rmsePF = sqrt(num/n);

o=hFracLM(keep);
p=pFracLM(keep);
n = length(o);
num=nansum((p-o).^2);
rmseLM = sqrt(num/n);

o=l10hL(keep);
p=l10pL(keep);
n = length(o);
num=nansum((p-o).^2);
rmseL = sqrt(num/n);

o=l10hHTL(keep);
p=l10pHTL(keep);
n = length(o);
num=nansum((p-o).^2);
rmseHTL = sqrt(num/n);

o=l10hLTL(keep);
p=l10pLTL(keep);
n = length(o);
num=nansum((p-o).^2);
rmseLTL = sqrt(num/n);

% Bias (Historic minus Climatol)
%average error = bias
%skill(3,c) = nansum(p-o) / n;
o=l10h(keep);
p=l10p(keep);
n = length(o);
bias = nansum(o-p) / n;

o=l10hF(keep);
p=l10pF(keep);
n = length(o);
biasF = nansum(o-p) / n;

o=l10hP(keep);
p=l10pP(keep);
n = length(o);
biasP = nansum(o-p) / n;

o=l10hD(keep);
p=l10pD(keep);
n = length(o);
biasD = nansum(o-p) / n;

o=l10hB(keep);
p=l10pB(keep);
n = length(o);
biasB = nansum(o-p) / n;

o=hFracPD(keep);
p=pFracPD(keep);
n = length(o);
biasPD = nansum(o-p) / n;

o=hFracPF(keep);
p=pFracPF(keep);
n = length(o);
biasPF = nansum(o-p) / n;

o=hFracLM(keep);
p=pFracLM(keep);
n = length(o);
biasLM = nansum(o-p) / n;

o=l10hL(keep);
p=l10pL(keep);
n = length(o);
biasL = nansum(o-p) / n;

o=l10hHTL(keep);
p=l10pHTL(keep);
n = length(o);
biasHTL = nansum(o-p) / n;

o=l10hLTL(keep);
p=l10pLTL(keep);
n = length(o);
biasLTL = nansum(o-p) / n;

%% Table
fish_stat(1,1) = rall;
fish_stat(2,1) = rF;
fish_stat(3,1) = rP;
fish_stat(4,1) = rD;
fish_stat(5,1) = rB;
fish_stat(6,1) = rPD;
fish_stat(7,1) = rPF;
fish_stat(8,1) = rLM;
fish_stat(9,1) = rL;
fish_stat(10,1) = rHTL;
fish_stat(11,1) = rLTL;

fish_stat(1,2) = rmse;
fish_stat(2,2) = rmseF;
fish_stat(3,2) = rmseP;
fish_stat(4,2) = rmseD;
fish_stat(5,2) = rmseB;
fish_stat(6,2) = rmsePD;
fish_stat(7,2) = rmsePF;
fish_stat(8,2) = rmseLM;
fish_stat(9,2) = rmseL;
fish_stat(10,2) = rmseHTL;
fish_stat(11,2) = rmseLTL;

fish_stat(1,3) = bias;
fish_stat(2,3) = biasF;
fish_stat(3,3) = biasP;
fish_stat(4,3) = biasD;
fish_stat(5,3) = biasB;
fish_stat(6,3) = biasPD;
fish_stat(7,3) = biasPF;
fish_stat(8,3) = biasLM;
fish_stat(9,3) = biasL;
fish_stat(10,3) = biasHTL;
fish_stat(11,3) = biasLTL;

fish_stat(1,4) = pall;
fish_stat(2,4) = pF;
fish_stat(3,4) = pP;
fish_stat(4,4) = pD;
fish_stat(5,4) = pB;
fish_stat(6,4) = pPD;
fish_stat(7,4) = pPF;
fish_stat(8,4) = pLM;
fish_stat(9,4) = pL;
fish_stat(10,4) = pHTL;
fish_stat(11,4) = pLTL;

% save
Fstat = array2table(fish_stat,'VariableNames',{'r','RMSE','Bias','p'},...
    'RowNames',{'All Fish','F','P','D','B','Frac Pel-Dem','Frac Pel-Forage',...
    'Frac Large-Med','TEeffL','TEeffHTL','TEeffLTL'});
writetable(Fstat,[dpath 'LME_clim_stats_' cfile2 '.csv'],'Delimiter',',',...
    'WriteRowNames',true)
save([dpath 'LME_clim_stats_' cfile2 '.mat'],'fish_stat')

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
scatter(l10hF(keep),l10pF(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.375 0.5 0.3 0.025],'orientation','horizontal')
text(-2.75,0.75,['r = ' sprintf('%2.2f',rF) ' (p = ' sprintf('%2.2f',pF) ')'])
text(-2.75,0.5,['RMSE = ' sprintf('%2.2f',rmseF)])
axis([-3 1 -3 1])
xlabel('NoNuUpdate')
ylabel('Clim')
title('Forage Fishes')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hP(keep),l10pP(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-5.5,1.0,['r = ' sprintf('%2.2f',rP) ' (p = ' sprintf('%2.2f',pP) ')'])
text(-5.5,0.5,['RMSE = ' sprintf('%2.2f',rmseP)])
axis([-6 2 -6 2])
xlabel('NoNuUpdate')
ylabel('Clim')
title('Large Pelagics')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hD(keep),l10pD(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-1.75,1.4,['r = ' sprintf('%2.2f',rD) ' (p = ' sprintf('%2.2f',pD) ')'])
text(-1.75,1.1,['RMSE = ' sprintf('%2.2f',rmseD)])
axis([-2 2 -2 2])
xlabel('NoNuUpdate')
ylabel('Clim')
title('Demersals')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10h(keep),l10p(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-1.75,1.4,['r = ' sprintf('%2.2f',rall) ' (p = ' sprintf('%2.2f',pall) ')'])
text(-1.75,1.1,['RMSE = ' sprintf('%2.2f',rmse)])
axis([-2 2 -2 2])
xlabel('NoNuUpdate')
ylabel('Clim')
title('All fishes')
stamp([harv '_' cfile2])
print('-dpng',[ppath 'Clim_',harv,'_comp_types_temp.png'])

%% Fractions
figure(2)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(hFracPF(keep),pFracPF(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.375 0.5 0.3 0.025],'orientation','horizontal')
text(-0.1,1.1,['r = ' sprintf('%2.2f',rPF) ' (p = ' sprintf('%2.2f',pPF) ')'])
text(-0.1,1.0,['RMSE = ' sprintf('%2.2f',rmsePF)])
axis([-0.2 1.2 -0.2 1.2])
xlabel('NoNuUpdate')
ylabel('Clim')
title('P / (P+F)')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(hFracPD(keep),pFracPD(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-0.1,1.1,['r = ' sprintf('%2.2f',rPD) ' (p = ' sprintf('%2.2f',pPD) ')'])
text(-0.1,1.0,['RMSE = ' sprintf('%2.2f',rmsePD)])
axis([-0.2 1.2 -0.2 1.2])
xlabel('NoNuUpdate')
ylabel('Clim')
title('P / (P+D)')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(hFracLM(keep),pFracLM(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-0.1,1.1,['r = ' sprintf('%2.2f',rLM) ' (p = ' sprintf('%2.2f',pLM) ')'])
text(-0.1,1.0,['RMSE = ' sprintf('%2.2f',rmseLM)])
axis([-0.2 1.2 -0.2 1.2])
xlabel('NoNuUpdate')
ylabel('Clim')
title('L / (L+M)')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hB(keep),l10pB(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-3.75,-0.3,['r = ' sprintf('%2.2f',rB) ' (p = ' sprintf('%2.2f',pB) ')'])
text(-3.75,-0.7,['RMSE = ' sprintf('%2.2f',rmseB)])
axis([-4 0 -4 0])
xlabel('NoNuUpdate')
ylabel('Clim')
title('Benthos')
stamp([harv '_' cfile2])
print('-dpng',[ppath 'Clim_',harv,'_comp_fracs_temp.png'])

% benthos figs look the same scale, so mistake somewhere else

%% TEeffs
figure(3)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hL(keep),l10pL(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.375 0.5 0.3 0.025],'orientation','horizontal')
text(-4.75,-1.25,['r = ' sprintf('%2.2f',rL) ' (p = ' sprintf('%2.2f',pL) ')'])
text(-4.75,-1.5,['RMSE = ' sprintf('%2.2f',rmseL)])
axis([-5 -1 -5 -1])
xlabel('NoNuUpdate')
ylabel('Clim')
title('log_1_0 TEeff ATL')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hHTL(keep),l10pHTL(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-3.75,-1.2,['r = ' sprintf('%2.2f',rHTL) ' (p = ' sprintf('%2.2f',pHTL) ')'])
text(-3.75,-1.4,['RMSE = ' sprintf('%2.2f',rmseHTL)])
axis([-4 -1 -4 -1])
xlabel('NoNuUpdate')
ylabel('Clim')
title('log_1_0 TEeff HTL')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hLTL(keep),l10pLTL(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-3.4,-0.65,['r = ' sprintf('%2.2f',rLTL) ' (p = ' sprintf('%2.2f',pLTL) ')'])
text(-3.4,-0.85,['RMSE = ' sprintf('%2.2f',rmseLTL)])
axis([-3.5 -0.5 -3.5 -0.5])
xlabel('NoNuUpdate')
ylabel('Clim')
title('log_1_0 TEeff LTL')
stamp([harv '_' cfile2])
print('-dpng',[ppath 'Clim_',harv,'_comp_TEeffs_temp.png'])

