% All correlations
clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

harv = 'All_fish03';

%% Clim grid
Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
cdir='/Volumes/FEISTY/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([gpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([gpath 'esm26_area_1deg.mat']);
load([gpath 'LME_clim_temp_zoop_det.mat']);

%% FEISTY No Nu Update
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A080_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/clim_complete/post_proc/pp_figs/',...
    cfile,'/NoNuUpdate_'];
fpath=['/Volumes/FEISTY/NC/Clim_comp_tests/' cfile '/NoNuUpdate_'];
load([fpath 'LME_clim_fished_',harv,'_' cfile '.mat']);
load([fpath 'TEeffDet_Climatol_All_fish03_' cfile '.mat'],'lme_te');

%%
hlme_ptemp = lme_ptemp;
hlme_area_km2 = lme_area * 1e-6;
hlme = lme_mask_onedeg;
hID = ID;

%%
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

pFracPD = hlme_Pmcatch ./ (hlme_Pmcatch + hlme_Dmcatch);
pFracPF = hlme_Pmcatch ./ (hlme_Pmcatch + hlme_Fmcatch);
pFracLM = hlme_Lmbio ./ (hlme_Lmbio + hlme_Mmbio);

l10p=log10(hlme_mcatch);
l10pF=log10(hlme_Fmcatch);
l10pP=log10(hlme_Pmcatch);
l10pD=log10(hlme_Dmcatch);
l10pB=log10(hlme_Bmbio);

l10ATL = log10(lme_te(:,2));
l10HTL = log10(lme_te(:,3));
l10LTL = log10(lme_te(:,4));

%% SAUP: All, F, P, D, Frac P
load([spath 'SAUP_LME_Catch_top10_Stock.mat']);
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

%r
[rall,pall]=corr(l10s(keep),l10p(keep));
[rF,pF]=corr(l10sF(keep),l10pF(keep));
[rP,pP]=corr(l10sP(keep),l10pP(keep));
[rD,pD]=corr(l10sD(keep),l10pD(keep));
[rPD,pPD]=corr(sFracPD(keep),pFracPD(keep));

%root mean square error
o=l10s(keep);
p=l10p(keep);
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);

o=l10sF(keep);
p=l10pF(keep);
n = length(o);
num=nansum((p-o).^2);
rmseF = sqrt(num/n);

o=l10sP(keep);
p=l10pP(keep);
n = length(o);
num=nansum((p-o).^2);
rmseP = sqrt(num/n);

o=l10sD(keep);
p=l10pD(keep);
n = length(o);
num=nansum((p-o).^2);
rmseD = sqrt(num/n);

o=sFracPD(keep);
p=pFracPD(keep);
n = length(o);
num=nansum((p-o).^2);
rmsePD = sqrt(num/n);

%Fmed
Fall=10^(median(l10s(keep)-l10p(keep)));
FF=10^(median(l10sF(keep)-l10pF(keep)));
FP=10^(median(l10sP(keep)-l10pP(keep)));
FD=10^(median(l10sD(keep)-l10pD(keep)));
FPD=10^(median(sFracPD(keep)-pFracPD(keep)));

%% DvD: Frac P
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/DanielVD_PelDem/Colleen_modeledfish_LME.mat')

did=[1:61,63];
did2 = notLELC(notLELC<=63);

%r
[rDvD,pDvD]=corr(FracLP(did),pFracPD(did));
[rDvD2,pDvD2]=corr(FracLP(did2),pFracPD(did2));

%root mean square error
o=FracLP(did);
p=pFracPD(did);
n = length(o);
num=nansum((p-o).^2);
rmseDvD = sqrt(num/n);

o=FracLP(did2);
p=pFracPD(did2);
n = length(o);
num=nansum((p-o).^2);
rmseDvD2 = sqrt(num/n);

%Fmed
FDvD=10^(median(FracLP(did)-pFracPD(did)));
FDvD2=10^(median(FracLP(did2)-pFracPD(did2)));

%% Stock: All
%Reconciling Fisheries Catch and Ocean Productivity
%***TEMPLATE FOR FEEDBACK, PENDING FINAL CHECKS***
%model: 4
%alpha: 0.14
%TE_p: 0.13
%TE_b: 0.40
%f_T: 0.74
%T_100,warm: 19.99
%All fluxes in g C m-2 yr-1, Temperature in degrees celsius
%cols = LME  NPP   MESOZP  FDET   TLeq     T  modcatch SAUcatch
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'])

% POEM in gC/m2
%gWW
wlme_mcatch = nansum(lme_mcatch,2);
%gWW/m2
wmlme_mcatch = wlme_mcatch ./ lme_area;
%gC/m2
glme_mcatch = wmlme_mcatch / 9;

%r
[rPNAS,pPNAS]=corr(StockPNAS(:,7),glme_mcatch(keep));

%root mean square error
o=StockPNAS(:,7);
p=glme_mcatch(keep);
n = length(o);
num=nansum((p-o).^2);
rmsePNAS = sqrt(num/n);

%Fmed
FPNAS=10^(median(StockPNAS(:,7)-glme_mcatch(keep)));

%% Mauread TE eff
spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/';
load([spath 'Maureaud_etal_2017_s002_ECI.mat']);

% ECI for clim years (1991-1995?)
mECI = mean(ECI(:,2:6),2);
mid = ECI(:,1);

% FEISTY LME TEeffs
pECI = lme_te(mid,3);

Lma = log10(mECI); %log10(mECI)
Lpo = log10(pECI);

%r
[rL,pL]=corr(Lma,Lpo);

%root mean square error
o=Lma;
p=Lpo;
n = length(o);
num=nansum((p-o).^2);
rmseL = sqrt(num/n);

%Fmed
FL=10^(median(Lma-Lpo));

%% Table
fish_stat(1,1) = rall;
fish_stat(2,1) = rF;
fish_stat(3,1) = rP;
fish_stat(4,1) = rD;
fish_stat(5,1) = rPD;
fish_stat(6,1) = rDvD;
fish_stat(7,1) = rPNAS;
fish_stat(1,2) = pall;
fish_stat(2,2) = pF;
fish_stat(3,2) = pP;
fish_stat(4,2) = pD;
fish_stat(5,2) = pPD;
fish_stat(6,2) = pDvD;
fish_stat(7,2) = pPNAS;
fish_stat(1,3) = rmse;
fish_stat(2,3) = rmseF;
fish_stat(3,3) = rmseP;
fish_stat(4,3) = rmseD;
fish_stat(5,3) = rmsePD;
fish_stat(6,3) = rmseDvD;
fish_stat(7,3) = rmsePNAS;

% save
Fstat = array2table(fish_stat,'VariableNames',{'r','p','RMSE'},...
    'RowNames',{'SAU All Fish','SAU F','SAU P','SAU D','SAU Frac Pelagic',...
    'vanD Frac Pelagic','Stock All Fish'});
writetable(Fstat,[fpath 'Clim_obs_LME_all_ms_stats_' cfile '.csv'],'Delimiter',',',...
    'WriteRowNames',true)
save([fpath 'Clim_obs_LME_all_ms_stats_' cfile '.mat'],'fish_stat')


