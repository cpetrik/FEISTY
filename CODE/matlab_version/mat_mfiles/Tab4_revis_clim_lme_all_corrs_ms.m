% All correlations
clear all
close all

mspath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/';
spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([cpath 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'esm26_1deg_shore_dist.mat']);
load([cpath 'LME_clim_temp_zoop_det.mat']);

AREA_OCN = max(area,1);
tlme = lme_mask_onedeg;

%% POEM
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath = [dp cfile '/'];
%load([fpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);
load([fpath 'LME_clim_fished_',harv,'_' cfile '.mat'],'lme_mcatch',...
    'lme_mbio','lme_sbio','lme_area');
lme_area_km2 = lme_area * 1e-6;

% POEM LME biomass in MT
plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
plme_Bmbio = lme_mbio(:,9) * 1e-6;
plme_Bsbio = lme_sbio(:,9) * 1e-6;
% MT/km2
plme_mcatch = plme_mcatch ./ lme_area_km2;
plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;
plme_Bmbio = plme_Bmbio ./ lme_area_km2;
plme_Bsbio = plme_Bsbio ./ lme_area_km2;

pFracPD = plme_Pmcatch ./ (plme_Pmcatch + plme_Dmcatch);

l10p=log10(plme_mcatch);
l10pF=log10(plme_Fmcatch);
l10pP=log10(plme_Pmcatch);
l10pD=log10(plme_Dmcatch);

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

% POEM
%tab  = MT/km2
%tab2 = gC/m2
%cols = 'LME','ZP','Det','Bent','T','Pcatch','Scatch'
load([fpath 'LME_prod_catch_SAUP_' cfile '.mat'],'tab','tab2')

%r
[rPNAS,pPNAS]=corr(StockPNAS(:,7),tab2(:,6));

%root mean square error
o=StockPNAS(:,7);
p=tab2(:,6);
n = length(o);
num=nansum((p-o).^2);
rmsePNAS = sqrt(num/n);

%Fmed
FPNAS=10^(median(StockPNAS(:,7)-tab2(:,6)));

%% Mauread TE eff
spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/';
load([spath 'Maureaud_etal_2017_s002_ECI.mat']);

% POEM file info
load([fpath 'TEeff_Climatol_All_fish03_' cfile '.mat']);

% ECI for clim years (1991-1995?)
mECI = mean(ECI(:,2:6),2);
mid = ECI(:,1);

% POEM LME TEeffs
pECI = lme_te(mid,4);

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
fish_stat(1,2) = rmse;
fish_stat(2,2) = rmseF;
fish_stat(3,2) = rmseP;
fish_stat(4,2) = rmseD;
fish_stat(5,2) = rmsePD;
fish_stat(6,2) = rmseDvD;
fish_stat(7,2) = rmsePNAS;
fish_stat(1,3) = pall;
fish_stat(2,3) = pF;
fish_stat(3,3) = pP;
fish_stat(4,3) = pD;
fish_stat(5,3) = pPD;
fish_stat(6,3) = pDvD;
fish_stat(7,3) = pPNAS;

% save
Fstat = array2table(fish_stat,'VariableNames',{'r','RMSE','p'},...
    'RowNames',{'SAU All Fish','SAU F','SAU P','SAU D','SAU Frac Pelagic',...
    'vanD Frac Pelagic','Stock All Fish'});
writetable(Fstat,[fpath 'Tab4_LME_all_ms_stats_' cfile '.csv'],'Delimiter',',',...
    'WriteRowNames',true)
writetable(Fstat,[mspath 'Tab4_LME_all_ms_stats_' cfile '.csv'],'Delimiter',',',...
    'WriteRowNames',true)
save([fpath 'Tab4_LME_all_ms_stats_' cfile '.mat'],'fish_stat')


