% All correlations
% CORE-forced
% Observed effort

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/SAUP/';

%% CORE-forced
load('/Volumes/MIP/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');
load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

ID = GRD.ID;

% ESM2M = same grid as CORE
gpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/cobalt_data/';
load([gpath 'hindcast_gridspec.mat'],'AREA_OCN');
load([gpath 'lme_mask_esm2m.mat']);
load([cpath 'LME_hist9095_temp_zoop_det.mat'],'lme_ptemp','lme_area');

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);
tlme = lme_mask_esm2m';

%% FEISTY LME biomass in MT/km2
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/CORE/'];
ppath = [pp cfile '/CORE/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

harv = 'fished_obs';

load([fpath 'LME_core_obs_fished_Catch_top10.mat'])

plme_mcatch = alme_mcatch10;
plme_Fmcatch = Flme_mcatch10;
plme_Pmcatch = Plme_mcatch10;
plme_Dmcatch = Dlme_mcatch10;

pFracPD = sFracPD;

l10p=log10(plme_mcatch);
l10pF=log10(plme_Fmcatch);
l10pP=log10(plme_Pmcatch);
l10pD=log10(plme_Dmcatch);

clear Flme_mcatch10 Plme_mcatch10 Dlme_mcatch10 sFracPD

%% SAUP: All, F, P, D, Frac P
load([spath 'SAUP_LME_Catch_top10_Stock.mat']);
load(['/Users/cpetrik/Dropbox/Princeton/FEISTY_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

%r
rall=corr(l10s(keep),l10p(keep));
rF=corr(l10sF(keep),l10pF(keep));
rP=corr(l10sP(keep),l10pP(keep));
rD=corr(l10sD(keep),l10pD(keep));
rPD=corr(sFracPD(keep),pFracPD(keep));

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
load('/Users/cpetrik/Dropbox/Princeton/FEISTY_other/DanielVD_PelDem/Colleen_modeledfish_LME.mat')

did=[1:61,63];
did2 = notLELC(notLELC<=63);

%r
rDvD=corr(FracLP(did),pFracPD(did));
rDvD2=corr(FracLP(did2),pFracPD(did2));

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
load(['/Users/cpetrik/Dropbox/Princeton/FEISTY_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'])

% FEISTY
%put back into = gC/m2
%1e6 gWW in 1 MT 
%(1/9) gC in 1 gWW
%1e6 m2 in 1 km2
% gWW    gC      m-2
% 1e6 * (1/9) * 1e-6
modcatch = plme_mcatch(StockPNAS(:,1),:) * (1/9);
test(:,1) = StockPNAS(:,7);
test(:,2) = modcatch;
 
%r
rPNAS=corr(StockPNAS(:,7),modcatch);

%root mean square error
o=StockPNAS(:,7);
p=modcatch;
n = length(o);
num=nansum((p-o).^2);
rmsePNAS = sqrt(num/n);

%Fmed
FPNAS=10^(median(StockPNAS(:,7)-modcatch));

%% Mauread TE eff
% spath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/poem_ms/';
% load([spath 'Maureaud_etal_2017_s002_ECI.mat']);
% 
% % FEISTY file info
% load([fpath 'TEeffDet_Climatol_All_fish03_' cfile '.mat']);
% 
% % ECI for clim years (1991-1995?)
% mECI = mean(ECI(:,2:6),2);
% mid = ECI(:,1);
% 
% % FEISTY LME TEeffs
% pECI = lme_te(mid,4);
% 
% Lma = log10(mECI); %log10(mECI)
% Lpo = log10(pECI);
% 
% %r
% rL=corr(Lma,Lpo);
% 
% %root mean square error
% o=Lma;
% p=Lpo;
% n = length(o);
% num=nansum((p-o).^2);
% rmseL = sqrt(num/n);
% 
% %Fmed
% FL=10^(median(Lma-Lpo));

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
fish_stat(1,3) = Fall;
fish_stat(2,3) = FF;
fish_stat(3,3) = FP;
fish_stat(4,3) = FD;
fish_stat(5,3) = FPD;
fish_stat(6,3) = FDvD;
fish_stat(7,3) = FPNAS;

% save
Fstat = array2table(fish_stat,'VariableNames',{'r','RMSE','Fmed'},...
    'RowNames',{'SAU All Fish','SAU F','SAU P','SAU D','SAU Frac Pelagic',...
    'vanD Frac Pelagic','Stock All Fish'});
writetable(Fstat,[fpath 'core_obs_fished_LME_all_ms_stats_' cfile '.csv'],'Delimiter',',',...
    'WriteRowNames',true)
save([fpath 'core_obs_fished_LME_all_ms_stats_' cfile '.mat'],'fish_stat')


