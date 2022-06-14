% CORE-forced
% Observed effort
%Use same methods as Stock et al. 2017 to calc top 10 catch yrs

clear all
close all

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/CORE/'];
ppath = [pp cfile '/CORE/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

harv = 'fished_obs';
tharv = 'Observed effort';

load([fpath 'LME_core_',harv,'_' cfile '.mat']);

%% 
lme_area_km2 = lme_area * 1e-6;

%mcatch = mf, mp, md, lp, ld
%MT/km2
Alme_catch_all = nansum(lme_mcatch,3) * 1e-6 ./ lme_area_km2;
Flme_catch_all = (lme_mcatch(:,:,1)) * 1e-6 ./ lme_area_km2;
Plme_catch_all = (lme_mcatch(:,:,2)+lme_mcatch(:,:,4)) * 1e-6 ./ lme_area_km2;
Dlme_catch_all = (lme_mcatch(:,:,3)+lme_mcatch(:,:,5)) * 1e-6 ./ lme_area_km2;

alme_mcatch10 = NaN*ones(66,1);
Flme_mcatch10 = NaN*ones(66,1);
Plme_mcatch10 = NaN*ones(66,1);
Dlme_mcatch10 = NaN*ones(66,1);

%Top 10 yrs by LME SAUP
for i=1:66
    [sort_Alme_catch,ix] = sort(Alme_catch_all(i,:),'descend');
    sort_Flme_catch = Flme_catch_all(i,ix);
    sort_Plme_catch = Plme_catch_all(i,ix);
    sort_Dlme_catch = Dlme_catch_all(i,ix);
    alme_mcatch10(i) = nanmean(sort_Alme_catch(1:10));
    Flme_mcatch10(i) = nanmean(sort_Flme_catch(1:10));
    Plme_mcatch10(i) = nanmean(sort_Plme_catch(1:10));
    Dlme_mcatch10(i) = nanmean(sort_Dlme_catch(1:10));
end

%% log10
l10s=log10(alme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

%% save
save([fpath 'LME_core_obs_fished_Catch_top10.mat'],'alme_mcatch10','Flme_mcatch10',...
    'Plme_mcatch10','Dlme_mcatch10','sFracPD')


