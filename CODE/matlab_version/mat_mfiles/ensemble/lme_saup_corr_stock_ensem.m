function [r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch)

%FEISTY catch vs. SAUP catch by LME
%Use same methods as Stock et al. 2017 to reduce SAUP dataset

dp = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/Dc_D075_Sm025_nmort1_noCC_RE00100/';
load([dp 'SAUP_Stock_top10.mat']);

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
dpath = ['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;

load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

%% SAUP
% MT/km2
sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

%% FEISTY  LME biomass in MT
plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
% MT/km2
plme_mcatch = plme_mcatch ./ lme_area_km2;
plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;

pFracPD = plme_Pmcatch ./ (plme_Pmcatch + plme_Dmcatch);

l10p=log10(plme_mcatch+eps);
l10pF=log10(plme_Fmcatch+eps);
l10pP=log10(plme_Pmcatch+eps);
l10pD=log10(plme_Dmcatch+eps);

%% Drop Arctic, Antarctic, Hawaii, Australia -------------------------
% Stats
mis = nan*ones(length(keep),5);
%r
r(1)=corr(l10s(keep),l10p(keep));
r(2)=corr(l10sF(keep),l10pF(keep));
r(3)=corr(l10sP(keep),l10pP(keep));
r(4)=corr(l10sD(keep),l10pD(keep));
r(5)=corr(sFracPD(keep),pFracPD(keep));

%root mean square error
o=l10s(keep);
p=l10p(keep);
n = length(o);
num=nansum((p-o).^2);
rmse(1) = sqrt(num/n);
ss(1) = num;
mis(:,1) = (p-o);

o=l10sF(keep);
p=l10pF(keep);
n = length(o);
num=nansum((p-o).^2);
rmse(2) = sqrt(num/n);
ss(2) = num;
mis(:,2) = (p-o);

o=l10sP(keep);
p=l10pP(keep);
n = length(o);
num=nansum((p-o).^2);
rmse(3) = sqrt(num/n);
ss(3) = num;
mis(:,3) = (p-o);

o=l10sD(keep);
p=l10pD(keep);
n = length(o);
num=nansum((p-o).^2);
rmse(4) = sqrt(num/n);
ss(4) = num;
mis(:,4) = (p-o);

o=sFracPD(keep);
p=pFracPD(keep);
n = length(o);
num=nansum((p-o).^2);
rmse(5) = sqrt(num/n);
ss(5) = num;
mis(:,5) = (p-o);

end
