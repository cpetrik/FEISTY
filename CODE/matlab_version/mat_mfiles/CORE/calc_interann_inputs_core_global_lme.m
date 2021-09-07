% calc interan means and anomalies of inputs from CORE
% globally and by LME

clear all
close all

fpath = '/Volumes/MIP/GCM_DATA/CORE-forced/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% grid info
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t','AREA_OCN');
load([cpath 'lme_mask_esm2m.mat']);

load([fpath 'CORE_interann_mean_forcings_anom.mat']);

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);
tlme = lme_mask_esm2m';
LAT = double(geolat_t);
LON = double(geolon_t);

%% time means 
[ni,nj,nt] = size(tp_100);
tp2 = reshape(tp,ni*nj,nyr);
tb2 = reshape(tb,ni*nj,nyr);
det2 = reshape(det,ni*nj,nyr);
mz = reshape(zoo,ni*nj,nyr);
mzl = reshape(zlos,ni*nj,nyr);

%% time mean anomaly
mtp = nanmean(tp2);
mtb = nanmean(tb2);
mdet = nanmean(det2);
mzoo = nanmean(mz);
mzlos = nanmean(mzl);

amtp = mtp - nanmean(mtp);
amtb = mtb - nanmean(mtb);
amdet = mdet - nanmean(mdet);
amzoo = mzoo - nanmean(mzoo);
amzlos = mzlos - nanmean(mzlos);

%% time means in LMEs
ltp  = NaN*ones(66,nyr);
ltb  = NaN*ones(66,nyr);
ldet = NaN*ones(66,nyr);
lmz  = NaN*ones(66,nyr);
lmzl = NaN*ones(66,nyr);

for L=1:66
    lid = find(tlme==L);
    
    ltp(L,:) = nanmean(tp2(lid,:));
    ltb(L,:) = nanmean(tb2(lid,:));
    ldet(L,:) = nanmean(det2(lid,:));
    lmz(L,:) = nanmean(mz(lid,:));
    lmzl(L,:) = mean(mzl(lid,:));
    
end

%% LME anomalies
altp = ltp - nanmean(ltp,2);
altb = ltb - nanmean(ltb,2);
aldet = ldet - nanmean(ldet,2);
alzm = lmz - nanmean(lmz,2);
alzml = lmzl - nanmean(lmzl,2);

%%
save([fpath 'CORE_interann_mean_forcings_anom.mat'],...
    'mtp','mtb','mdet','mzoo','mzlos',...
    'amtp','amtb','amdet','amzoo','amzlos',...
    'ltp','ltb','ldet','lmz','lmzl',...
    'altp','altb','aldet','alzm','alzml','-append')
