% FEISTY output at all locations
% Yield per year calcs & mean

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];

% MP
ncid = netcdf.open([fpath 'Forecast_' harv '_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.yield = yield;
clear biomass yield

% MF
ncid = netcdf.open([fpath 'Forecast_' harv '_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.yield = yield;
clear biomass yield

% MD
ncid = netcdf.open([fpath 'Forecast_' harv '_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.yield = yield;
clear biomass yield

% LP
ncid = netcdf.open([fpath 'Forecast_' harv '_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.yield = yield;
clear biomass yield

% LD
ncid = netcdf.open([fpath 'Forecast_' harv '_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.yield = yield;
clear biomass yield

% Benthic material
ncid = netcdf.open([fpath 'Forecast_' harv '_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

clear biomass

%% Take means and totals
% Totals only in lmes
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
load([cpath 'hindcast_gridspec.mat'],'AREA_OCN');
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

load([cpath 'lme_mask_esm2m.mat']);
tlme = lme_mask_esm2m';
tlme(~isnan(tlme)) = 1;
lme_grid = tlme(ID);

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);
area = AREA_OCN(ID);
area_km2 = area * 1e-6;

[ni,nt] = size(LD.yield);
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
nyr = nt/12;
mos = repmat(MNTH,ni,nyr);
mns = repmat(MNTH,ni,1);
area_mat = repmat(area_km2,1,nt);
lme_mat = repmat(lme_grid,1,nt);

%% Units
units_yield = 'g_m2_day';
units_catch = 'g_km2_mo';

MF.catch = MF.yield .*mos .*area_mat .*lme_mat;
MP.catch = MP.yield .*mos .*area_mat .*lme_mat;
MD.catch = MD.yield .*mos .*area_mat .*lme_mat;
LP.catch = LP.yield .*mos .*area_mat .*lme_mat;
LD.catch = LD.yield .*mos .*area_mat .*lme_mat;

%% Time
% %mean yield per mo
% mf_tmc=nanmean(MF.catch,1);
% mp_tmc=nanmean(MP.catch,1);
% md_tmc=nanmean(MD.catch,1);
% lp_tmc=nanmean(LP.catch,1);
% ld_tmc=nanmean(LD.catch,1);
%
% %total yield per mo
% mf_ttc=nansum(MF.catch,1);
% mp_ttc=nansum(MP.catch,1);
% md_ttc=nansum(MD.catch,1);
% lp_ttc=nansum(LP.catch,1);
% ld_ttc=nansum(LD.catch,1);
%
% %% 50 yrs (2051-2100)
% y = 2005+(1/12):(1/12):2100;
% yr50=find(y>=2051 & y<2100);
% mf_mc50=nanmean(MF.catch(:,yr50),2);
% mp_mc50=nanmean(MP.catch(:,yr50),2);
% md_mc50=nanmean(MD.catch(:,yr50),2);
% lp_mc50=nanmean(LP.catch(:,yr50),2);
% ld_mc50=nanmean(LD.catch(:,yr50),2);
% mf_tc50=nansum(MF.catch(:,yr50),2);
% mp_tc50=nansum(MP.catch(:,yr50),2);
% md_tc50=nansum(MD.catch(:,yr50),2);
% lp_tc50=nansum(LP.catch(:,yr50),2);
% ld_tc50=nansum(LD.catch(:,yr50),2);
%
% %% Every 5 years
% st=1:60:length(time);
% en=60:60:length(time);
%
% for n=1:length(st)
%
%     mp_mc(:,n)=nanmean(MP.catch(:,st(n):en(n)),2);
%     mf_mc(:,n)=nanmean(MF.catch(:,st(n):en(n)),2);
%     md_mc(:,n)=nanmean(MD.catch(:,st(n):en(n)),2);
%     lp_mc(:,n)=nanmean(LP.catch(:,st(n):en(n)),2);
%     ld_mc(:,n)=nanmean(LD.catch(:,st(n):en(n)),2);
%
%     mp_tc(:,n)=nansum(MP.catch(:,st(n):en(n)),2);
%     mf_tc(:,n)=nansum(MF.catch(:,st(n):en(n)),2);
%     md_tc(:,n)=nansum(MD.catch(:,st(n):en(n)),2);
%     lp_tc(:,n)=nansum(LP.catch(:,st(n):en(n)),2);
%     ld_tc(:,n)=nansum(LD.catch(:,st(n):en(n)),2);
% end

%% Every year
st=1:12:length(time);
en=12:12:length(time);

mf_tyc = nan*ones(ni,nyr);
mp_tyc = nan*ones(ni,nyr);
md_tyc = nan*ones(ni,nyr);
lp_tyc = nan*ones(ni,nyr);
ld_tyc = nan*ones(ni,nyr);
for n=1:length(st)
    
    mp_tyc(:,n)=nansum(MP.catch(:,st(n):en(n)),2);
    mf_tyc(:,n)=nansum(MF.catch(:,st(n):en(n)),2);
    md_tyc(:,n)=nansum(MD.catch(:,st(n):en(n)),2);
    lp_tyc(:,n)=nansum(LP.catch(:,st(n):en(n)),2);
    ld_tyc(:,n)=nansum(LD.catch(:,st(n):en(n)),2);
end

tmn = mf_tyc + mp_tyc + md_tyc + lp_tyc + ld_tyc;
stmn = sum(tmn);

mp_tsyc = nansum(mp_tyc);
mf_tsyc = nansum(mf_tyc);
md_tsyc = nansum(md_tyc);
lp_tsyc = nansum(lp_tyc);
ld_tsyc = nansum(ld_tyc);

%%
save([fpath 'Means_fore_' harv '_' cfile '.mat'],'time',...
    'mf_tyc','mp_tyc','md_tyc','lp_tyc','ld_tyc',...
    'mf_tsyc','mp_tsyc','md_tsyc','lp_tsyc','ld_tsyc',...
    'units_yield','units_catch','-append');







