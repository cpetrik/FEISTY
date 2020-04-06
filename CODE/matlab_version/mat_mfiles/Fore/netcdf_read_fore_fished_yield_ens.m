function netcdf_read_fore_fished_yield_ens(fname,simname)

% FEISTY output at all locations
% .catch = .yield .*mos .*area_mat .*lme_mat;

%% MP
ncid = netcdf.open([fname '_med_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fname '_med_f.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fname '_med_d.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fname '_lrg_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fname '_lrg_d.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fname '_bent.nc'],'NC_NOWRITE');
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

% tmn = mf_tyc + mp_tyc + md_tyc + lp_tyc + ld_tyc;
% stmn = sum(tmn);

mp_tsyc = nansum(mp_tyc);
mf_tsyc = nansum(mf_tyc);
md_tsyc = nansum(md_tyc);
lp_tsyc = nansum(lp_tyc);
ld_tsyc = nansum(ld_tyc);

%%
save([fname '_Means_' simname '.mat'],...
    'mf_tyc','mp_tyc','md_tyc','lp_tyc','ld_tyc',...
    'mf_tsyc','mp_tsyc','md_tsyc','lp_tsyc','ld_tsyc',...
    'units_yield','units_catch','-append');

end

