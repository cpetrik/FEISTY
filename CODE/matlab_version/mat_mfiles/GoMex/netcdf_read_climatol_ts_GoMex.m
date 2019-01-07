% POEM output at all locations

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);

%% GoMex cells
[ni,nj]=size(lon);
ocean = lmask;
ocean(ID) = ID;
lid = find(lme_mask_onedeg==5);
[G,ia,ib] = intersect(ID,lid);
%ia are vector rows for GoMex

%% SP
ncid = netcdf.open([fpath 'Climatol_' harv '_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.bio = biomass;
SP.prod = prod;
clear biomass prod

% SF
ncid = netcdf.open([fpath 'Climatol_' harv '_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass;
SF.prod = prod;
clear biomass prod

% SD
ncid = netcdf.open([fpath 'Climatol_' harv '_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
SD.prod = prod;
clear biomass prod

% MP
ncid = netcdf.open([fpath 'Climatol_' harv '_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
MP.prod = prod;
MP.yield = yield;
clear biomass prod yield

%% MF
ncid = netcdf.open([fpath 'Climatol_' harv '_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
MF.prod = prod;
MF.yield = yield;
clear biomass prod yield

% MD
ncid = netcdf.open([fpath 'Climatol_' harv '_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
MD.prod = prod;
MD.yield = yield;
clear biomass prod yield

% LP
ncid = netcdf.open([fpath 'Climatol_' harv '_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
LP.prod = prod;
LP.yield = yield;
clear biomass prod yield

% LD
ncid = netcdf.open([fpath 'Climatol_' harv '_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
LD.prod = prod;
LD.yield = yield;
clear biomass prod yield

% Benthic material
ncid = netcdf.open([fpath 'Climatol_' harv '_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

BENT.bio = biomass;
clear biomass

%% Take means and totals
nt = length(time);
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];


%Time
sp_tmean=mean(SP.bio(ia,1:nt),1);
sf_tmean=mean(SF.bio(ia,1:nt),1);
sd_tmean=mean(SD.bio(ia,1:nt),1);
mp_tmean=mean(MP.bio(ia,1:nt),1);
mf_tmean=mean(MF.bio(ia,1:nt),1);
md_tmean=mean(MD.bio(ia,1:nt),1);
lp_tmean=mean(LP.bio(ia,1:nt),1);
ld_tmean=mean(LD.bio(ia,1:nt),1);
b_tmean=mean(BENT.bio(ia,1:nt),1);

sp_tmprod=mean(SP.prod(ia,1:nt),1);
sf_tmprod=mean(SF.prod(ia,1:nt),1);
sd_tmprod=mean(SD.prod(ia,1:nt),1);
mp_tmprod=mean(MP.prod(ia,1:nt),1);
mf_tmprod=mean(MF.prod(ia,1:nt),1);
md_tmprod=mean(MD.prod(ia,1:nt),1);
lp_tmprod=mean(LP.prod(ia,1:nt),1);
ld_tmprod=mean(LD.prod(ia,1:nt),1);

mf_tmy=mean(MF.yield(ia,1:nt),1);
mp_tmy=mean(MP.yield(ia,1:nt),1);
md_tmy=mean(MD.yield(ia,1:nt),1);
lp_tmy=mean(LP.yield(ia,1:nt),1);
ld_tmy=mean(LD.yield(ia,1:nt),1);


%%
save([fpath 'GoMex_ts_means_bio_prod_fish_Climatol_' harv '_' cfile '.mat'],...
    'time','ia',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean','b_tmean',...
    'sf_tmprod','sp_tmprod','sd_tmprod',...
    'mf_tmprod','mp_tmprod','md_tmprod',...
    'lp_tmprod','ld_tmprod',...
    'mf_tmy','mp_tmy','md_tmy','lp_tmy','ld_tmy');





