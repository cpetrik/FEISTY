% POEM output at all locations

clear all
close all

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D100_nmort0_BE05_CC050_RE1000';

fpath=['/Volumes/GFDL/NC/Matlab_runs/' cfile '/'];


%% SP
ncid = netcdf.open([fpath 'Spinup_pristine_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.bio = biomass;
SP.clev = clev;
SP.con = con;
SP.die = die;
SP.gamma = gamma;
SP.nu = nu;
SP.prod = prod;
SP.rec = rec;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% SF
ncid = netcdf.open([fpath 'Spinup_pristine_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass;
SF.clev = clev;
SF.con = con;
SF.die = die;
SF.gamma = gamma;
SF.nu = nu;
SF.prod = prod;
SF.rec = rec;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% SD
ncid = netcdf.open([fpath 'Spinup_pristine_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
SD.clev = clev;
SD.con = con;
SD.die = die;
SD.gamma = gamma;
SD.nu = nu;
SD.prod = prod;
SD.rec = rec;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% MP
ncid = netcdf.open([fpath 'Spinup_pristine_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
MP.clev = clev;
MP.con = con;
MP.die = die;
MP.gamma = gamma;
MP.nu = nu;
MP.prod = prod;
MP.rec = rec;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% MF
ncid = netcdf.open([fpath 'Spinup_pristine_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
MF.clev = clev;
MF.con = con;
% MF.DD = DD;
MF.die = die;
% MF.egg = egg;
MF.gamma = gamma;
MF.nu = nu;
MF.prod = prod;
MF.rec = rec;
MF.rep = rep;
% MF.S = S;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% MD
ncid = netcdf.open([fpath 'Spinup_pristine_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
MD.clev = clev;
MD.con = con;
MD.die = die;
MD.gamma = gamma;
MD.nu = nu;
MD.prod = prod;
MD.rec = rec;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% LP
ncid = netcdf.open([fpath 'Spinup_pristine_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
LP.clev = clev;
LP.con = con;
% LP.DD = DD;
LP.die = die;
% LP.egg = egg;
LP.gamma = gamma;
LP.nu = nu;
LP.prod = prod;
LP.rec = rec;
LP.rep = rep;
% LP.S = S;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% LD
ncid = netcdf.open([fpath 'Spinup_pristine_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
LD.clev = clev;
LD.con = con;
% LD.DD = DD;
LD.die = die;
% LD.egg = egg;
LD.gamma = gamma;
LD.nu = nu;
LD.prod = prod;
LD.rec = rec;
LD.rep = rep;
% LD.S = S;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% Benthic material
ncid = netcdf.open([fpath 'Spinup_pristine_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

BENT.bio = biomass;
clear biomass 

%% Take means

%Time
sp_tmean=mean(SP.bio,1);
sf_tmean=mean(SF.bio,1);
sd_tmean=mean(SD.bio,1);
mp_tmean=mean(MP.bio,1);
mf_tmean=mean(MF.bio,1);
md_tmean=mean(MD.bio,1);
lp_tmean=mean(LP.bio,1);
ld_tmean=mean(LD.bio,1);
b_tmean=mean(BENT.bio,1);

% Last year
lyr=time((end-12+1):end);
sp_mean=mean(SP.bio(:,lyr),2);
sf_mean=mean(SF.bio(:,lyr),2);
sd_mean=mean(SD.bio(:,lyr),2);
mp_mean=mean(MP.bio(:,lyr),2);
mf_mean=mean(MF.bio(:,lyr),2);
md_mean=mean(MD.bio(:,lyr),2);
lp_mean=mean(LP.bio(:,lyr),2);
ld_mean=mean(LD.bio(:,lyr),2);
b_mean=mean(BENT.bio(:,lyr),2);


%
save([fpath 'Means_spinup_' cfile '.mat'],...
    'sf_mean','sp_mean','sd_mean','mf_mean','mp_mean','md_mean','b_mean',...
    'lp_mean','ld_mean','sf_tmean','sp_tmean','sd_tmean','mf_tmean','mp_tmean',...
    'md_tmean','b_tmean','lp_tmean','ld_tmean','time','lyr');


