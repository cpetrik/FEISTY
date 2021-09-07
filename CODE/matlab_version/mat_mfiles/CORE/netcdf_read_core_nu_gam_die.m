% POEM output at all locations

clear all
close all


cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
vars = '_nu_gam_die';

fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/'];

%% SP
ncid = netcdf.open([fpath 'Core_' harv vars '_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.nu = nu;
SP.gamma = gamma;
SP.die = die;
clear clev gamma die

% SF
ncid = netcdf.open([fpath 'Core_' harv vars '_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.nu = nu;
SF.gamma = gamma;
SF.die = die;
clear clev gamma die

% SD
ncid = netcdf.open([fpath 'Core_' harv vars '_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.nu = nu;
SD.gamma = gamma;
SD.die = die;
clear clev gamma die

% MP
ncid = netcdf.open([fpath 'Core_' harv vars '_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.nu = nu;
MP.gamma = gamma;
MP.die = die;
clear clev gamma die

% MF
ncid = netcdf.open([fpath 'Core_' harv vars '_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.nu = nu;
MF.gamma = gamma;
MF.die = die;
clear clev gamma die nu

% MD
ncid = netcdf.open([fpath 'Core_' harv vars '_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.nu = nu;
MD.gamma = gamma;
MD.die = die;
clear clev gamma die

% LP
ncid = netcdf.open([fpath 'Core_' harv vars '_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.nu = nu;
LP.gamma = gamma;
LP.die = die;
clear clev gamma die nu

% LD
ncid = netcdf.open([fpath 'Core_' harv vars '_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.nu = nu;
LD.gamma = gamma;
LD.die = die;
clear clev gamma die nu


%% Take means
[ids,nt] = size(LD.die);

%Time
sp_tmgamma=mean(SP.gamma,1);
sf_tmgamma=mean(SF.gamma,1);
sd_tmgamma=mean(SD.gamma,1);
mp_tmgamma=mean(MP.gamma,1);
mf_tmgamma=mean(MF.gamma,1);
md_tmgamma=mean(MD.gamma,1);
lp_tmgamma=mean(LP.gamma,1);
ld_tmgamma=mean(LD.gamma,1);

sp_tmdie=mean(SP.die,1);
sf_tmdie=mean(SF.die,1);
sd_tmdie=mean(SD.die,1);
mp_tmdie=mean(MP.die,1);
mf_tmdie=mean(MF.die,1);
md_tmdie=mean(MD.die,1);
lp_tmdie=mean(LP.die,1);
ld_tmdie=mean(LD.die,1);

sp_tmnu=mean(SP.nu,1);
sf_tmnu=mean(SF.nu,1);
sd_tmnu=mean(SD.nu,1);
mp_tmnu=mean(MP.nu,1);
mf_tmnu=mean(MF.nu,1);
md_tmnu=mean(MD.nu,1);
lp_tmnu=mean(LP.nu,1);
ld_tmnu=mean(LD.nu,1);

%% Spatial mean over all time
time=1:nt;
lyr=time((end-12+1):end);

sp_mgamma=mean(SP.gamma,2);
sf_mgamma=mean(SF.gamma,2);
sd_mgamma=mean(SD.gamma,2);
mp_mgamma=mean(MP.gamma,2);
mf_mgamma=mean(MF.gamma,2);
md_mgamma=mean(MD.gamma,2);
lp_mgamma=mean(LP.gamma,2);
ld_mgamma=mean(LD.gamma,2);

sp_mdie=mean(SP.die,2);
sf_mdie=mean(SF.die,2);
sd_mdie=mean(SD.die,2);
mp_mdie=mean(MP.die,2);
mf_mdie=mean(MF.die,2);
md_mdie=mean(MD.die,2);
lp_mdie=mean(LP.die,2);
ld_mdie=mean(LD.die,2);

sp_mnu=mean(SP.nu,2);
sf_mnu=mean(SF.nu,2);
sd_mnu=mean(SD.nu,2);
mp_mnu=mean(MP.nu,2);
mf_mnu=mean(MF.nu,2);
md_mnu=mean(MD.nu,2);
lp_mnu=mean(LP.nu,2);
ld_mnu=mean(LD.nu,2);

%% Every mo
sp_die=SP.die;
sf_die=SF.die;
sd_die=SD.die;
mp_die=MP.die;
mf_die=MF.die;
md_die=MD.die;
lp_die=LP.die;
ld_die=LD.die;

sp_gamma=SP.gamma;
sf_gamma=SF.gamma;
sd_gamma=SD.gamma;
mp_gamma=MP.gamma;
mf_gamma=MF.gamma;
md_gamma=MD.gamma;
lp_gamma=LP.gamma;
ld_gamma=LD.gamma;

sp_nu=SP.nu;
sf_nu=SF.nu;
sd_nu=SD.nu;
mp_nu=MP.nu;
mf_nu=MF.nu;
md_nu=MD.nu;
lp_nu=LP.nu;
ld_nu=LD.nu;

%%

save([fpath 'Means_nu_gam_die_Core_' harv '_' cfile '.mat'],...
    'sf_tmgamma','sp_tmgamma','sd_tmgamma','mf_tmgamma','mp_tmgamma',...
    'md_tmgamma','lp_tmgamma','ld_tmgamma',...
    'sf_tmdie','sp_tmdie','sd_tmdie','mf_tmdie','mp_tmdie',...
    'md_tmdie','lp_tmdie','ld_tmdie',...
    'mf_tmnu','lp_tmnu','ld_tmnu',...
    'sf_mgamma','sp_mgamma','sd_mgamma','mf_mgamma','mp_mgamma',...
    'md_mgamma','lp_mgamma','ld_mgamma',...
    'sf_mdie','sp_mdie','sd_mdie','mf_mdie','mp_mdie',...
    'md_mdie','lp_mdie','ld_mdie',...
    'sf_mnu','sp_mnu','sd_mnu','mf_mnu','mp_mnu','md_mnu',...
    'lp_mnu','ld_mnu',...
    'time','lyr',...
    'sf_die','sp_die','sd_die','mf_die','mp_die','md_die',...
    'lp_die','ld_die','sf_gamma','sp_gamma','sd_gamma',...
    'mf_gamma','mp_gamma','md_gamma','lp_gamma','ld_gamma',...
    'sf_nu','sp_nu','sd_nu',...
    'mf_nu','mp_nu','md_nu','lp_nu','ld_nu');







