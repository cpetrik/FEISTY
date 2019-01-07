% POEM output at all locations

clear all
close all


cfile = 'Dc_enc70-b200_m4-b175-k08_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
vars = '_nu_gam_die_clev';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];

%% SP
ncid = netcdf.open([fpath 'Climatol_' harv vars '_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.clev = clev;
SP.nu = nu;
SP.gamma = gamma;
SP.die = die;
clear clev gamma die

% SF
ncid = netcdf.open([fpath 'Climatol_' harv vars '_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.clev = clev;
SF.nu = nu;
SF.gamma = gamma;
SF.die = die;
clear clev gamma die

% SD
ncid = netcdf.open([fpath 'Climatol_' harv vars '_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.clev = clev;
SD.nu = nu;
SD.gamma = gamma;
SD.die = die;
clear clev gamma die

% MP
ncid = netcdf.open([fpath 'Climatol_' harv vars '_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.clev = clev;
MP.nu = nu;
MP.gamma = gamma;
MP.die = die;
clear clev gamma die

% MF
ncid = netcdf.open([fpath 'Climatol_' harv vars '_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.clev = clev;
MF.nu = nu;
MF.gamma = gamma;
MF.die = die;
clear clev gamma die nu

% MD
ncid = netcdf.open([fpath 'Climatol_' harv vars '_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.clev = clev;
MD.nu = nu;
MD.gamma = gamma;
MD.die = die;
clear clev gamma die

% LP
ncid = netcdf.open([fpath 'Climatol_' harv vars '_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.clev = clev;
LP.nu = nu;
LP.gamma = gamma;
LP.die = die;
clear clev gamma die nu

% LD
ncid = netcdf.open([fpath 'Climatol_' harv vars '_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.clev = clev;
LD.nu = nu;
LD.gamma = gamma;
LD.die = die;
clear clev gamma die nu


%% Take means
[ids,nt] = size(LD.clev);

%Time
sp_tmclev=mean(SP.clev,1);
sf_tmclev=mean(SF.clev,1);
sd_tmclev=mean(SD.clev,1);
mp_tmclev=mean(MP.clev,1);
mf_tmclev=mean(MF.clev,1);
md_tmclev=mean(MD.clev,1);
lp_tmclev=mean(LP.clev,1);
ld_tmclev=mean(LD.clev,1);

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

%% Last year
time=1:nt;
lyr=time((end-12+1):end);
sp_mclev=mean(SP.clev(:,lyr),2);
sf_mclev=mean(SF.clev(:,lyr),2);
sd_mclev=mean(SD.clev(:,lyr),2);
mp_mclev=mean(MP.clev(:,lyr),2);
mf_mclev=mean(MF.clev(:,lyr),2);
md_mclev=mean(MD.clev(:,lyr),2);
lp_mclev=mean(LP.clev(:,lyr),2);
ld_mclev=mean(LD.clev(:,lyr),2);

sp_mgamma=mean(SP.gamma(:,lyr),2);
sf_mgamma=mean(SF.gamma(:,lyr),2);
sd_mgamma=mean(SD.gamma(:,lyr),2);
mp_mgamma=mean(MP.gamma(:,lyr),2);
mf_mgamma=mean(MF.gamma(:,lyr),2);
md_mgamma=mean(MD.gamma(:,lyr),2);
lp_mgamma=mean(LP.gamma(:,lyr),2);
ld_mgamma=mean(LD.gamma(:,lyr),2);

sp_mdie=mean(SP.die(:,lyr),2);
sf_mdie=mean(SF.die(:,lyr),2);
sd_mdie=mean(SD.die(:,lyr),2);
mp_mdie=mean(MP.die(:,lyr),2);
mf_mdie=mean(MF.die(:,lyr),2);
md_mdie=mean(MD.die(:,lyr),2);
lp_mdie=mean(LP.die(:,lyr),2);
ld_mdie=mean(LD.die(:,lyr),2);

sp_mnu=mean(SP.nu(:,lyr),2);
sf_mnu=mean(SF.nu(:,lyr),2);
sd_mnu=mean(SD.nu(:,lyr),2);
mp_mnu=mean(MP.nu(:,lyr),2);
mf_mnu=mean(MF.nu(:,lyr),2);
md_mnu=mean(MD.nu(:,lyr),2);
lp_mnu=mean(LP.nu(:,lyr),2);
ld_mnu=mean(LD.nu(:,lyr),2);

%%
save([fpath 'Means_Climatol_' harv '_' cfile '.mat'],'sf_mclev',...
    'sf_tmclev','sp_tmclev','sd_tmclev','mf_tmclev','mp_tmclev',...
    'md_tmclev','lp_tmclev','ld_tmclev',...
    'sf_tmgamma','sp_tmgamma','sd_tmgamma','mf_tmgamma','mp_tmgamma',...
    'md_tmgamma','lp_tmgamma','ld_tmgamma',...
    'sf_tmdie','sp_tmdie','sd_tmdie','mf_tmdie','mp_tmdie',...
    'md_tmdie','lp_tmdie','ld_tmdie',...
    'sf_tmnu','sp_tmnu','sd_tmnu','mf_tmnu','mp_tmnu','md_tmnu',...
    'lp_tmnu','ld_tmnu',...
    'sf_mclev','sp_mclev','sd_mclev','mf_mclev','mp_mclev',...
    'md_mclev','lp_mclev','ld_mclev',...
    'sf_mgamma','sp_mgamma','sd_mgamma','mf_mgamma','mp_mgamma',...
    'md_mgamma','lp_mgamma','ld_mgamma',...
    'sf_mdie','sp_mdie','sd_mdie','mf_mdie','mp_mdie',...
    'md_mdie','lp_mdie','ld_mdie',...
    'sf_mnu','sp_mnu','sd_mnu','mf_mnu','mp_mnu','md_mnu',...
    'lp_mnu','ld_mnu',...
    'time','lyr','-append');

save([fpath 'Means_nu_gam_die_clev_Climatol_' harv '_' cfile '.mat'],'sf_mclev',...
    'sf_tmclev','sp_tmclev','sd_tmclev','mf_tmclev','mp_tmclev',...
    'md_tmclev','lp_tmclev','ld_tmclev',...
    'sf_tmgamma','sp_tmgamma','sd_tmgamma','mf_tmgamma','mp_tmgamma',...
    'md_tmgamma','lp_tmgamma','ld_tmgamma',...
    'sf_tmdie','sp_tmdie','sd_tmdie','mf_tmdie','mp_tmdie',...
    'md_tmdie','lp_tmdie','ld_tmdie',...
    'mf_tmnu','lp_tmnu','ld_tmnu',...
    'sf_mclev','sp_mclev','sd_mclev','mf_mclev','mp_mclev',...
    'md_mclev','lp_mclev','ld_mclev',...
    'sf_mgamma','sp_mgamma','sd_mgamma','mf_mgamma','mp_mgamma',...
    'md_mgamma','lp_mgamma','ld_mgamma',...
    'sf_mdie','sp_mdie','sd_mdie','mf_mdie','mp_mdie',...
    'md_mdie','lp_mdie','ld_mdie',...
    'sf_mnu','sp_mnu','sd_mnu','mf_mnu','mp_mnu','md_mnu',...
    'lp_mnu','ld_mnu',...
    'time','lyr');







