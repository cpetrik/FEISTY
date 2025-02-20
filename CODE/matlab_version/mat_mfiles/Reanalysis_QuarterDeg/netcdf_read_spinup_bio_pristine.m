% Output at all locations

clear 
close all

cfile = 'Dc_enc70-b200_m400-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/QuarterDeg/'];

nt=12*200;

%% Benthic material
ncid = netcdf.open([fpath 'Spinup_pristine_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

B.bio = biomass;

lyr=time((end-12+1):end);

b_tmean=mean(B.bio,1,'omitnan');
b_mean=mean(B.bio(:,lyr),2,'omitnan');
save([fpath 'Means_spinup_' cfile '.mat'],...
    'b_mean','b_tmean','time','lyr','-append');
%%
clear biomass B

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

% Last year
sp_tmean=mean(SP.bio,1,'omitnan');
sp_mean=mean(SP.bio(:,lyr),2,'omitnan');

clear biomass SP

%%
save([fpath 'Means_spinup_' cfile '.mat'],...
    'sp_mean','sp_tmean','-append');

%% SF
ncid = netcdf.open([fpath 'Spinup_pristine_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);

sf_tmean=mean(SF.bio,1,'omitnan');
sf_mean=mean(SF.bio(:,lyr),2,'omitnan');
save([fpath 'Means_spinup_' cfile '.mat'],...
    'sf_mean','sf_tmean','-append');

clear biomass SF

%% SD
ncid = netcdf.open([fpath 'Spinup_pristine_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;

sd_tmean=mean(SD.bio,1,'omitnan');
sd_mean=mean(SD.bio(:,lyr),2,'omitnan');
save([fpath 'Means_spinup_' cfile '.mat'],...
    'sd_mean','sd_tmean','-append');

clear biomass SD

%% MP
ncid = netcdf.open([fpath 'Spinup_pristine_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;

mp_tmean=mean(MP.bio,1,'omitnan');
mp_mean=mean(MP.bio(:,lyr),2,'omitnan');

%%
save([fpath 'Means_spinup_' cfile '.mat'],...
    'mp_mean','mp_tmean','-append');

clear biomass MP

%% MF
ncid = netcdf.open([fpath 'Spinup_pristine_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;

mf_tmean=mean(MF.bio,1,'omitnan');
mf_mean=mean(MF.bio(:,lyr),2,'omitnan');

save([fpath 'Means_spinup_' cfile '.mat'],...
    'mf_mean','mf_tmean','-append');

clear biomass MF

%% MD
ncid = netcdf.open([fpath 'Spinup_pristine_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;

md_tmean=mean(MD.bio,1,'omitnan');
md_mean=mean(MD.bio(:,lyr),2,'omitnan');
save([fpath 'Means_spinup_' cfile '.mat'],...
    'md_mean','md_tmean','-append');

clear biomass MD


%% LP
ncid = netcdf.open([fpath 'Spinup_pristine_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;

lp_tmean=mean(LP.bio,1,'omitnan');
lp_mean=mean(LP.bio(:,lyr),2,'omitnan');
save([fpath 'Means_spinup_' cfile '.mat'],...
    'lp_mean','lp_tmean','-append');

clear biomass LP

%% LD
ncid = netcdf.open([fpath 'Spinup_pristine_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;

ld_tmean=mean(LD.bio,1,'omitnan');
ld_mean=mean(LD.bio(:,lyr),2,'omitnan');
save([fpath 'Means_spinup_' cfile '.mat'],...
    'ld_mean','ld_tmean','-append');

clear biomass LD

%%





