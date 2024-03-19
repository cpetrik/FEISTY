% FEISTY output at all locations

clear 
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';
vers = 'HAD';
harv = 'All_fishobs';
%Spinup_IPSL_All_fishobs_bent.nc

fpath=['/Volumes/petrik-lab/Feisty/NC/NEMURO/',cfile,'/',vers,'/'];

%% Benthic material
ncid = netcdf.open([fpath 'Spinup_' vers '_' harv '_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

nt = length(time);

Bent.bio = biomass;
BENT.bio = biomass(:,nt);
clear biomass

%% SP
ncid = netcdf.open([fpath 'Spinup_' vers '_' harv '_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.bio = biomass;
Sml_p.bio = biomass(:,nt);
clear biomass 

% SF
ncid = netcdf.open([fpath 'Spinup_' vers '_' harv '_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
Sml_f.bio = biomass(:,nt);
clear biomass 

% SD
ncid = netcdf.open([fpath 'Spinup_' vers '_' harv '_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
Sml_d.bio = biomass(:,nt);
clear biomass 

%% MP
ncid = netcdf.open([fpath 'Spinup_' vers '_' harv '_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
Med_p.bio = biomass(:,nt);
clear biomass yield

% MF
ncid = netcdf.open([fpath 'Spinup_' vers '_' harv '_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
Med_f.bio = biomass(:,nt);
clear biomass yield

% MD
ncid = netcdf.open([fpath 'Spinup_' vers '_' harv '_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
Med_d.bio = biomass(:,nt);
clear biomass yield

% LP
ncid = netcdf.open([fpath 'Spinup_' vers '_' harv '_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
Lrg_p.bio = biomass(:,nt);
clear biomass yield

% LD
ncid = netcdf.open([fpath 'Spinup_' vers '_' harv '_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
Lrg_d.bio = biomass(:,nt);
clear biomass yield

%% Take means
%Time
sp_tmean=mean(SP.bio,1,"omitnan");
sf_tmean=mean(SF.bio,1,"omitnan");
sd_tmean=mean(SD.bio,1,"omitnan");
mp_tmean=mean(MP.bio,1,"omitnan");
mf_tmean=mean(MF.bio,1,"omitnan");
md_tmean=mean(MD.bio,1,"omitnan");
lp_tmean=mean(LP.bio,1,"omitnan");
ld_tmean=mean(LD.bio,1,"omitnan");
b_tmean=mean(Bent.bio,1,"omitnan");


% Last year
lyr=time((end-12+1):end);
%lyr=1:12;
sp_mean=mean(SP.bio(:,lyr),2,"omitnan");
sf_mean=mean(SF.bio(:,lyr),2,"omitnan");
sd_mean=mean(SD.bio(:,lyr),2,"omitnan");
mp_mean=mean(MP.bio(:,lyr),2,"omitnan");
mf_mean=mean(MF.bio(:,lyr),2,"omitnan");
md_mean=mean(MD.bio(:,lyr),2,"omitnan");
lp_mean=mean(LP.bio(:,lyr),2,"omitnan");
ld_mean=mean(LD.bio(:,lyr),2,"omitnan");
b_mean=mean(Bent.bio(:,lyr),2,"omitnan");

%
save([fpath 'Means_Spinup_' vers '_' harv '_' cfile '.mat'],...
    'sf_mean','sp_mean','sd_mean','mf_mean','mp_mean','md_mean','b_mean',...
    'lp_mean','ld_mean','sf_tmean','sp_tmean','sd_tmean','mf_tmean','mp_tmean',...
    'md_tmean','b_tmean','lp_tmean','ld_tmean','time','lyr');

% Save last year for initializing forecast runs
save([fpath 'Last_mo_Spinup_' vers '_' harv '_' cfile '.mat'],'Sml_f','Sml_p','Sml_d',... 
    'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')







