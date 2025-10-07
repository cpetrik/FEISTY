% Output at all locations

clear 
close all

%%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

%fpath=['/Volumes/petrik-lab/Feisty/NC/spawning/' cfile '/CMIP6/'];
fpath=['/project/Feisty/NC/spawning/' cfile '/CMIP6/'];

nt=12*200;
vers = 'Spinup1950_const_spawning_All_fish03';

%% SP
ncid = netcdf.open([fpath vers '_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.bio = biomass;
clear biomass 

%% SF
ncid = netcdf.open([fpath vers '_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
clear biomass 

%% SD
ncid = netcdf.open([fpath vers '_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
clear biomass 

%% MP
ncid = netcdf.open([fpath vers '_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
clear biomass 

%% MF
ncid = netcdf.open([fpath vers '_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
clear biomass 

%% MD
ncid = netcdf.open([fpath vers '_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
clear biomass 

%% LP
ncid = netcdf.open([fpath vers '_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
clear biomass 

%% LD
ncid = netcdf.open([fpath vers '_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
clear biomass 

%% Benthic material
ncid = netcdf.open([fpath vers '_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

B.bio = biomass;
clear biomass

%% Save last month for initializing forecast runs
% Save last year for initializing forecast runs   .bio(:,nt-11:nt)

Sml_f.bio = mean(SF.bio(:,end),2,'omitnan');
Sml_p.bio = mean(SP.bio(:,end),2,'omitnan');
Sml_d.bio = mean(SD.bio(:,end),2,'omitnan');
Med_f.bio = mean(MF.bio(:,end),2,'omitnan');
Med_p.bio = mean(MP.bio(:,end),2,'omitnan');
Med_d.bio = mean(MD.bio(:,end),2,'omitnan');
Lrg_p.bio = mean(LP.bio(:,end),2,'omitnan');
Lrg_d.bio = mean(LD.bio(:,end),2,'omitnan');
BENT.bio  = mean(B.bio(:,end),2,'omitnan');

save([fpath 'Last_mo_' vers '_' cfile '.mat'],'Sml_f','Sml_p','Sml_d',...
    'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')

%% Take means

%Time
sp_tmean=mean(SP.bio,1,'omitnan');
sf_tmean=mean(SF.bio,1,'omitnan');
sd_tmean=mean(SD.bio,1,'omitnan');
mp_tmean=mean(MP.bio,1,'omitnan');
mf_tmean=mean(MF.bio,1,'omitnan');
md_tmean=mean(MD.bio,1,'omitnan');
lp_tmean=mean(LP.bio,1,'omitnan');
ld_tmean=mean(LD.bio,1,'omitnan');
b_tmean=mean(B.bio,1,'omitnan');

% Last year
[nid,nt] = size(LD.bio);
time=1:nt;
lyr=time((end-12+1):end);
sp_mean=mean(SP.bio(:,lyr),2,'omitnan');
sf_mean=mean(SF.bio(:,lyr),2,'omitnan');
sd_mean=mean(SD.bio(:,lyr),2,'omitnan');
mp_mean=mean(MP.bio(:,lyr),2,'omitnan');
mf_mean=mean(MF.bio(:,lyr),2,'omitnan');
md_mean=mean(MD.bio(:,lyr),2,'omitnan');
lp_mean=mean(LP.bio(:,lyr),2,'omitnan');
ld_mean=mean(LD.bio(:,lyr),2,'omitnan');
b_mean=mean(B.bio(:,lyr),2,'omitnan');


%%
save([fpath 'Means_' vers '_' cfile '.mat'],...
    'sf_mean','sp_mean','sd_mean','mf_mean','mp_mean','md_mean','b_mean',...
    'lp_mean','ld_mean','sf_tmean','sp_tmean','sd_tmean','mf_tmean','mp_tmean',...
    'md_tmean','b_tmean','lp_tmean','ld_tmean','time','lyr');


%%
Sml_f.bio = mean(SF.bio(:,lyr),2,'omitnan');
Sml_p.bio = mean(SP.bio(:,lyr),2,'omitnan');
Sml_d.bio = mean(SD.bio(:,lyr),2,'omitnan');
Med_f.bio = mean(MF.bio(:,lyr),2,'omitnan');
Med_p.bio = mean(MP.bio(:,lyr),2,'omitnan');
Med_d.bio = mean(MD.bio(:,lyr),2,'omitnan');
Lrg_p.bio = mean(LP.bio(:,lyr),2,'omitnan');
Lrg_d.bio = mean(LD.bio(:,lyr),2,'omitnan');
BENT.bio  = mean(B.bio(:,lyr),2,'omitnan');

save([fpath 'Last_yr_' vers '_' cfile '.mat'],'Sml_f','Sml_p','Sml_d',...
    'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')


