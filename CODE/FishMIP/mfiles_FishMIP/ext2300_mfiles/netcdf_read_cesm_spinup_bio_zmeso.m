% FEISTY output at all locations

clear
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

mod = 'CESM2-WACCM';
exper = [mod '_spinup_zmeso_pristine'];

%fpath=['/Volumes/petrik-lab/Feisty/NC/WG2300/',cfile,'/',mod,'/'];
fpath=['/project/Feisty/NC/WG2300/',cfile,'/',mod,'/'];

%% SP
ncid = netcdf.open([fpath exper '_empHP_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%%
[ni,nt] = size(biomass);

SP.bio = biomass;
clear biomass

%% SF
ncid = netcdf.open([fpath exper '_empHP_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
clear biomass

% SD
ncid = netcdf.open([fpath exper '_empHP_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
clear biomass

% MP
ncid = netcdf.open([fpath exper '_empHP_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
clear biomass

% MF
ncid = netcdf.open([fpath exper '_empHP_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
clear biomass

% MD
ncid = netcdf.open([fpath exper '_empHP_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
clear biomass

% LP
ncid = netcdf.open([fpath exper '_empHP_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
clear biomass

% LD
ncid = netcdf.open([fpath exper '_empHP_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
clear biomass

% Benthic material
ncid = netcdf.open([fpath exper '_empHP_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

Bent.bio = biomass;
clear biomass

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
b_tmean=mean(Bent.bio,1,'omitnan');

%Space
t=time;
mo=t/12;
mo=mo+1850;
yrP=find(mo>1949 & mo<=1950);

sp_mean=mean(SP.bio(:,yrP),2);
sf_mean=mean(SF.bio(:,yrP),2);
sd_mean=mean(SD.bio(:,yrP),2);
mp_mean=mean(MP.bio(:,yrP),2);
mf_mean=mean(MF.bio(:,yrP),2);
md_mean=mean(MD.bio(:,yrP),2);
lp_mean=mean(LP.bio(:,yrP),2);
ld_mean=mean(LD.bio(:,yrP),2);
b_mean =mean(Bent.bio(:,yrP),2);

save([fpath 'Means_' exper '_' cfile '.mat'],'time','mo',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean','b_tmean',...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean')


%% Save last year for initializing forecast runs
Sml_f.bio = mean(SF.bio(:,nt-11:nt),2,'omitnan');
Sml_p.bio = mean(SP.bio(:,nt-11:nt),2,'omitnan');
Sml_d.bio = mean(SD.bio(:,nt-11:nt),2,'omitnan');
Med_f.bio = mean(MF.bio(:,nt-11:nt),2,'omitnan');
Med_p.bio = mean(MP.bio(:,nt-11:nt),2,'omitnan');
Med_d.bio = mean(MD.bio(:,nt-11:nt),2,'omitnan');
Lrg_p.bio = mean(LP.bio(:,nt-11:nt),2,'omitnan');
Lrg_d.bio = mean(LD.bio(:,nt-11:nt),2,'omitnan');
BENT.bio  = mean(Bent.bio(:,nt-11:nt),2,'omitnan');

save([fpath 'Last_mo_' exper '_' cfile '.mat'],'Sml_f','Sml_p','Sml_d',...
    'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')


%%
% figure
% plot(time,mf_tmean,'r'); hold on
% plot(time,lp_tmean,'b'); hold on
% plot(time,ld_tmean,'k')


