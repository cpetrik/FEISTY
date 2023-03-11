% FEISTY output at all locations

clear 
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];

mod = 'v3.2';

%% SP
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_',mod,'_empHP_sml_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_',mod,'_empHP_sml_f.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_',mod,'_empHP_sml_d.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_',mod,'_empHP_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
MP.yield = yield;
clear yield
clear biomass

% MF
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_',mod,'_empHP_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
MF.yield = yield;
clear yield
clear biomass

% MD
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_',mod,'_empHP_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
MD.yield = yield;
clear yield
clear biomass

% LP 
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_',mod,'_empHP_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
LP.yield = yield;
clear yield
clear biomass

%% LD
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_',mod,'_empHP_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
LD.yield = yield;
clear yield
clear biomass

% Benthic material
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_',mod,'_empHP_bent.nc'],'NC_NOWRITE');
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

mp_tmcatch=mean(MP.yield,1,'omitnan');
mf_tmcatch=mean(MF.yield,1,'omitnan');
md_tmcatch=mean(MD.yield,1,'omitnan');
lp_tmcatch=mean(LP.yield,1,'omitnan');
ld_tmcatch=mean(LD.yield,1,'omitnan');

%Space
sp_smean=mean(SP.bio,2,'omitnan');
sf_smean=mean(SF.bio,2,'omitnan');
sd_smean=mean(SD.bio,2,'omitnan');
mp_smean=mean(MP.bio,2,'omitnan');
mf_smean=mean(MF.bio,2,'omitnan');
md_smean=mean(MD.bio,2,'omitnan');
lp_smean=mean(LP.bio,2,'omitnan');
ld_smean=mean(LD.bio,2,'omitnan');
b_smean =mean(Bent.bio,2,'omitnan');

mp_mcatch=mean(MP.yield,2,'omitnan');
mf_mcatch=mean(MF.yield,2,'omitnan');
md_mcatch=mean(MD.yield,2,'omitnan');
lp_mcatch=mean(LP.yield,2,'omitnan');
ld_mcatch=mean(LD.yield,2,'omitnan');

% Each year
st=1:12:length(time);
en=12:12:length(time);

for n=1:length(st)
    sp_mean(:,n)=nanmean(SP.bio(:,st(n):en(n)),2);
    sf_mean(:,n)=nanmean(SF.bio(:,st(n):en(n)),2);
    sd_mean(:,n)=nanmean(SD.bio(:,st(n):en(n)),2);
    mp_mean(:,n)=nanmean(MP.bio(:,st(n):en(n)),2);
    mf_mean(:,n)=nanmean(MF.bio(:,st(n):en(n)),2);
    md_mean(:,n)=nanmean(MD.bio(:,st(n):en(n)),2);
    lp_mean(:,n)=nanmean(LP.bio(:,st(n):en(n)),2);
    ld_mean(:,n)=nanmean(LD.bio(:,st(n):en(n)),2);
    b_mean(:,n)=nanmean(Bent.bio(:,st(n):en(n)),2);
    
    mp_my(:,n)=nanmean(MP.yield(:,st(n):en(n)),2);
    mf_my(:,n)=nanmean(MF.yield(:,st(n):en(n)),2);
    md_my(:,n)=nanmean(MD.yield(:,st(n):en(n)),2);
    lp_my(:,n)=nanmean(LP.yield(:,st(n):en(n)),2);
    ld_my(:,n)=nanmean(LD.yield(:,st(n):en(n)),2);
end

%
save([fpath 'Means_Hist_ctrlclim_All_fishobs_',mod,'_' cfile '.mat'],'time',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean','b_tmean',...
    'mf_tmcatch','mp_tmcatch','md_tmcatch',...
    'lp_tmcatch','ld_tmcatch',...
    'sf_smean','sp_smean','sd_smean',...
    'mf_smean','mp_smean','md_smean',...
    'lp_smean','ld_smean','b_smean',...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean',...
    'mf_mcatch','mp_mcatch','md_mcatch',...
    'lp_mcatch','ld_mcatch',...
    'mf_my','mp_my','md_my',...
    'lp_my','ld_my')







