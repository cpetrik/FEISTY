% FEISTY output at all locations

clear 
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];

mod = 'v3.2';

%% time
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_time.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

st=1:12:length(time);
en=12:12:length(time);

%% MF
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
clear biomass

%Time
mf_tmean=mean(MF.bio,1,'omitnan');
%Space
mf_smean=mean(MF.bio,2,'omitnan');
% Each year
for n=1:length(st)
    mf_mean(:,n)=nanmean(MF.bio(:,st(n):en(n)),2);
end

clear MF

%% Benthic material
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

Bent.bio = biomass;
clear biomass

%Time
b_tmean=mean(Bent.bio,1,'omitnan');
b_smean =mean(Bent.bio,2,'omitnan');
for n=1:length(st)
    b_mean(:,n)=nanmean(Bent.bio(:,st(n):en(n)),2);
end

clear Bent


%%
save([fpath 'Means_PI_ctrlclim_All_fishobs_',mod,'_' cfile '.mat'],...
    'mf_tmean','b_tmean',...
    'mf_smean','b_smean',...
    'mf_mean','b_mean','-append')



