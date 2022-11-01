% FEISTY output at all locations

clear all
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];

%% loop over cycles
nid = 670589;
ntc = 240;
nta = ntc*6;
yrs = 1841:1960;
st = 1:ntc:nta;
en = ntc:ntc:nta;

Bent.bio = nan*ones(nid,nta);

%%
for c=1:10

    % Benthic material
    ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_empHP_bent_cycle',...
        num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    Bent.bio(:,st(c):en(c)) = biomass;
    clear biomass

end

%% time
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_empHP_time.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

st=1:12:length(time);
en=12:12:length(time);

%%
%Time
b_tmean=mean(Bent.bio,1,'omitnan');
b_smean =mean(Bent.bio,2,'omitnan');
for n=1:length(st)
    b_mean(:,n)=nanmean(Bent.bio(:,st(n):en(n)),2);
end

%%
save([fpath 'Means_PI_ctrlclim_All_fishobs_' cfile '.mat'],...
    'b_tmean','b_smean','b_mean','-append')



