% Read Fish-MIP CESM Historic netcdfs
% Integrate over top 100 m for small zooplankton

clear all
close all

fpath='/Volumes/GFDL/Fish-MIP/CESM/Hist/';
cpath='/Volumes/GFDL/Fish-MIP/GFDL/Hist/';


%% CESM
%mmolC/m2
ncid = netcdf.open([fpath 'cesm_hist_szoo_zint_monthly_',num2str(tstart(t)),'-',...
    num2str(tend(t)),'.nc4'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%NaNs on land cells
szoo(szoo>=9e+36)=nan;
%top 100 m
sz100 = szoo(:,:,1:10,:);
%integrated over top 100m
isz100 = squeeze(nansum(sz100,3));



%% COBALT
%molC/m2
ncid = netcdf.open([fpath 'cesm_hist_szoo_zint_monthly_',num2str(tstart(t)),'-',...
    num2str(tend(t)),'.nc4'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%NaNs on land cells
szoo(szoo>=9e+36)=nan;
%top 100 m
sz100 = szoo(:,:,1:10,:);
%integrated over top 100m
isz100 = squeeze(nansum(sz100,3));




