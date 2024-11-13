% Read NASA GISS output netcdfs

clear 
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/GISS/VolMIP/';

%%
ncdisp([fpath 'Data2Colleen.nc'])


%%
ncid = netcdf.open([fpath 'APR1925.oijlVolMIPCarbANL1924Control.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Time
yr = (time-time(1)+1)/365;

%%
%save([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat']);