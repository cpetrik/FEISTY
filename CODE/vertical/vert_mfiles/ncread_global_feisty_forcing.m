% 0.5 degree global MOM6-COBALTv3 only


clear
close all

fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%%
ncdisp([fpath 'ocean_cobalt_feisty_forcing_z.199001-199412.nmdz.nc'])


%% 
ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_z.199001-199412.nmdz.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
netcdf.close(ncid);










