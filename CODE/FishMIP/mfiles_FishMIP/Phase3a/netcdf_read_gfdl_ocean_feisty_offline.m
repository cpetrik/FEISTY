% Read GFDL 1/4 netcdfs
% Thetao 3D
% Raw units provided by GFDL
% Before conversion by Matthias

clear 
close all

%fpath='/Volumes/MIP/Fish-MIP/Phase3/GFDL_reanalysis/';
fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/BATS_1D/';

%% one file
ncdisp([fpath '20040101.ocean_feisty_forcing_offline.nc'])

%%

%%
ncid = netcdf.open([fpath '20040101.ocean_feisty_forcing_offline.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

netcdf.close(ncid);

%% 
save([fpath '20040101.ocean_feisty_forcing_offline.mat']);

%% viz means
[ni, nj, nk, nt] = size(nmdz);

medz = reshape(nmdz,ni*nj,nk,nt);
test = squeeze(double(medz(1,:,:)));

