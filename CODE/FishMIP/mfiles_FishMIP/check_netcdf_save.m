clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'pristine';
fpath=['/Volumes/GFDL/NC/FishMIP/CESM1-BEC/' cfile '/'];

fname1 = 'feisty_cesm1-bgc_nobc_historical_nosoc_co2_';
fname2 = '_global_annual_1850-2005.nc4';

file_tsb = [fpath fname1 'tsb' fname2];
file_tcb = [fpath fname1 'tcb' fname2];
file_b10 = [fpath fname1 'b10cm' fname2];
file_b30 = [fpath fname1 'b30cm' fname2];

%%
ncid = netcdf.open(file_tsb,'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);

test=tsb(:,:,10);

%%
ncid = netcdf.open(file_tcb,'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);

test2=tcb(:,:,50);

%%
ncid = netcdf.open(file_b10,'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);

test3=b10cm(:,:,100);

%%
ncid = netcdf.open(file_b30,'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);

test4=b30cm(:,:,150);