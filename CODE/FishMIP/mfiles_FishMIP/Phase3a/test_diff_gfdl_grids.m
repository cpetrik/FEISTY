% Zoop files from Xiao seem to be on diff grid than ISIMIP 1/4 files
% Xiao: 1440x1080 
% ISIMIP: 1440x720
% orig seems like there is "land" in the arctic for zoop, but is also there
% for phys vars so not concerned
% but how to regrid to match other isimip vars?

gpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/GFDL_reanalysis/';

load([gpath 'gridspec_gfdl-mom6-cobalt2_obsclim_15arcmin_orig.mat'])

ncdisp([gpath '20000101.ocean_scalar_month_FishMIP.nc']) %SW rad

ncdisp([gpath '20000101.ocean_month_z_FishMIP.nc']) %sal, vel

ncdisp([gpath '20000101.ocean_month_z_FishMIP_b.nc']) %thetao

%%
scid = netcdf.open([gpath '20000101.ocean_month_z_FishMIP.nc'],'NC_NOWRITE');
[sdims,svars,sgatts,unlimdimid] = netcdf.inq(scid);

for s = 1:svars
    varname = netcdf.inqVar(scid, s-1);
    eval([ varname ' = netcdf.getVar(scid,s-1);']);
 
end
netcdf.close(scid);

figure
sest = double(squeeze(so(:,:,1,6)));
pcolor(sest); shading flat

figure
uest = double(squeeze(uo(:,:,1,6)));
pcolor(uest); shading flat

%%
tcid = netcdf.open([gpath '20000101.ocean_month_z_FishMIP_b.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

for t = 1:tvars
    varname = netcdf.inqVar(tcid, t-1);
    eval([ varname ' = netcdf.getVar(tcid,t-1);']);
   
end
netcdf.close(tcid);

figure
test = double(squeeze(thetao(:,:,1,6)));
pcolor(test); shading flat

%%
ncid = netcdf.open([gpath '20000101.ocean_scalar_month_FishMIP.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    
end
netcdf.close(ncid);
