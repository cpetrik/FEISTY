% Netcdf of area for MOM6 COBALT2 sims
% prepared by Matthias at ISIMIP

clear all
close all

%%
Pdir = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';

ncid = netcdf.open([Pdir 'cellarea_onedeg.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%%
[LAT,LON] = meshgrid(lat,lon);

figure
pcolor(LON,LAT,cell_area); shading flat
colorbar

figure
pcolor(LON); shading flat
colorbar

figure
pcolor(LAT); shading flat
colorbar

%%
save([Pdir 'cellarea_onedeg.mat'],'cell_area','LAT','LON');


