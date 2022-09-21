% Netcdf of LME mask for MOM6 COBALT2 sims
% prepared by Matthias at ISIMIP

clear all
close all

%%
Pdir = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';

ncid = netcdf.open([Pdir 'LMEs66_masks_0.25deg.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([Pdir 'cellarea_15arcmin.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);


%% 
% fuck this, I have to loop through all 66 masks to put on one grid???


%%
save([Pdir 'lme_gfdl-mom6-cobalt2_onedeg.mat'],'tlme');%,'AREA_OCN');


