% Write netcdf file of GFDL 1/4 degree
% mdz & lgz hploss 100m int
% Regridded by code merging Xiao's with 1 line from Matthias

clear 
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';
qpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';

%% Inputs
load([fpath '19590101-20101231.ocean_cobalt_fluxes_int_FishMIP_CP_remapped.mat']);

load([qpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'])
load([qpath 'lme_gfdl-mom6-cobalt2_15arcmin.mat'],'tlme')

%% netcdf write
% nans to a large number
missing_value = 1.000000020040877e+20;

hploss_nmdz_100(isnan(hploss_nmdz_100)) = missing_value;
hploss_nlgz_100(isnan(hploss_nlgz_100)) = missing_value;

%% molN/m2/s --> gC/m2/s
hploss_nmdz_100 = (hploss_nmdz_100 * (106/16) * 12.01);
hploss_nlgz_100 = (hploss_nlgz_100 * (106/16) * 12.01);

hploss_lz100_units        = 'gC/m2/s';
hploss_mz100_units        = 'gC/m2/s';
hploss_lz100_long_name    = 'large zooplankton hploss rate integrated over top 100m';
hploss_mz100_long_name    = 'medium zooplankton hploss rate integrated over top 100m';

%% Quick look
mz = hploss_nmdz_100(:,:,50) *60*60*24*365; 
lz = hploss_nlgz_100(:,:,50) *60*60*24*365;

figure(1)
pcolor(log10(mz))
shading flat
cmocean('matter')
colorbar
caxis([-2 2])
title('MZ')

figure(2)
pcolor(log10(lz))
shading flat
cmocean('matter')
colorbar
caxis([-2 2])
title('LZ')

figure(3)
pcolor((deptho))
shading flat
colorbar
title('depth')

figure(4)
pcolor((tlme))
shading flat
colorbar
title('LME')

%%
close all

%% Setup netcdf path to store to
file_tfb = [fpath '19590101-20101231.ocean_cobalt_hploss_int100.nc'];

[ni,nj,nt] = size(hploss_nmdz_100);

%%
% LAT = single(LAT);
% LON = single(LON);

%Use Netcdf4 classic
% cmode = netcdf.getConstant('NETCDF4');
% cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));
cmode = netcdf.getConstant('64BIT_OFFSET');

%% 
ncidFB = netcdf.create(file_tfb,cmode);

time_dim = netcdf.defDim(ncidFB,'time',nt);
lon_dim = netcdf.defDim(ncidFB,'nlon',ni);
lat_dim = netcdf.defDim(ncidFB,'nlat',nj);

vidtFB = netcdf.defVar(ncidFB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidFB,vidtFB,'long_name','time');
netcdf.putAtt(ncidFB,vidtFB,'standard_name','time');
netcdf.putAtt(ncidFB,vidtFB,'calendar','365_day');
netcdf.putAtt(ncidFB,vidtFB,'axis','T');
netcdf.putAtt(ncidFB,vidtFB,'units',time_units);

vidlon = netcdf.defVar(ncidFB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidFB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidFB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidFB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidFB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidFB,vidlat,'axis','Y');

vidMZ = netcdf.defVar(ncidFB,'hploss_nmdz_100','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
%netcdf.defVarChunking(ncidFB,vidMZ,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidFB,vidMZ,'long_name',hploss_mz100_long_name);
netcdf.putAtt(ncidFB,vidMZ,'units',hploss_mz100_units);
netcdf.defVarFill(ncidFB,vidMZ,false,missing_value);
netcdf.putAtt(ncidFB,vidMZ,'missing value',missing_value);

vidMZL = netcdf.defVar(ncidFB,'hploss_nlgz_100','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
%netcdf.defVarChunking(ncidFB,vidMZL,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidFB,vidMZL,'long_name',hploss_lz100_long_name);
netcdf.putAtt(ncidFB,vidMZL,'units',hploss_lz100_units);
netcdf.defVarFill(ncidFB,vidMZL,false,missing_value);
netcdf.putAtt(ncidFB,vidMZL,'missing value',missing_value);

vidFrac = netcdf.defVar(ncidFB,'lme_mask','NC_FLOAT',[lon_dim,lat_dim]);
%netcdf.defVarChunking(ncidFB,vidFrac,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidFB,vidFrac,'long_name','LME mask');
netcdf.putAtt(ncidFB,vidFrac,'units','');
netcdf.defVarFill(ncidFB,vidFrac,false,missing_value);
netcdf.putAtt(ncidFB,vidFrac,'missing value',missing_value);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidFB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidFB,varid,'_FillValue',missing_value);
netcdf.putAtt(ncidFB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidFB,varid,'institution','UCSD');

netcdf.endDef(ncidFB);

netcdf.putVar(ncidFB,vidlat,LAT);
netcdf.putVar(ncidFB,vidlon,LON);
netcdf.putVar(ncidFB,vidtFB,time);
netcdf.putVar(ncidFB,vidMZ,hploss_nmdz_100);
netcdf.putVar(ncidFB,vidMZL,hploss_nlgz_100);
netcdf.putVar(ncidFB,vidFrac,tlme);

netcdf.close(ncidFB);

%%
ncdisp(file_tfb)

