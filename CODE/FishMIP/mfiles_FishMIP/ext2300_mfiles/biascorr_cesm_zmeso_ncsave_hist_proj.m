% Calc Fish-MIP outputs saved as NetCDF
% Historical time period

clear 
close all

%%
fpath='/project/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/hist/';

load([fpath 'cesm2_hist_temp_btm_monthly_1850_2014.mat'],'time');
load([fpath 'cesm2_hist_zooc_150_monthly_1850_2014.mat'],'units_vint');
load([fpath 'cesm2_hist_biascorr_zooc_zmeso_diatfraconly_molC_monthly_1850_2014.mat'],...
    'zmeso_corr','CLAT','CLON');

time_units = 'months since 1601-1-1';

%% Setup netcdf path to store to
fname1 = 'cesm2-waccm_r1i1p1f1_historical_biascorr_zmeso_60arcmin_global_monthly_1850_2014.nc';

file_tpb = [fpath fname1];

[ni,nj,nt] = size(zmeso_corr);

%% Lat & Lon should be vectors
LAT = (CLAT(1,:));
LON = CLON(:,1);

%Use Netcdf4 classic
cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

%% zmeso
ncidSB = netcdf.create(file_tpb,cmode);

time_dim = netcdf.defDim(ncidSB,'time',nt);
lon_dim = netcdf.defDim(ncidSB,'lon',ni);
lat_dim = netcdf.defDim(ncidSB,'lat',nj);

vidtSB = netcdf.defVar(ncidSB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidSB,vidtSB,'long_name','time');
netcdf.putAtt(ncidSB,vidtSB,'standard_name','time');
netcdf.putAtt(ncidSB,vidtSB,'units','months since 1601-01-01' );
netcdf.putAtt(ncidSB,vidtSB,'calendar','gregorian');
netcdf.putAtt(ncidSB,vidtSB,'axis','T');

vidlon = netcdf.defVar(ncidSB,'lon','NC_DOUBLE',lon_dim);
netcdf.putAtt(ncidSB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidSB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidSB,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidSB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidSB,'lat','NC_DOUBLE',lat_dim);
netcdf.putAtt(ncidSB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidSB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidSB,vidlat,'units','degrees_north');
netcdf.putAtt(ncidSB,vidlat,'axis','Y');

vidbioSB = netcdf.defVar(ncidSB,'tpb','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidSB,vidbioSB,'long_name','Mesozooplankton Carbon Concentration 150m Vertical Integral');
netcdf.putAtt(ncidSB,vidbioSB,'units','g m-2' );
netcdf.defVarFill(ncidSB,vidbioSB,false,1.000000020040877e+20);
netcdf.putAtt(ncidSB,vidbioSB,'missing value',1.000000020040877e+20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSB,varid,'_FillValue',1.000000020040877e+20);
netcdf.putAtt(ncidSB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidSB,varid,'institution','UCSD');
netcdf.putAtt(ncidSB,varid,'mesozooplankton substitute','diatom fraction of zooc');
netcdf.putAtt(ncidSB,varid,'bias correction','Heneghan obsGLMM');

netcdf.endDef(ncidSB);

zmeso_corr = single(zmeso_corr);
netcdf.putVar(ncidSB,vidlat,LAT);
netcdf.putVar(ncidSB,vidlon,LON);
netcdf.putVar(ncidSB,vidbioSB,zmeso_corr);
netcdf.putVar(ncidSB,vidtSB,time);

netcdf.close(ncidSB);

%%
ncdisp(file_tpb)

