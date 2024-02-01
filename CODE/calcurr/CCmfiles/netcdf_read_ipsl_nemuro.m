% Read

clear
close all

%%
%fpath='/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/IPSLdown/';
fpath='/Users/cpetrik/Documents/NEMURO/';

ncdisp([fpath 'new_feisty_ipsl_1980-2100.nc'])

%%
% Global Attributes:
% history     = 'FERRET V7.4  29-Oct-23'
% Conventions = 'CF-1.6'
% missing_value = -9.999999999999999e+33
% _FillValue    = -9.999999999999999e+33

% Dimensions:
% XI_RHO   = 186
% ETA_RHO  = 181
% S_RHO1_1 = 1
% TMO      = 1452  (UNLIMITED)
% bnds     = 2

% Variables:
% LON
lon_long_name     = 'LON_RHO';

% LAT
lat_long_name     = 'LAT_RHO';
                       
% BATHY
bathy_long_name     = 'H';
bathy_units         = 'm';

% S_RHO1_1
srho_long_name     = 'S-coordinate at RHO-points';
% positive      = 'up'
% point_spacing = 'even'
% axis          = 'Z'
% standard_name = 'altitude'

% TMO
tmo_units         = 'days since 1900-01-01 00:00:00';
tmo_time_origin   = '01-JAN-1900';
tmo_standard_name = 'time';

% TMO_bnds

% TEMP_BOT
% Size:       186x181x1x1452
tb_long_name = 'TEMP[K=1,GT=TMO@ASN]';
tb_units     = 'C';

% TEMP_AVG_200M
% Size:       186x181x1452
tp_long_name = 'TEMPZ[meanZ=-212.5:-2.5]';
tp_units     = 'C';

% LZOO_INT_200M
lz_long_name     = 'LZOO_200M[sumZ=-1:4.1633E-17]';
lz_units         = 'mmolN/m2';

% PZOO_INT_200M
% Size:       186x181x1452
pz_long_name    = 'PZOO_200M[sumZ=-1:4.1633E-17]';
pz_units        = 'mmolN/m2';

% PON_FLX_100M
% Size:       186x181x1452
pon_long_name   = '(-40)*(PONZ[Z=-90@AVE]-PONZ[Z=-110@AVE])/20';
pon_units       = 'mmolN/m3/d';

%Zooplankton units are mmolN/m2. 
%LZ is mesozooplankton and PZ is euphausiids.
%PON flux units are mmolN/m3/day. 
%The remineralization rate of PON is 0.02 day-1 at 0C with a Q10 relationship.

%%
ncid = netcdf.open([fpath 'new_feisty_ipsl_1980-2100.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
TEMP_BOT = squeeze(TEMP_BOT);

%NaNs on land cells
TEMP_AVG_200M(TEMP_AVG_200M<=-9.9e+33)=nan;
TEMP_BOT(TEMP_BOT<=-9.9e+33)=nan;
LZOO_INT_200M(LZOO_INT_200M<=-9.9e+33)=nan;
PZOO_INT_200M(PZOO_INT_200M<=-9.9e+33)=nan;
PON_FLX_100M(PON_FLX_100M<=-9.9e+33)=nan;

%%
save([fpath 'feisty_ipsl_gridspec.mat'],'BATHY','bathy_units',...
    'bathy_long_name','S_RHO1_1','TMO','TMO_bnds','srho_long_name',...
    'tmo_units','tmo_time_origin','tmo_standard_name','LAT','LON');

save([fpath 'feisty_ipsl_lzoo_1980-2100.mat'],'LZOO_INT_200M',...
    'lz_long_name','lz_units',...
    'S_RHO1_1','TMO','TMO_bnds','srho_long_name','tmo_units',...
    'tmo_time_origin','tmo_standard_name','LAT','LON');

save([fpath 'feisty_ipsl_pon_1980-2100.mat'],'PON_FLX_100M',...
    'pon_long_name','pon_units',...
    'S_RHO1_1','TMO','TMO_bnds','srho_long_name','tmo_units',...
    'tmo_time_origin','tmo_standard_name','LAT','LON');

save([fpath 'feisty_ipsl_pzoo_1980-2100.mat'],'PZOO_INT_200M',...
    'pz_long_name','pz_units',...
    'S_RHO1_1','TMO','TMO_bnds','srho_long_name','tmo_units',...
    'tmo_time_origin','tmo_standard_name','LAT','LON');

save([fpath 'feisty_ipsl_tb_1980-2100.mat'],'TEMP_BOT','tb_long_name',...
    'S_RHO1_1','TMO','TMO_bnds','srho_long_name','tmo_units',...
    'tmo_time_origin','tmo_standard_name','LAT','LON','tb_units');

save([fpath 'feisty_ipsl_tp_1980-2100.mat'],'TEMP_AVG_200M',...
    'tp_long_name','tp_units',...
    'S_RHO1_1','TMO','TMO_bnds','srho_long_name','tmo_units',...
    'tmo_time_origin','tmo_standard_name','LAT','LON');


