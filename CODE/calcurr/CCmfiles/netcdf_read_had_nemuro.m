% Read

clear
close all

%%
fpath='/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/HADdown/';
%fpath='/Users/cpetrik/Documents/NEMURO/';

ncdisp([fpath 'feisty_hadley_1980-2100.nc'])

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

% TMO
tmo_units         = 'days since 1900-01-01 00:00:00';
tmo_time_origin   = '01-JAN-1900';
tmo_standard_name = 'time';

% TMO_bnds

% TEMP_BOT
% Size:       186x181x1452
tb_long_name = 'TEMP[K=1,GT=TMO@ASN]';
tb_units     = 'C';

% TEMP_AVG_200M
% Size:       186x181x1452
tp_long_name = 'TEMPZ[meanZ=-212.5:-2.5]';
tp_units     = 'C';

% LZOO_INT_200M
% Size:       186x181x1452
lz_long_name     = 'LZOO_200M[sumZ=-1:4.1633E-17]';
lz_units         = 'mmolN/m2';

% PZOO_INT_200M
% Size:       186x181x1452
pz_long_name    = 'PZOO_200M[sumZ=-1:4.1633E-17]';
pz_units        = 'mmolN/m2';

% PON_BOT
% Size:       186x181x1452
% missing_value = -9.999999999999999e+33;
% _FillValue    = -9.999999999999999e+33;
pon_long_name     = 'PON[K=1@AVE,GT=TMO@ASN]';
pon_units       = 'mmolN/m3';
% history       = 'From /data01/fiechter/WC12_OffBio/HAD_TV_1980-2100/Monthly/wc12_avg_1980.nc'


%Zooplankton units are mmolN/m2. 
%LZ is mesozooplankton and PZ is euphausiids.
%PON flux units are mmolN/m3/day. 

%%
ncid = netcdf.open([fpath 'feisty_hadley_1980-2100.nc'],'NC_NOWRITE');

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
PON_BOT(PON_BOT<=-9.9e+33)=nan;

%%
save([fpath 'feisty_hadley_gridspec.mat'],'BATHY','bathy_units',...
    'bathy_long_name','TMO','TMO_bnds',...
    'tmo_units','tmo_time_origin','tmo_standard_name','LAT','LON');

save([fpath 'feisty_hadley_lzoo_1980-2100.mat'],'LZOO_INT_200M',...
    'lz_long_name','lz_units',...
    'TMO','TMO_bnds','tmo_units',...
    'tmo_time_origin','tmo_standard_name','LAT','LON');

save([fpath 'feisty_hadley_pon_1980-2100.mat'],'PON_BOT',...
    'pon_long_name','pon_units',...
    'TMO','TMO_bnds','tmo_units',...
    'tmo_time_origin','tmo_standard_name','LAT','LON');

save([fpath 'feisty_hadley_pzoo_1980-2100.mat'],'PZOO_INT_200M',...
    'pz_long_name','pz_units',...
    'TMO','TMO_bnds','tmo_units',...
    'tmo_time_origin','tmo_standard_name','LAT','LON');

save([fpath 'feisty_hadley_tb_1980-2100.mat'],'TEMP_BOT','tb_long_name',...
    'TMO','TMO_bnds','tmo_units',...
    'tmo_time_origin','tmo_standard_name','LAT','LON','tb_units');

save([fpath 'feisty_hadley_tp_1980-2100.mat'],'TEMP_AVG_200M',...
    'tp_long_name','tp_units',...
    'TMO','TMO_bnds','tmo_units',...
    'tmo_time_origin','tmo_standard_name','LAT','LON');


