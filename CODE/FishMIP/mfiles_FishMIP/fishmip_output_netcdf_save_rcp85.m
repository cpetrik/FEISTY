% Calc Fish-MIP outputs saved as NetCDF
% Historic simulation

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'pristine';
fpath=['/Volumes/GFDL/NC/FishMIP/CESM1-BEC/' cfile '/'];

%%  Outputs
%Total system carbon biomass, tsb, gCm-2, All primary producers and consumers
%Total consumer carbon biomass density, tcb, gCm-2 All consumers (trophic level >1, vertebrates and invertebrates)
%Carbon biomass density of consumers > 10cm, b10cm, gCm-2 If asymptotic length (Linf) is > 10cm, include in > 10cm class
%Carbon biomass density of consumers > 30cm, b30cm, gCm-2 If asymptotic length (Linf) is > 30cm, include in > 30cm class
%Monthly or annual

%tsb = phy + zoo + F + P + D + B
%tcb = zoo + F + P + D + B
%b10cm = M + L + (0.1)*B
%b30cm = L

load([fpath 'FishMIP_output_Forecast_pristine_' cfile '.mat']);
[ni,nj,nt] = size(tsb);

%NaN fills
tsb(isnan(tsb))=1e36;
tcb(isnan(tcb))=1e36;
b10cm(isnan(b10cm))=1e36;
b30cm(isnan(b30cm))=1e36;

%% Output data naming conventions
% For Runs 1,2 and 3 (preindustrial, historical and rcp85)
% ? Model: e.g., ?apecosm?
% ? Forcing GCM/reanalysis model: ?cesm1-bgc? or ?gfdl-esm2m?
% ? Bias correction (none): ?nobc?
% ? Climate Scenario: ?pre-industrial? or ?historical? or ?rcp85?
% ? Human impacts (fishing ? none): ?nosoc?
% ? Default CO2 scenario: ?co2?
% ? Variable name (see table): e.g., ?b10cm?
% ? Region: ?global?, or if applicable the region e.g., ?baltic-sea?
% ? Temporal resolution: ?monthly? or ?annual?
% ? First year of reporting period - Last year of reporting period: e.g., 1860-1869
% ? Netcdf-4: ?.nc4?
% e.g. apecosm_cesm1-bgc_nobc_historical_nosoc_co2_b10cm_global_monthly_1860-1869.nc4

% For Runs 4 and 5 (NPP control, temperature control)
% ? Model: e.g., ?apecosm?
% ? Forcing GCM/reanalysis model: ?cesm1-bgc? or ?gfdl-esm2m?
% ? Bias correction (none): ?nobc?
% ? Scenario: ?npp-control? or ?temperature-control?
% ? Human impacts (fishing ? none): ?nosoc?
% ? Default CO2 scenario: ?co2?
% ? Variable name (see table): e.g., ?b10cm?
% ? Region: ?global?, or if applicable the region e.g., ?baltic-sea?
% ? Temporal resolution: ?monthly? or ?annual?
% ? First year of reporting period - Last year of reporting period: e.g., 1860-1869
% ? Netcdf-4: ?.nc4?
% e.g apecosm_cesm1-bgc_nobc_npp-control_nosoc_co2_b10cm_global_monthly_1860- 1869.nc4

% Output netcdf metadata and other information
% ? In the metadata for the output data, include lines called ?pH input used (acidification)? and ?Diazotroph input used:?.
% -> For pH input used (acidification): ?yes? or ?no?
% -> For diazotroph input used: ?no?, or ?yes, diazotroph carbon?, ?yes, diaztroph production?, or ?yes, diazotroph carbon and production?.
% If you use integrated primary production and integrated phytoplankton biomass, say ?only if included with integrated phytoplankton production (and/or) biomass?
% ? Don?t omit your default values of specifiers; always give a value

%% Setup netcdf path to store to
fname1 = 'feisty_cesm1-bgc_nobc_rcp85_nosoc_co2_';
fname2 = '_global_annual_2006-2100.nc4';

file_tsb = [fpath fname1 'tsb' fname2];
file_tcb = [fpath fname1 'tcb' fname2];
file_b10 = [fpath fname1 'b10cm' fname2];
file_b30 = [fpath fname1 'b30cm' fname2];

%% tsb
ncidSB = netcdf.create(file_tsb,'NETCDF4');

lon_dim = netcdf.defDim(ncidSB,'longitude',ni);
lat_dim = netcdf.defDim(ncidSB,'latitude',nj);
time_dim = netcdf.defDim(ncidSB,'time',nt);

vidlat = netcdf.defVar(ncidSB,'LAT','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidSB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidSB,vidlat,'units','degrees');

vidlon = netcdf.defVar(ncidSB,'LON','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidSB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidSB,vidlon,'units','degrees' );

vidtSB = netcdf.defVar(ncidSB,'time','double',time_dim);
netcdf.putAtt(ncidSB,vidtSB,'long_name','time');
netcdf.putAtt(ncidSB,vidtSB,'units','years' );

vidbioSB = netcdf.defVar(ncidSB,'tsb','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidSB,vidbioSB,'long_name','total system carbon biomass density');
netcdf.putAtt(ncidSB,vidbioSB,'units','grams Carbon m-2' );
netcdf.defVarFill(ncidSB,vidbioSB,false,1.0e36);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSB,varid,'pH_input_used_(acidification):','no');
netcdf.putAtt(ncidSB,varid,'Diazotroph_input_used:',...
    'only if included with integrated phytoplankton biomass');

netcdf.endDef(ncidSB);

netcdf.putVar(ncidSB,vidlat,LAT);
netcdf.putVar(ncidSB,vidlon,LON);
netcdf.putVar(ncidSB,vidbioSB,tsb);
netcdf.putVar(ncidSB,vidtSB,Years);

netcdf.close(ncidSB);

%% tcb
ncidCB = netcdf.create(file_tcb,'NETCDF4');

lon_dim = netcdf.defDim(ncidCB,'longitude',ni);
lat_dim = netcdf.defDim(ncidCB,'latitude',nj);
time_dim = netcdf.defDim(ncidCB,'time',nt);

vidlat = netcdf.defVar(ncidCB,'LAT','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidCB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidCB,vidlat,'units','degrees');

vidlon = netcdf.defVar(ncidCB,'LON','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidCB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidCB,vidlon,'units','degrees' );

vidtCB = netcdf.defVar(ncidCB,'time','double',time_dim);
netcdf.putAtt(ncidCB,vidtCB,'long_name','time');
netcdf.putAtt(ncidCB,vidtCB,'units','years' );

vidbioCB = netcdf.defVar(ncidCB,'tcb','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidCB,vidbioCB,'long_name','total consumer carbon biomass density');
netcdf.putAtt(ncidCB,vidbioCB,'units','grams Carbon m-2' );
netcdf.defVarFill(ncidCB,vidbioCB,false,1.0e36);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidCB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidCB,varid,'pH_input_used_(acidification):','no');
netcdf.putAtt(ncidCB,varid,'Diazotroph_input_used:',...
    'only if included with integrated phytoplankton biomass');

netcdf.endDef(ncidCB);

netcdf.putVar(ncidCB,vidlat,LAT);
netcdf.putVar(ncidCB,vidlon,LON);
netcdf.putVar(ncidCB,vidbioCB,tcb);
netcdf.putVar(ncidCB,vidtCB,Years);

netcdf.close(ncidCB);

%% b10cm
ncid10 = netcdf.create(file_b10,'NETCDF4');

lon_dim = netcdf.defDim(ncid10,'longitude',ni);
lat_dim = netcdf.defDim(ncid10,'latitude',nj);
time_dim = netcdf.defDim(ncid10,'time',nt);

vidlat = netcdf.defVar(ncid10,'LAT','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncid10,vidlat,'long_name','latitude');
netcdf.putAtt(ncid10,vidlat,'units','degrees');

vidlon = netcdf.defVar(ncid10,'LON','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncid10,vidlon,'long_name','longitude');
netcdf.putAtt(ncid10,vidlon,'units','degrees' );

vidt10 = netcdf.defVar(ncid10,'time','double',time_dim);
netcdf.putAtt(ncid10,vidt10,'long_name','time');
netcdf.putAtt(ncid10,vidt10,'units','years' );

vidbio10 = netcdf.defVar(ncid10,'b10cm','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncid10,vidbio10,'long_name','carbon biomass density of consumers > 10cm');
netcdf.putAtt(ncid10,vidbio10,'units','grams Carbon m-2' );
netcdf.defVarFill(ncid10,vidbio10,false,1.0e36);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid10,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid10,varid,'pH_input_used_(acidification):','no');
netcdf.putAtt(ncid10,varid,'Diazotroph_input_used:',...
    'only if included with integrated phytoplankton biomass');

netcdf.endDef(ncid10);

netcdf.putVar(ncid10,vidlat,LAT);
netcdf.putVar(ncid10,vidlon,LON);
netcdf.putVar(ncid10,vidbio10,b10cm);
netcdf.putVar(ncid10,vidt10,Years);

netcdf.close(ncid10);

%% b30cm
ncid30 = netcdf.create(file_b30,'NETCDF4');

lon_dim = netcdf.defDim(ncid30,'longitude',ni);
lat_dim = netcdf.defDim(ncid30,'latitude',nj);
time_dim = netcdf.defDim(ncid30,'time',nt);

vidlat = netcdf.defVar(ncid30,'LAT','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncid30,vidlat,'long_name','latitude');
netcdf.putAtt(ncid30,vidlat,'units','degrees');

vidlon = netcdf.defVar(ncid30,'LON','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncid30,vidlon,'long_name','longitude');
netcdf.putAtt(ncid30,vidlon,'units','degrees' );

vidt30 = netcdf.defVar(ncid30,'time','double',time_dim);
netcdf.putAtt(ncid30,vidt30,'long_name','time');
netcdf.putAtt(ncid30,vidt30,'units','years' );

vidbio30 = netcdf.defVar(ncid30,'b30cm','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncid30,vidbio30,'long_name','carbon biomass density of consumers > 30cm');
netcdf.putAtt(ncid30,vidbio30,'units','grams Carbon m-2' );
netcdf.defVarFill(ncid30,vidbio30,false,1.0e36);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid30,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid30,varid,'pH_input_used_(acidification):','no');
netcdf.putAtt(ncid30,varid,'Diazotroph_input_used:',...
    'only if included with integrated phytoplankton biomass');

netcdf.endDef(ncid30);

netcdf.putVar(ncid30,vidlat,LAT);
netcdf.putVar(ncid30,vidlon,LON);
netcdf.putVar(ncid30,vidbio30,b30cm);
netcdf.putVar(ncid30,vidt30,Years);

netcdf.close(ncid30);











