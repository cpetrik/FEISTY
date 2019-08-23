% Calc Fish-MIP outputs and save as NetCDF

clear all
close all

%% Fish data
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP_CESM/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'pristine';
tharv = 'F=0';
fpath=['/Volumes/GFDL/NC/FishMIP/CESM1-BEC/' cfile '/'];

load([fpath 'Means_Historic_' harv '_' cfile '.mat']);

hpath=['/Volumes/GFDL/Fish-MIP/CESM/hist/'];
load([hpath 'cesm_hist_phy_zint_monthly_185001-200512.mat'])
load([hpath 'cesm_hist_zoo_zint_monthly_185001-200512.mat'])

Years = 1850:2005;

%% Outputs
%Total system carbon biomass, tsb, gCm-2, All primary producers and consumers
%Total consumer carbon biomass density, tcb, gCm-2 All consumers (trophic level >1, vertebrates and invertebrates)
%Carbon biomass density of consumers > 10cm, b10cm, gCm-2 If asymptotic length (Linf) is > 10cm, include in > 10cm class
%Carbon biomass density of consumers > 30cm, b30cm, gCm-2 If asymptotic length (Linf) is > 30cm, include in > 30cm class
%Monthly or annual

%tsb = phy + zoo + F + P + D + B
%tcb = zoo + F + P + D + B
%b10cm = M + L + (0.1-0.25)*B
%b30cm = L

All  = sp_tot + sf_tot + sd_tot + mp_tot + mf_tot + md_tot + lp_tot + ld_tot;
AllF = sf_tot + mf_tot;
AllP = sp_tot + mp_tot + lp_tot;
AllD = sd_tot + md_tot + ld_tot;
AllS = sp_tot + sf_tot + sd_tot;
AllM = mp_tot + mf_tot + md_tot;
AllL = lp_tot + ld_tot;

vb10cm = AllM + AllL + (0.1)*b_tot;
vb30cm = AllL;

%% Units
%fish: gWW/m^2
%zoo_int: mmol C/m^2
%phy_int: mmol C/m^2
%tp: degC
%tb: degC

% Convert to gC m-2 from gWW m-2
vb10cm = (1/9) * vb10cm;
vb10cm = (1/9) * vb30cm;

% Convert to gC m-2 from mmol C m-2
% 1e-3 mol in 1 mmol
% 12.01 g C in 1 mol C
zint = nz_int * 1e-3 * 12.01;
pint = np_int * 1e-3 * 12.01;

%% Reshape to ocean cells only
cpath = '/Volumes/GFDL/Fish-MIP/CESM/';
load([cpath 'gridspec_cesm.mat'],'LAT','LON');
load([cpath 'Data_grid_cesm.mat']);

[ni,nj,NT] = size(nz_int);

zoo = reshape(zint,ni*nj,NT);
phy = reshape(pint,ni*nj,NT);

vzoo = zoo(GRD.ID,:);
vphy = phy(GRD.ID,:);

%% Annual totals of phy and zoo
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
t=time;
mo=t/12;
mo=mo+1850;
st=1:12:length(time);
en=12:12:length(time);

phy_tot = nan(size(sd_tot));
zoo_tot = nan(size(sd_tot));
for m=1:length(en)
    yr1 = st(m):en(m);
    phy_tot(:,m)=sum(vphy(:,yr1).*MNTH,2);
    zoo_tot(:,m)=sum(vzoo(:,yr1).*MNTH,2);
end

%%
vtsb = phy_tot + zoo_tot + AllF + AllP + AllD + b_tot;
vtcb = zoo_tot + AllF + AllP + AllD + b_tot;

%% Reshape to lat,lon,yr
[nid,nt] = size(sd_tot);

tsb = 1.0e36*ones(ni,nj,nt);
tcb = 1.0e36*ones(ni,nj,nt);
b10cm = 1.0e36*ones(ni,nj,nt);
b30cm = 1.0e36*ones(ni,nj,nt);

for y=1:length(en)
    gtsb = nan(ni,nj);
    ttsb = vtsb(:,y);
    gtsb(GRD.ID) = ttsb;
    tsb(:,:,y) = gtsb;
    
    gtcb = nan(ni,nj);
    ttcb = vtcb(:,y);
    gtcb(GRD.ID) = ttcb;
    tcb(:,:,y) = gtcb;
    
    gb10cm = nan(ni,nj);
    tb10cm = vb10cm(:,y);
    gb10cm(GRD.ID) = tb10cm;
    b10cm(:,:,y) = gb10cm;
    
    gb30cm = nan(ni,nj);
    tb30cm = vb30cm(:,y);
    gb30cm(GRD.ID) = tb30cm;
    b30cm(:,:,y) = gb30cm;
end

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
fname1 = 'feisty_cesm1-bgc_nobc_historical_nosoc_co2_';
fname2 = '_global_annual_1850-2005.nc';

file_tsb = [fpath fname1 'tsb' fname2];
file_tcb = [fpath fname1 'tcb' fname2];
file_b10 = [fpath fname1 'b10cm' fname2];
file_b30 = [fpath fname1 'b30cm' fname2];

%% A snippet of code to write a netCDF file.  You need to read
% the time variable and unit information from another file.

% xgrid, ygrid, and level are vectors that are defined elsewhere
% fnin, fnout are filenames for existing and new files, respectively

% var and varname are the variable to be written and its name, respectively

nccreate( fnout, varname, 'Dimensions', { 'lon', length(xgrid), ...
  'lat', length(ygrid), 'level', length(level), 'time', inf } );

ncwrite(    'fnout.nc', 'w', w );
ncwriteatt( 'fnout.nc', 'w', 'long_name', 'vertical velocity' );
ncwriteatt( 'fnout.nc', 'w', 'units', 'Pa/s' );

nccreate(   'fnout.nc', 'lon', 'Dimensions', { 'x', 211, 'y', 121 }, ...
  'DataType', 'single' );
ncwrite(    'fnout.nc', 'lon', xgrid2 );    % "xgrid2" is the 2-dimensional matrix of longitudes
ncwriteatt( 'fnout.nc', 'lon', 'long_name', 'longitude' );
ncwriteatt( 'fnout.nc', 'lon', 'units', 'degrees_east' );

nccreate(   'fnout.nc', 'lat', 'Dimensions', { 'x', 211, 'y', 121 }, ...
  'DataType', 'single' );
ncwrite(    'fnout.nc', 'lat', ygrid2 );    % "ygrid2" is the 2-dimensional matrix of latitudes 
ncwriteatt( 'fnout.nc', 'lat', 'long_name', 'latitude' );
ncwriteatt( 'fnout.nc', 'lat', 'units', 'degrees_north' );

nccreate(   'fnout.nc', 'time', 'Dimensions', { 'time', inf } );
ncwrite(    'fnout.nc', 'time', time_value );   % a vector of time values relative to a reference time 
ncwriteatt( 'fnout.nc', 'time', 'long_name', 'time' );
ncwriteatt( 'fnout.nc', 'time', 'units', ...
  'seconds since 1970-01-01 00:00:00.0 0:00' );

%% tsb
nccreate( file_tsb, tsb, 'Dimensions', { 'lon', ni, ...
  'lat', nj, 'time', nt } );
ncwrite(    file_tsb, 'tsb', tsb );
ncwriteatt( file_tsb, 'tsb', 'long_name', 'total system biomass' );
ncwriteatt( file_tsb, 'tsb', 'units', 'gramsCarbon/m2' );

nccreate(   file_tsb, 'lon', 'Dimensions', { 'x', ni, 'y', nj }, ...
  'DataType', 'single' );
ncwrite(    file_tsb, 'lon', LON );     
ncwriteatt( file_tsb, 'lon', 'long_name', 'longitude' );
ncwriteatt( file_tsb, 'lon', 'units', 'degrees_east' );

nccreate(   file_tsb, 'lat', 'Dimensions', { 'x', ni, 'y', nj }, ...
  'DataType', 'single' );
ncwrite(    file_tsb, 'lat', LAT );    
ncwriteatt( file_tsb, 'lat', 'long_name', 'latitude' );
ncwriteatt( file_tsb, 'lat', 'units', 'degrees_north' );

nccreate(   file_tsb, 'time', 'Dimensions', { 'time', nt } );
ncwrite(    file_tsb, 'time', Years );   % a vector of time values relative to a reference time 
ncwriteatt( file_tsb, 'time', 'long_name', 'time' );
ncwriteatt( file_tsb, 'time', 'units', 'year' );

%If the variable has missing values, you define the FillValue in 
%nccreate: " 'FillValue', value ". Do not use NAN as the FillValue. 
%The Matlab documentation uses NaN, but netCDF operators (NCOs) may or may 
%not interpret the NAN correctly. Instead use a number for the missing value 
%that won't be found in your data.


%% Def vars of netcdf file
%netcdf.setDefaultFormat('NC_FORMAT_64BIT');
ncidSB = netcdf.create(file_tsb,'NC_WRITE');
ncidCB = netcdf.create(file_tcb,'NC_WRITE');
ncid10 = netcdf.create(file_b10,'NC_WRITE');
ncid30 = netcdf.create(file_b30,'NC_WRITE');

lon_dim     = netcdf.defDim(ncidSB,'longitude',ni);
lat_dim     = netcdf.defDim(ncidSB,'latitude',nj);
time_dim    = netcdf.defDim(ncidSB,'time',nt);
%fill_dim    = netcdf.defDim(ncidSB,'fill',1);
fill_dim    = netcdf.defDim(ncidSB,'fill value',1.0e36);
%vidfill     = netcdf.defVar(ncidSB,'fill value',1.0e36);
vidbioSB    = netcdf.defVar(ncidSB,'biomass gCm-2','double',[lon_dim,lat_dim,time_dim]);
vidtSB      = netcdf.defVar(ncidSB,'year','double',time_dim);
netcdf.endDef(ncidSB);

lon_dim     = netcdf.defDim(ncidCB,'longitude',ni);
lat_dim     = netcdf.defDim(ncidCB,'latitude',nj);
time_dim    = netcdf.defDim(ncidCB,'time',nt);
fill_dim    = netcdf.defDim(ncidCB,'fill value',1.0e36);
vidbioCB    = netcdf.defVar(ncidCB,'biomass gCm-2','double',[lon_dim,lat_dim,time_dim]);
vidtCB      = netcdf.defVar(ncidCB,'year','double',time_dim);
netcdf.endDef(ncidCB);

lon_dim     = netcdf.defDim(ncid10,'longitude',ni);
lat_dim     = netcdf.defDim(ncid10,'latitude',nj);
time_dim    = netcdf.defDim(ncid10,'time',nt);
fill_dim    = netcdf.defDim(ncid10,'fill value',1.0e36);
vidbio10    = netcdf.defVar(ncid10,'biomass gCm-2','double',[lon_dim,lat_dim,time_dim]);
vidt10      = netcdf.defVar(ncid10,'year','double',time_dim);
netcdf.endDef(ncid10);

lon_dim     = netcdf.defDim(ncid30,'longitude',ni);
lat_dim     = netcdf.defDim(ncid30,'latitude',nj);
time_dim    = netcdf.defDim(ncid30,'time',nt);
fill_dim    = netcdf.defDim(ncid30,'fill',1);
vidfill     = netcdf.defVar(ncid30,'fill value',1.0e36);
vidbio30    = netcdf.defVar(ncid30,'biomass gCm-2','double',[lon_dim,lat_dim,time_dim]);
vidt30      = netcdf.defVar(ncid30,'year','double',time_dim);
netcdf.endDef(ncid30);

%% Put
years = 1850:2005;

netcdf.putVar(ncidSB,vidbioSB,tsb);
netcdf.putVar(ncid10,vidt10,years);

netcdf.putVar(ncidCB,vidbioCB,tcb);
netcdf.putVar(ncid10,vidt10,years);

netcdf.putVar(ncid10,vidbio10,b10cm);
netcdf.putVar(ncid10,vidt10,years);

netcdf.putVar(ncid30,vidbio30,b30cm);
netcdf.putVar(ncid10,vidt10,years);

%% Close save
netcdf.close(ncidSB);
netcdf.close(ncidCB);
netcdf.close(ncid10);
netcdf.close(ncid30);

