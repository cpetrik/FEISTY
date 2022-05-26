% Read CMIP6 netcdfs
% GFDL-ESM4 Hist
% zooc-vint
% 1/4 grid from ISIMIP

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/Phase3/QuarterDeg/';
ppath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/Phase3/quarter_degree_vetting/';

%% zooc zint
ncdisp([fpath 'gfdl-mom6-cobalt2_obsclim_zooc-vint_15arcmin_global_monthly_1961_2010.nc'])

%%
% Dimensions:
% time = 600   (UNLIMITED)
% lon  = 1440
% lat  = 720
% Variables:

% time
% units         = 'months since 1901-1-1 00:00:00'
% calendar      = '360_day'

% zmeso-vint
% Size:       1440x720x600
% Dimensions: lon,lat,time
% Datatype:   single
% Attributes:
FillValue     = 1.000000020040877e+20;
missing_value = 1.000000020040877e+20;
standard_name = 'mole_content_of_zooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Zooplankton Carbon Content';
units         = 'mol m-2';
% comment       = 'vertically integrated over all ocean levels by ISIMIP data management team'

%%
ncid = netcdf.open([fpath 'gfdl-mom6-cobalt2_obsclim_zooc-vint_15arcmin_global_monthly_1961_2010.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zooc_vint = netcdf.getVar(ncid,n-1);
    
end
netcdf.close(ncid);

%% mean over time
zooc_vint(zooc_vint>1e19) = nan;
mzooc_vint = double(mean(zooc_vint,3));

min(zooc_vint(:));
max(zooc_vint(:));
quantile(zooc_vint(:),[0.0 1.0])

%%
save([fpath 'gfdl-mom6-cobalt2_obsclim_zooc-vint_15arcmin_global_monthly_1961_2010.mat'],...
    'zooc_vint','time',...
    'long_name','standard_name','units','lat','lon');

%% Map
[gLAT,gLON] = meshgrid(lat,lon);

clatlim=[-90 90];
clonlim=[-280 80];

mlonlim=[-180 180];

load coastlines;

%% MAPS

figure(1)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,mzooc_vint)
%cmocean('matter')
colormap('jet')
caxis([0.009 0.28])
colorbar
title('zooc-vint')
%text(-2.5,2.25,num2str(qint_poc),'HorizontalAlignment','left')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_GFDL_15arcmin_jet_zoocvint.png'])




