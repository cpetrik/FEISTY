% Read CMIP6 netcdfs
% GFDL-ESM4 Hist
% Intpoc
% 1/4 grid from ISIMIP

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/Phase3/QuarterDeg/';
ppath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/Phase3/quarter_degree_vetting/';

%% intpoc zall
ncdisp([fpath 'gfdl-mom6-cobalt2_obsclim_intpoc_15arcmin_global_monthly_1961_2010.nc'])

%%
% Dimensions:
% time = 600   (UNLIMITED)
% lon  = 1440
% lat  = 720
% Variables:

% time
% units         = 'months since 1901-1-1 00:00:00'
% calendar      = '360_day'

% intpoc
% Size:       1440x720x600
% Dimensions: lon,lat,time
standard_name = 'cean_mass_content_of_particulate_organic_matter_expressed_as_carbon';
long_name     = 'Particulate Organic Carbon Content';
units         = 'kg m-2';
FillValue     = 1.000000020040877e+20;
missing_value = 1.000000020040877e+20;
GFDL_variable = 'wc_vert_int_on';

%%
ncid = netcdf.open([fpath 'gfdl-mom6-cobalt2_obsclim_intpoc_15arcmin_global_monthly_1961_2010.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% mean over time
intpoc(intpoc>1e19) = nan;
mintpoc = double(mean(intpoc,3));

min(intpoc(:))
max(intpoc(:))

%%
save([fpath 'gfdl-mom6-cobalt2_obsclim_intpoc_15arcmin_global_monthly_1961_2010.mat'],...
    'intpoc','time',...
    'long_name','units','lat','lon');

%% Map
[gLAT,gLON] = meshgrid(lat,lon);

clatlim=[-90 90];
clonlim=[-280 80];

mlonlim=[-180 180];

load coastlines;

%% MAPS
%1 - Intpoc
figure(1)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,mintpoc)
%cmocean('matter')
colormap('jet')
caxis([6e-5 0.007])
colorbar
title('int-poc')
%text(-2.5,2.25,num2str(qint_poc),'HorizontalAlignment','left')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_GFDL_15arcmin_jet_intpoc.png'])




