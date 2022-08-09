% Read GFDL 1/4 netcdfs
% Total organic nitrogen vertical integral
% Raw units provided by GFDL
% Before conversion by Matthias

clear all
close all

%fpath='/Volumes/MIP/Fish-MIP/Phase3/GFDL_reanalysis/';
fpath='/Volumes/petrik-lab/Fish-MIP/Phase3/GFDL_reanalysis/';

%% one file
ncdisp([fpath '20000101.ocean_n_budget_FishMIP.nc'])

%%
% Global Attributes:
% filename         = '20000101.ocean_cobalt_n_budget.nc'
% title            = 'OM4p25_LM3TAN_jra_from1959_init234yrs_dynrivernut_O1756_new'
% associated_files = 'areacello: 20000101.ocean_static.nc'
% grid_type        = 'regular'
% grid_tile        = 'N/A'
% history          = 'Thu Mar 10 09:35:03 2022: ncks -v wc_vert_int_on 20000101.ocean_cobalt_n_budget.nc 20000101.ocean_cobalt_n_budget_FishMIP.nc'
% NCO              = '"4.5.4"'

% Dimensions:
% time = 12    (UNLIMITED)
% yh   = 1080
% xh   = 1440


% time
% Size:       12x1
% Dimensions: time
% Datatype:   double
% Attributes:
% long_name      = 'time'
% units          = 'days since 1959-01-01 00:00:00'
% cartesian_axis = 'T'
% calendar_type  = 'JULIAN'
% calendar       = 'JULIAN'
% bounds         = 'time_bnds'

% wc_vert_int_on
% Size:       1440x1080x12
% Dimensions: xh,yh,time
% Datatype:   single
% Attributes:
% long_name     = 'Total organic nitrogen vertical integral'
% units         = 'mol m-2'
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20
% cell_methods  = 'area:mean yh:mean xh:mean time: mean'
% cell_measures = 'area: areacello'
% time_avg_info = 'average_T1,average_T2,average_DT'

% xh
% Size:       1440x1
% Dimensions: xh
% Datatype:   double
% Attributes:
% long_name      = 'h point nominal longitude'
% units          = 'degrees_east'
% cartesian_axis = 'X'

% yh
% Size:       1080x1
% Dimensions: yh
% Datatype:   double
% Attributes:
% long_name      = 'h point nominal latitude'
% units          = 'degrees_north'
% cartesian_axis = 'Y'


%%
ncid = netcdf.open([fpath '20000101.ocean_n_budget_FishMIP.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end


%% Time




% save([fpath 'gfdl_hist_zmeso_100_monthly_1950_2014.mat'],'zmeso_100','yr',...
%     'long_name','standard_name','units_orig','units_vint','lat','lon',...
%     'runs','z100','lev');
%
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
print('-dpng',[ppath 'Map_GFDL_CMIP6_1deg_from_ESGF_jet_intpoc.png'])



