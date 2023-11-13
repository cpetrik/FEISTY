% Read BOATS tc
% forced with GFDL 1/4 netcdfs
% obsclim

clear 
close all

%fpath='/Volumes/MIP/Fish-MIP/Phase3/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/BOATS/';

%% one file
ncdisp([fpath 'boats_gfdl-mom6-cobalt2_obsclim_histsoc_default_tc_global_monthly_1961_2010.nc'])

%%
% Format: netcdf4_classic

% Global Attributes:
% contact                 = 'Daniele Bianchi <dbianchi@atmos.ucla.edu>, Eric Galbraith <eric.galbraith@mcgill.ca>, Jerome Guiet <jerome.c.guiet@gmail.com>'
% institution             = 'AOS - University of California Los Angeles'
% comment                 = 'Impact model output for ISIMIP3a'
% date_created            = '15/02/23-15:27'
% isimip_qc_pass_date     = '2023-05-22 11:07 UTC'
% isimip_qc_version       = '3.0.2'

% Dimensions:
% lon  = 1440
% lat  = 720
% time = 600   (UNLIMITED)

% Variables:
time_units = 'days since 1901-01-01 00:00:00';
% calendar      = 'gregorian'
% axis          = 'T'
% tcb
% Size:       1440x720x600
% Dimensions: lon,lat,time
% Datatype:   single
% Attributes:
% _FillValue    = 1.000000020040877e+20
% tc
% Size:       1440x720x600
% Dimensions: lon,lat,time
% Datatype:   single
% Attributes:
% _FillValue    = 1.000000020040877e+20
% long_name     = 'Total Catch (all commercial functional groups / size classes)'
tc_units = 'g m-2';
% missing_value = 1e+20


%%
ncid = netcdf.open([fpath 'boats_gfdl-mom6-cobalt2_obsclim_histsoc_default_tc_global_monthly_1961_2010.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.0000000e+20) = NaN;']);
end
netcdf.close(ncid);

%%
tc(tc>1e19) = nan;
tc = double(tc);

[ni,nj,nt]=size(tc);

%% Time
yr = 1901 + (time/365);

%% LME catch totals --------------------------------------------

%% Map data
cpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
load([cpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat']);
load([cpath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'], 'GRD');
load([cpath 'lme_gfdl-mom6-cobalt2_15arcmin.mat'],'tlme');
load([cpath 'cellarea_15arcmin.mat'],'cell_area');

%% Make annual ts for comp with fishing
% g/m2/mo? --> g/mo
tc = tc .* repmat(cell_area,1,1,nt);

%% Sum months each year to get annual total
st = 1:12:nt;
en = 12:12:nt;
nyr = length(st);

Cann=NaN*ones(ni,nj,nyr);

for t=1:nyr
    Cann(:,:,t) = sum(tc(:,:,st(t):en(t)),3,'omitnan');
end

%% Sum LMEs
lme_mcatch = NaN*ones(66,nyr);

for t=1:nyr
    % g/mo? --> total g
    Cyr = Cann(:,:,t);

    for L=1:66
        lid = find(tlme==L);
        %total catch g
        lme_mcatch(L,t) = sum(Cyr(lid),'omitnan');
    end
end

lme_area = NaN*ones(66,1);
for L=1:66
    lid = find(tlme==L);
    %total area of LME
    lme_area(L,1) = sum(cell_area(lid),'omitnan');
end

%%
save([fpath 'boats_gfdl-mom6-cobalt2_obsclim_histsoc_tc_LME_annual_1961_2010.mat'],...
    'lme_mcatch','lme_area');

 
