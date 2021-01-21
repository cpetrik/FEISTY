% FEISTY output at all locations

close all
clear all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/MIP/NC/FishMIP/GFDL_CMIP6/' cfile '/'];

%%
ncdisp([fpath 'feisty_gfdl-esm4_nobasd_historical_nat_default_bd90cm_global_monthly_1950_2014.nc'])

% time  
%            Size:       780x1
%            Dimensions: time
%            Datatype:   double
%            Attributes:
%                        standard_name = 'time'
%                        long_name     = 'time'
%                        units         = 'months since 1601-01-01'
%                        calendar      = '360_day'
%                        axis          = 'T'
                       
%%
ncdisp([fpath 'feisty_gfdl-esm4_nobasd_historical_nat_default_bp30cm_global_monthly_1950_2014.nc'])

% bp30cm
%            Size:       360x180x780
%            Dimensions: lon,lat,time
%            Datatype:   single
%            Attributes:
%                        _FillValue    = 1.000000020040877e+20
%                        missing_value = 1.000000020040877e+20
%                        long_name     = 'Biomass Density of Small Pelagics <30cm'
%                        units         = 'g m-2'

%%
ncdisp([fpath 'feisty_gfdl-esm4_nobasd_historical_nat_default_bp90cm_global_monthly_1950_2014.nc'])

%%
ncdisp([fpath 'feisty_gfdl-esm4_nobasd_historical_nat_default_tcb_global_monthly_1950_2014.nc'])

%%
ncdisp([fpath 'feisty_gfdl-esm4_nobasd_historical_nat_default_tdb_global_monthly_1950_2014.nc'])

%%
ncdisp([fpath 'feisty_gfdl-esm4_nobasd_historical_nat_default_tpb_global_monthly_1950_2014.nc'])

%total consumber biomass in log10 bins tcblog10 = 360x180xMOsx6
%(1g, 10g, 100g, 1kg, 10kg, 100kg)
%Small <1g      (0.001-0.5g; 0.02g mean)    (0.46-3.68cm; mean 1.3cm)
%Med 1-100g     (0.5-250g; 11.2g mean)      (3.68-29.24cm; mean 10.4cm)
%Lrg 1-100kg    (250-125000g; 5600g mean)   (29.24-232.08cm; mean 82.4cm)



