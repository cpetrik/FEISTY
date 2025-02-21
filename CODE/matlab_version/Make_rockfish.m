%% FEISTY Make file

clear 
close all

%%%%!! EXPERIMENTS
spinup_pristine = false;
climatol        = false;
core_fished     = false;
hindcast        = true;
pre_industrial  = false;
historic_pristine   = false;
historic_fished     = false;
forecast_pristine   = false;
forecast_fished     = false;

tic

if spinup_pristine
    %Spinup_pristine()
    %Spinup_CORE_pristine()
    %Spinup_CORE_fished_obs_data()
    Spinup_fished_gfdl_mom6_cobalt2_15arcmin_2meso_v2()
end
if climatol
%    Climatol_pristine()
    Climatol_fished()
%    Climatol_be2()
end
if core_fished
    %CORE_fished()
    CORE_fished_obs()
    %CORE_rec_rep_nmort()
    %CORE_nu_gam_die()
end
if hindcast
    %Hindcast_fished_gfdl_mom6_cobalt2_15arcmin_2meso_v2()
    Hindcast_fished_gfdl_mom6_cobalt2_15arcmin_2meso_prod()
    %Hindcast_fished_gfdl_mom6_cobalt2_15arcmin_2meso_yield()
end
if pre_industrial
    %Pre_industrial()
    %Pre_industrial_long()
    Pre_industrial_long_nu()
    Pre_industrial_long_gamma()
    Pre_industrial_long_pred()
    Pre_industrial_long_rep()
end
if historic_pristine
    Historic_pristine()
end
if historic_fished
    Historic_fished()
    %Historic_fished_prod()
    %Historic_fished_met()
end
if forecast_pristine
    Forecast_pristine()
end
if forecast_fished
    Forecast_fished()
end

toc
