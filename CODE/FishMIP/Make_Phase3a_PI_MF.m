%% FEISTY Make file

clear all
close all

%%%%!! EXPERIMENTS
pre_ctrlclim_15arcmin   = true;
pre_ctrlclim_onedeg     = false;

tic

if pre_ctrlclim_15arcmin  
    PI_fishing_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim_server_B()
    netcdf_read_pi_ctrlclim_fishing_15arcmin_MF
    netcdf_read_pi_ctrlclim_pristine_15arcmin_bio
end

toc
