%% FEISTY Make file

clear all
close all

%%%%!! EXPERIMENTS
spinup                  = false;
pre_ctrlclim_onedeg     = false;
hist_onedeg             = false;
spinup15                = true;
hist_15                 = true;

tic

if spinup 
    Spinup_pristine_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
    Spinup_fishing_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
end
if pre_ctrlclim_onedeg   
    PI_fishing_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
    netcdf_read_pi_ctrlclim_fishing_15arcmin_bio
    PI_pristine_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
end
if hist_onedeg  
    Hist_pristine_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
    Hist_fishing_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server() 
    Hist_pristine_empHP_gfdl_mom6_cobalt2_onedeg_obsclim_server()
    Hist_fishing_empHP_gfdl_mom6_cobalt2_onedeg_obsclim_server()
end
if spinup15 
    PI_fishing_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim_server()
    netcdf_read_pi_ctrlclim_fishing_15arcmin_bio
end
if hist_15  
    Hist_fishing_empHP_gfdl_mom6_cobalt2_15arcmin_obsclim_server()
    Hist_fishing_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim_server() 
    netcdf_read_hist_obsclim_fishing_15arcmin_bio
    netcdf_read_hist_ctrlclim_fishing_15arcmin_bio
    %Hist_pristine_empHP_gfdl_mom6_cobalt2_15arcmin_obsclim_server()
    %Hist_pristine_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
end

toc
