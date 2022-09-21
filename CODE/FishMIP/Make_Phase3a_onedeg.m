%% FEISTY Make file

clear all
close all

%%%%!! EXPERIMENTS
spinup                  = false;
spinup_onedeg           = false;
pre_ctrlclim_onedeg     = true;
hist_onedeg             = true;

tic

if spinup 
    Spinup_pristine_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
    Spinup_fishing_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
end
if spinup_onedeg
    Spinup_fishingv1_2_gfdl_onedeg_ctrlclim_server()
    Spinup_fishingv2_gfdl_onedeg_ctrlclim_server()
end
if pre_ctrlclim_onedeg   
    PI_fishingv1_2_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
    PI_fishingv2_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
    netcdf_read_pi_ctrlclim_fishing_onedeg_bio
%     PI_pristine_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
end
if hist_onedeg 
    %Hist_pristine_empHP_gfdl_mom6_cobalt2_onedeg_obsclim_server()
    Hist_fishingv1_2_empHP_gfdl_mom6_cobalt2_onedeg_obsclim_server()
    Hist_fishingv2_empHP_gfdl_mom6_cobalt2_onedeg_obsclim_server()
    %Hist_pristine_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
    Hist_fishingv1_2_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server() 
    Hist_fishingv2_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server() 
    netcdf_read_pi_ctrlclim_fishing_15arcmin_bio
end

toc
