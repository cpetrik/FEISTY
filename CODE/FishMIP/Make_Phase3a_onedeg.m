%% FEISTY Make file

clear 
close all

%%%%!! EXPERIMENTS
spinup_onedeg           = false;
pre_ctrlclim_onedeg     = false;
hist_onedeg_obs         = true;
hist_onedeg_ctrl        = true;

tic

if spinup_onedeg
%     Spinup_pristine_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
%     Spinup_fishing_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
%     Spinup_fishingv1_2_gfdl_onedeg_ctrlclim_server()
%     Spinup_fishingv2_gfdl_onedeg_ctrlclim_server()
    Spinup_fishingv3_2_gfdl_onedeg_ctrlclim_server()
end
if pre_ctrlclim_onedeg  
%     PI_pristine_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
%     PI_fishing_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
%     PI_fishingv1_2_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
%     PI_fishingv2_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
    PI_fishingv3_2_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
end
if hist_onedeg_obs 
%     Hist_pristine_empHP_gfdl_mom6_cobalt2_onedeg_obsclim_server()
%     Hist_fishing_empHP_gfdl_mom6_cobalt2_onedeg_obsclim_server()
%     Hist_fishingv1_2_empHP_gfdl_mom6_cobalt2_onedeg_obsclim_server()
%     Hist_fishingv2_empHP_gfdl_mom6_cobalt2_onedeg_obsclim_server()
    Hist_fishingv3_2_empHP_gfdl_mom6_cobalt2_onedeg_obsclim_server()
end
if hist_onedeg_ctrl 
%     Hist_pristine_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
%     Hist_fishing_empHP_gfdl_mom6_cobalt2_onedeg_obsclim_server()
%     Hist_fishingv1_2_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server() 
%     Hist_fishingv2_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
    Hist_fishingv3_2_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
end

toc
