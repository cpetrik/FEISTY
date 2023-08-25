%% FEISTY Make file

clear 
close all

%%%%!! EXPERIMENTS
spinup_onedeg           = false;
pre_ctrlclim_onedeg     = false;
hist_onedeg_obs         = true;

tic

if spinup_onedeg
    Spinup_fishing_assess_gfdl_onedeg_ctrlclim_server()
    Spinup_fishing_effec_gfdl_onedeg_ctrlclim_server()
    Spinup_fishing_nom_gfdl_onedeg_ctrlclim_server()
end
if pre_ctrlclim_onedeg  
    PI_fishing_assess_gfdl_mom6_cobalt2_onedeg_ctrlclim()
    PI_fishing_effec_gfdl_mom6_cobalt2_onedeg_ctrlclim()
    PI_fishing_nom_gfdl_mom6_cobalt2_onedeg_ctrlclim()
end
if hist_onedeg_obs 
    Hist_fishing_assess_gfdl_mom6_cobalt2_onedeg_obsclim()
    Hist_fishing_effec_gfdl_mom6_cobalt2_onedeg_obsclim()
    Hist_fishing_nom_gfdl_mom6_cobalt2_onedeg_obsclim()
end


toc
