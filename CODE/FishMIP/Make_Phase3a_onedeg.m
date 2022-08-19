%% FEISTY Make file

clear all
close all

%%%%!! EXPERIMENTS
spinup                  = false;
pre_ctrlclim_onedeg     = true;
hist_ctrlclim_onedeg    = false;
hist_obsclim_onedeg     = false;

tic

if spinup 
    Spinup_pristine_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
    Spinup_fishing_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
end
if pre_ctrlclim_onedeg   
    PI_pristine_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
    PI_fishing_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim_server()
end
if hist_ctrlclim_onedeg  
    Hist_pristine_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim()
    Hist_fishing_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim()
end
if hist_obsclim_onedeg   
    Hist_pristine_empHP_gfdl_mom6_cobalt2_onedeg_obsclim()
    Hist_fishing_empHP_gfdl_mom6_cobalt2_onedeg_obsclim()
end

toc
