%% FEISTY Make file

clear all
close all

%%%%!! EXPERIMENTS
spinup                  = false;
pre_ctrlclim_15arcmin   = true;
pre_ctrlclim_onedeg     = false;
hist_ctrlclim_15arcmin  = false;
hist_obsclim_15arcmin   = false;
hist_ctrlclim_onedeg    = false;
hist_obsclim_onedeg     = false;

tic

if spinup 
    Spinup_pristine_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim()
    Spinup_fishing_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim()
end
if pre_ctrlclim_15arcmin  
    PI_pristine_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim()
    %PI_fishing_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim()
end
if pre_ctrlclim_onedeg   
    PI_pristine_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim()
    PI_fishing_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim()
end
if hist_ctrlclim_15arcmin 
    Hist_pristine_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim()
    Hist_fishing_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim()
end
if hist_obsclim_15arcmin  
    Hist_pristine_empHP_gfdl_mom6_cobalt2_15arcmin_obsclim()
    Hist_fishing_empHP_gfdl_mom6_cobalt2_15arcmin_obsclim()
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
