%% FEISTY Make file

clear 
close all

%%%%!! EXPERIMENTS
spinup_onedeg           = false;
spinup_15arcmin         = false;
pre_ctrlclim_15arcmin   = false;
pre_ctrlclim_onedeg     = false;
hist_ctrlclim_15arcmin  = true;
hist_obsclim_15arcmin   = true;
hist_ctrlclim_onedeg    = false;
hist_obsclim_onedeg     = false;

tic

if spinup_onedeg
    %Spinup_pristine_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim()
    Spinup_fishingv1_2_gfdl_onedeg_ctrlclim_server()
    Spinup_fishingv2_gfdl_onedeg_ctrlclim_server()
end
if spinup_15arcmin
%     Spinup_pristine_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim()
%     Spinup_fishing_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim()
    Spinup_fishingv3_2_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim() %151 hrs
end
if pre_ctrlclim_15arcmin  
%     PI_pristine_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim_server()
    PI_fishingv3_2_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim_server() %48 hrs
end
if pre_ctrlclim_onedeg   
    PI_pristine_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim()
    PI_fishing_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim()
end
if hist_obsclim_15arcmin  
%     Hist_fishing_empHP_gfdl_mom6_cobalt2_15arcmin_obsclim_server()
%     Hist_pristine_empHP_gfdl_mom6_cobalt2_15arcmin_obsclim_server()
    Hist_fishingv3_2_empHP_gfdl_mom6_cobalt2_15arcmin_obsclim()
end
if hist_ctrlclim_15arcmin 
%     Hist_pristine_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim_server()
%     Hist_fishing_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim()
    Hist_fishingv3_2_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim()
end
if hist_obsclim_onedeg   
    Hist_pristine_empHP_gfdl_mom6_cobalt2_onedeg_obsclim()
    Hist_fishing_empHP_gfdl_mom6_cobalt2_onedeg_obsclim()
end
if hist_ctrlclim_onedeg  
    Hist_pristine_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim()
    Hist_fishing_empHP_gfdl_mom6_cobalt2_onedeg_ctrlclim()
end

toc
