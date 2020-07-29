%% FEISTY Make file

clear all
close all

%%%%!! EXPERIMENTS
spinup_cesm = false;
pre_industrial_cesm = false;
pre_industrial_gfdl = false;
historic_cesm = false;
historic_gfdl = true;
forecast_cesm = false;
forecast_gfdl = false;
temp_cont_cesm = false;
temp_cont_gfdl = false;
npp_cont_cesm = false;
npp_cont_gfdl = false;

tic
if spinup_cesm
    Spinup_cesm_noD()
end
if pre_industrial_cesm
    Pre_industrial_cesm_noD()
end
if pre_industrial_gfdl
    Pre_industrial_gfdl()
end
if historic_cesm
    Historic_cesm_noD()
end
if historic_gfdl
    Historic_fished_gfdl_1zoo()
end
if forecast_cesm
    Forecast_cesm_noD()
end
if forecast_gfdl
    Forecast_gfdl()
end
if temp_cont_cesm
    Temp_cont_cesm_noD()
end
if temp_cont_gfdl
    Temp_cont_gfdl()
end
if npp_cont_cesm
    NPP_cont_cesm_noD()
end
if npp_cont_gfdl
    NPP_cont_gfdl()
end
toc
