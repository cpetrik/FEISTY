%% FEISTY Make file

clear all
close all

%%%%!! EXPERIMENTS
spinup_ipsl = false;
spinup_gfdl = false;
pre_industrial_ipsl = false;
pre_industrial_gfdl = false;
historic_ipsl = false;
historic_gfdl = false;
ssp126_ipsl = true;
ssp126_gfdl = false;
ssp585_ipsl = true;
ssp585_gfdl = false;

forecast_cesm = false;
forecast_gfdl = false;
temp_cont_cesm = false;
temp_cont_gfdl = false;
npp_cont_cesm = false;
npp_cont_gfdl = false;

tic
if spinup_ipsl
%     Spinup_pristine_ipsl()
    Spinup_pristine_empHP_ipsl()
end
if spinup_gfdl
%     Spinup_pristine_gfdl()
    Spinup_pristine_empHP_gfdl()
end
if historic_ipsl
%     Historic_pristine_ipsl()
    Historic_pristine_empHP_ipsl()
end
if historic_gfdl
%     Historic_pristine_gfdl()
    Historic_pristine_empHP_gfdl()
end
if ssp126_ipsl
%     SSP126_pristine_ipsl()
    SSP126_pristine_empHP_ipsl()
end
if ssp126_gfdl
%     SSP126_pristine_gfdl()
    SSP126_pristine_empHP_gfdl()
end
if ssp585_ipsl
%     SSP585_pristine_ipsl()
    SSP585_pristine_empHP_ipsl()
end
if ssp585_gfdl
%     SSP585_pristine_gfdl()
    SSP585_pristine_empHP_gfdl()
end
if pre_industrial_ipsl
%     Preindust_pristine_ipsl()
    Preindust_pristine_empHP_ipsl()
end
if pre_industrial_gfdl
%     Preindust_pristine_gfdl()
    Preindust_pristine_empHP_gfdl()
end


if forecast_cesm
    Forecast_cesm_noD()
end
if forecast_gfdl
    Forecast_gfdl()
end
if temp_cont_cesm
    Temp_cont_ipsl_noD()
end
if temp_cont_gfdl
    Temp_cont_gfdl()
end
if npp_cont_cesm
    NPP_cont_ipsl_noD()
end
if npp_cont_gfdl
    NPP_cont_gfdl()
end
toc
