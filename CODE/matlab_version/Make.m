%% POEM Make file

clear all
close all

%%%%!! EXPERIMENTS
testoneloc = false;
testlocs = false;
histlocs = false;
climlocs = false;
oneloc_fore_pristine = false;
spinup_pristine = false;
climatol_loop = false;
climatol_param = false;
climatol_ens = false;
climatol = false;
climatol_crr = false;
climatol_con = false;
climatol_ngdc = false;
pre_industrial = false;
historic_pristine = false;
historic_fished = false;
historic_fished_ens = false;
forecast_pristine = false;
forecast_fished = false;
forecast_fished_ens = true;

tic
if testoneloc
    Testoneloc()
end
if testlocs
    Locs_clim_param_acmax_aenc_search()
end
if histlocs
    Locs_hist()
end
if climlocs
    Locs_clim()
end
if spinup_pristine
    Spinup_pristine()
end
if climatol_loop
    Climatol_pristine_search()
end
if climatol_param
    Climatol_param_sens_v3()
end
if climatol_ens
    Climatol_param_ensemble6()
end
if climatol
    Climatol_pristine()
end
if climatol_con
    Climatol_con_types()
end
if climatol_crr
    Climatol_con_rec_rep()
end
if climatol_ngdc
    Climatol_nu_gam_die_clev()
end
if pre_industrial
    Pre_industrial()
end
if historic_pristine
    Historic_pristine()
end
if historic_fished
    Historic_fished()
end
if historic_fished_ens
    Historic_fished_ensemble6()
end
if forecast_pristine
    Forecast_pristine()
end
if forecast_fished
    Forecast_fished()
end
if forecast_fished_ens
    Forecast_fished_ensem6()
end
toc
