%% POEM Make file

clear all
close all

%%%%!! EXPERIMENTS
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
historic_fished_ens = true;
forecast_pristine = false;
forecast_fished = false;
forecast_fished_ens = true;

tic
if climatol_loop
    Climatol_pristine_search()
end
if climatol_param
    Climatol_param_sens_v3()
end
if climatol_ens
    Climatol_param_ensemble6_samek()
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
    %Historic_fished_prod()
end
if historic_fished_ens
    Historic_fished_nu_ensemble6()
end
if forecast_pristine
    Forecast_pristine()
end
if forecast_fished
    Forecast_fished()
end
if forecast_fished_ens
    Forecast_fished_nu_ensem6()
end
toc
