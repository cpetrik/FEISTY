%% POEM Make file

clear all
close all

%%%%!! EXPERIMENTS
testoneloc = false;
testlocs = false;
histlocs = false;
forelocs = false;
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
historic_fished = true;
historic_fished_ens = false;
forecast_pristine = false;
forecast_fished = false;
forecast_fished_ens = false;

tic
if testoneloc
    Testoneloc()
end
if testlocs
    Locs_clim_param_acmax_aenc_search()
end
if forelocs
    Locs_fore()
end
if histlocs
    Locs_hist()
end
if climlocs
    Locs_clim()
end
if spinup_pristine
    %Spinup_pristine()
    Spinup_CORE_pristine()
end
if climatol_loop
    Climatol_fishing_RE_search()
end
if climatol_param
    Climatol_param_sens_v3()
end
if climatol_ens
    Climatol_param_ensemble6_samek()
end
if climatol
%    Climatol_pristine()
    Climatol_fished()
%    Climatol_be2()
end
if climatol_con
    Climatol_con_types()
end
if climatol_crr
    Climatol_con_rec_rep()
end
if climatol_ngdc
    Climatol_nu_gam_die_clev()
    %Climatol_death_vars()
end
if pre_industrial
    Pre_industrial()
end
if historic_pristine
    Historic_pristine()
end
if historic_fished
    %Historic_fished()
    %Historic_fished_prod()
    %CORE_fished()
    CORE_rec_rep_nmort()
    CORE_nu_gam_die()
end
if historic_fished_ens
    Historic_fished_ensemble6_samek()
    Historic_fished_prod_ensemble6_samek()
end
if forecast_pristine
    Forecast_pristine()
end
if forecast_fished
    Forecast_fished()
end
if forecast_fished_ens
    Forecast_fished_ensem6_samek()
end
toc
