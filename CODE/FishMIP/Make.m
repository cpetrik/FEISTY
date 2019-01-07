%% POEM Make file

clear all
close all

%%%%!! EXPERIMENTS
testoneloc = false;
testlocs = false;
histlocs = false;
climlocs = false;
oneloc_fishing = false;
oneloc_hind_pristine = false;
oneloc_fore_pristine = false;
spinup_pristine = false;
climatol_loop = false;
climatol_IC_day = false;
climatol_IC_biom = false;
climatol = false;
climatol_crr = false;
climatol_con = false;
climatol_ngdc = false;
pre_industrial = false;
historic_pristine = false;
historic_fished = false;
forecast_pristine = false;
forecast_fished = true;

tic
if testoneloc
    Testoneloc()
end
if testlocs
    Testlocs()
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
if climatol_IC_day
    Climatol_ICday_loop()
end
if climatol_IC_biom
    Climatol_ICbiom_loop()
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
if forecast_pristine
    Forecast_pristine()
end
if forecast_fished
    Forecast_fished()
end
toc
