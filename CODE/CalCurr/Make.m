%% POEM Make file

clear all
close all

%%%%!! EXPERIMENTS
testlocs = false;
histlocs = false;
climlocs = false;
spinup_3 = false;
spinup_10 = false;
climatol = false;
pre_industrial = false;
historic_pristine = false;
historic_fished3_loop = true;
historic_fished3 = false;
historic_fished10 = false;
forecast_pristine = false;
forecast_fished = false;

tic
if testlocs
    Testlocs()
end
if histlocs
    Locs_hist()
end
if climlocs
    Locs_clim()
end
if spinup_3
    Spinup_pristine_3km()
end
if spinup_10
    Spinup_pristine_10km()
end
if climatol
    Climatol_pristine()
end
if pre_industrial
    Pre_industrial()
end
if historic_pristine
    Historic_pristine()
end
if historic_fished3_loop
    Historic_fished_3km_search()
end
if historic_fished3
    Historic_fished_3km()
end
if historic_fished10
    Historic_fished_10km()
end
if forecast_pristine
    Forecast_pristine()
end
if forecast_fished
    Forecast_fished()
end
toc
