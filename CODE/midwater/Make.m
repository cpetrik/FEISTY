%% FEISTY Make file

clear all
close all

%%%%!! EXPERIMENTS
locs = false;
climatol = true;
pre_industrial = false;
historic = false;
forecast = false;

tic
if locs
    %Locs_clim()
    Locs_clim_fishvec()
end
if climatol
    %Climatol()
    Climatol_fishvec()
end
if pre_industrial
    Pre_industrial()
end
if historic
    Historic_fished()
    Historic_fished_prod()
end
if forecast
    Forecast_fished()
end

toc
