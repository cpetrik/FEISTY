%% POEM Make file

clear all
close all

%%%%!! EXPERIMENTS
climatol = true;
historic_fished = false;
forecast_fished = false;

tic
if climatol
    Climatol_fished_Zloss()
end
if historic_fished
    Historic_fished_Zloss()
end
if forecast_fished
    Forecast_fished_Zloss()
end
toc
