%% POEM Make file

clear all
close all

%%%%!! EXPERIMENTS
climatol = false;
historic_fished = true;
forecast_fished = true;

tic
if climatol
    Climatol_pristine()
end
if historic_fished
    Historic_fished_Zloss()
end
if forecast_fished
    Forecast_fished_Zloss()
end
toc
