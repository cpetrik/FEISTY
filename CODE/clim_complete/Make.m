%% FEISTY Make file

clear all
close all

%%%%!! EXPERIMENTS
climatol = true;
pre_industrial = false;
historic = false;
forecast = false;

tic
if climatol
    Climatol()
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
