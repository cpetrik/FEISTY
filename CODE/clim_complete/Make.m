%% FEISTY Make file

clear all
close all

%%%%!! EXPERIMENTS
climatol = false;
climatol_loop = false;
climatol_ensem = true;
pre_industrial = false;
historic_pristine = false;
historic_fished = false;
historic_fished_ens = false;
forecast_pristine = false;
forecast_fished = false;
forecast_fished_ens = false;

%tic
if climatol
    Climatol()
end
if climatol_loop
    Climatol_loop()
end
if climatol_ensem
    Climatol_param_ensemble6_samek()
end
if pre_industrial
    Pre_industrial()
end
if historic_pristine
    Historic_pristine()
end
if historic_fished
    Historic_fished()
    Historic_fished_prod()
end
if historic_fished_ens
    Historic_fished_ensemble6()
    Historic_fished_prod_ensemble6()
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
%toc
