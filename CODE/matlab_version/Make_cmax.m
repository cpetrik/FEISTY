%% POEM Make file

clear all
close all

%%%%!! EXPERIMENTS
historic_fished = false;
historic_fished_ens = true;
forecast_fished = false;
forecast_fished_ens = true;

tic
if historic_fished
%     Historic_fished()
%     Historic_fished_prod()
end
if historic_fished_ens
%     Historic_fished_ensemble6()
%     Historic_fished_prod_ensemble6()
    Historic_fished_cmax_ensemble6()
end
if forecast_fished
%     Forecast_fished()
end
if forecast_fished_ens
%     Forecast_fished_ensem6()
    Forecast_fished_cmax_ensem6()
end
toc
