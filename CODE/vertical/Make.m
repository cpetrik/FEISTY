%% FEISTY Make file

clear 
close all

%%%%!! EXPERIMENTS
testoneloc = true;

tic

if testoneloc
    %BATS_1D_spinup()
    BATS_1D_spinup_fromIC()
end


toc
