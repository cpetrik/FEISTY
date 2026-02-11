%% FEISTY Make file

clear 
close all

%%%%!! EXPERIMENTS
testoneloc = false;
testglobal = true;

tic

if testoneloc
    %BATS_1D_spinup()
    OSP_1D_spinup()
    %BATS_1D_spinup_fromIC()
end
if testglobal
    Global_vert_spinup()
end

toc
