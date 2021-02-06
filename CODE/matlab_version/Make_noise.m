%% FEISTY variance exper Make file

clear all
close all

%%%%!! EXPERIMENTS
Ctrl = true;
TP = false;
TB = false;
Det = false;
MZ = false;
LZ = false;
AllPel = false;
AllBtm = false;
AllZ = false;

tic
if Ctrl
    Noise_exper_control()
end
if TP
    Noise_exper_tp()
end
if TB
    Noise_exper_tb()
end
if Det
    Noise_exper_det()
end
if MZ
    Noise_exper_mz()
end
if LZ
    Noise_exper_lz()
end
if AllPel
    Noise_exper_all_pel()
end
if AllBtm
    Noise_exper_all_btm()
end
if AllZ
    Noise_exper_allZ()
end
toc
