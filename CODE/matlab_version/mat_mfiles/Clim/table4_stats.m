% Stats from Table 4 in Petrik et al. 2019
% SAUP, DvD, Stock PNAS
% save as csv files

clear 
close all

spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY_other/SAUP/';
cpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY_other/grid_cobalt/';
dp = '/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/';

%% DvD on grid
load('/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY_other/DvD_PelDem/Colleen_modeledfish_LME.mat')

%% SAUP
load([spath 'SAUP_LME_Catch_top10_Stock.mat']);

sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);


%% Comparison stats
did=[1:61,63];
load(['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'])
did2 = notLELC(notLELC<=63);

diffD = rPD_catch - dlme_Pfrac;
diffS = rPD_catch - sFracPD_grid;

%r
[rall,pall]=corr(FracLP(did),plme_rPDcatch(did));
[rall2,pall2]=corr(FracLP(did2),plme_rPDcatch(did2));
[rPD,pPD]=corr(sFracPD(notLELC),plme_rPDcatch(notLELC));

%root mean square error
o=FracLP(did);
p=plme_rPDcatch(did);
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);

o=FracLP(did2);
p=plme_rPDcatch(did2);
n = length(o);
num=nansum((p-o).^2);
rmse2 = sqrt(num/n);

o=sFracPD(notLELC);
p=plme_rPDcatch(notLELC);
n = length(o);
num=nansum((p-o).^2);
rmsePD = sqrt(num/n);

%Fmed
Fall=10^(median(FracLP(did)-plme_rPDcatch(did)));
Fall2=10^(median(FracLP(did2)-plme_rPDcatch(did2)));
FPD=10^(median(sFracPD(notLELC)-plme_rPDcatch(notLELC)));


% Table
fish_stat(1,1) = rall;
fish_stat(2,1) = rmse;
fish_stat(3,1) = Fall;
fish_stat(1,2) = rall2;
fish_stat(2,2) = rmse2;
fish_stat(3,2) = Fall2;
fish_stat(1,3) = rPD;
fish_stat(2,3) = rmsePD;
fish_stat(3,3) = FPD;

Fstat = array2table(fish_stat,'RowNames',{'r','RMSE','Fmed'},...
    'VariableNames',{'DvDAllLMEs','DvDnoLELC','SAUnoLELC'});
%writetable(Fstat,[dpath 'LME_DvD_SAU_stats_' cfile '.csv'],'Delimiter',',','WriteRowNames',true)
%save([dpath 'LME_DvD_SAU_stats_' cfile '.mat'],'fish_stat')