%% Table 7 of TEs by Location and LME compared to Mauread

clear all
close all

datap = '/Volumes/GFDL/NC/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/';

% Locations
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'clim_grid_180x360_id_locs_area_dep.mat'],'ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
spots=spots';

% POEM
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
BE = 0.075;
sname = 'Climatol_';
harv = 'All_fish03';
dpath = [datap char(cfile) '/'];
fpath = [figp char(cfile) '/'];

% Load saved data
%Locs
load([dpath sname harv '_locs_lastyr_effTEs.mat'],'TEeff');
load([dpath 'Locs_TE_clim_fished_',harv,'_' cfile '.mat'],'Tab2');

%TEeff: 'TE_Mb','TE_HTLb','TE_Lb','TE_LTLb',...
%       'TE_Md','TE_HTLd','TE_Ld','TE_LTLd'
TEeff_Md   = real(TEeff(5,:));
TEeff_HTLd = real(TEeff(6,:));
TEeff_Ld   = real(TEeff(7,:));
TEeff_LTLd = real(TEeff(8,:));

TE_Md   = real(TEeff(5,:).^(1/2));
TE_HTLd = real(TEeff(6,:).^(1/3));
TE_Ld   = real(TEeff(7,:).^(1/4));
TE_LTLd = real(TEeff(8,:));

%LMEs
load([dpath 'LME_TEeff_Mauread_comp_' cfile '.mat'],'tab','keep',...
    'mECI','pECI')
% tab(:,1)=keep;
% tab(:,2)=mECI;
% tab(:,3)=pECI;

%% Put into table for ms
%col: EBS, PUP, HOT
%row: TEL, TELTL, TEHTL, LMEhtl, ECI 
locs=[8;19;22;21;20;1;51;49;10;6;13];
spot2=[3;5;6;7;8;10;11;12;13;14;16];
name = {'SS','GS','NS','NwS','BS','EBS','OyaCur','KurCur','Haw','SEUS',...
    'Humb'};
Tab=table(TEeff_Ld(spot2)',TEeff_LTLd(spot2)',TEeff_HTLd(spot2)',...
    pECI(locs),mECI(locs),locs,name',...
    'VariableNames',{'TEeffL','TEeffLTL','TEeffHTL',...
    'POEM','Mauread','LME','Name',});
writetable(Tab,[spath 'LME_locs_TEeff_Mauread_',harv,'_' cfile '.csv'],'Delimiter',',');
save([dpath 'LME_locs_TEeff_Mauread_',harv,'_' cfile '.mat'],'Tab');



