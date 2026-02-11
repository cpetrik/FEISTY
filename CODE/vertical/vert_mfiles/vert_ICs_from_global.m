% ICs for vertical from global run

clear 
close all

%%
gpath = '/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_CC80_RE00100/CORE/';

cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CORE-forced/';
load([cpath 'core_grid_360x200_id_locs_area_dep.mat'],'ids','abbrev','T');

cname = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath = ['/Volumes/petrik-lab/Feisty/NC/MOM6-1D/BATS_vert/offline_feisty/' cname '/'];

%%
for s = 1:length(ids)
    load([gpath 'IC_means_1988_no_move_All_fish03_' abbrev{s} '.mat']);
    save([fpath 'IC_means_1988_no_move_All_fish03_' abbrev{s} '.mat'])
end

%% Do in init subroutine instead
