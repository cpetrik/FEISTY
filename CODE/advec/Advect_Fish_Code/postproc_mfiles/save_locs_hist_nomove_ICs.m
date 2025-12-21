% Output of FEISTY Historic-CORE at single locations
% 1988-2007, monthly means saved
% Save biomass as IC file

clear
close all

%%
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/CORE/'];
ppath = [pp cfile '/CORE/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

vpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CORE-forced/';
load([vpath 'ocean_cobalt_grid.mat'],'geolon_t','geolat_t');
load([vpath 'Data_grid_ocean_cobalt_ESM2Mcore.mat'],'GRD');
load([vpath 'core_grid_360x200_id_locs_area_dep.mat'],'ids','abbrev','T');

spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
spots=spots';
longnames = T.Location;

nspots = length(spots);

%%
exper = '1988_no_move';
load([fpath 'CORE_Hindcast',exper,'_All_fish03_locs.mat'])


%%
for s=1:nspots
    loc = spots{s};
    lname = [loc '_'];
    loclong = longnames{s};


    %% Save mean for initializing vertical runs

    Sml_f.bio = mean(Mo_Sml_f(:,1,s),1,'omitnan');
    Sml_p.bio = mean(Mo_Sml_p(:,1,s),1,'omitnan');
    Sml_d.bio = mean(Mo_Sml_d(:,1,s),1,'omitnan');
    Med_f.bio = mean(Mo_Med_f(:,1,s),1,'omitnan');
    Med_p.bio = mean(Mo_Med_p(:,1,s),1,'omitnan');
    Med_d.bio = mean(Mo_Med_d(:,1,s),1,'omitnan');
    Lrg_p.bio = mean(Mo_Lrg_p(:,1,s),1,'omitnan');
    Lrg_d.bio = mean(Mo_Lrg_d(:,1,s),1,'omitnan');
    BENT.bio  = mean(Mo_Cobalt(:,1,s),1,'omitnan');

    save([fpath 'IC_means_' exper '_All_fish03_' abbrev{s} '.mat'],'Sml_f','Sml_p','Sml_d',...
        'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')
end
