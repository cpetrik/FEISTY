% FEISTY output for Ensemble analysis

% 1. Total biomass of sardine
% 2. Growth/maturation from one size class to the next = rec
% 3. Energy available for growth scaled by the biomass = prod
% 4. Population prod = biomass growing/maturing from the juvenile size
% class to the adult size class at time=T divided by the
% biomass of adults at time=(T-1)
% 5. Catch N of 42nd parallel
% 6. Catch S of 42nd parallel
% 7. Annual predation mortality rate
% 8. Ecosystem prod = Forage biomass / Zooplankton biomass

clear
close all

%%
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';
harv = 'All_fishobs';
esm1 = {'IPSL','GFDL','HAD'};
esm2 = {'ipsl','gfdl','hadley'};

%%
for k=1:3

    %%
    vers = esm1{k};
    vers2 = esm2{k};
    fpath=['/Volumes/petrik-lab/Feisty/NC/NEMURO/',cfile,'/',vers,'/'];
    Cdir = ['/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/',vers,'down/'];


    %% Create csv table
    modyrs = 1980:2010;

    % 1. Total biomass of sardine
    TotBiom(:,1) = modyrs';
    TotBiom(:,2) = f_tbio';

    % 2. Growth/maturation from one size class to the next = rec
    Growth(:,1) = modyrs';
    Growth(:,2) = f_arec';

    % 3. Energy available for growth scaled by the biomass = prod
    AvailEn(:,1) = modyrs';
    AvailEn(:,2) = f_amprod';

    % 4. Population prod
    PopProd(:,1) = modyrs';
    PopProd(:,2) = f_pp';

    % 5. Catch N of 42nd parallel
    CatchN(:,1) = modyrs';
    CatchN(:,2) = f_Nc';

    % 6. Catch S of 42nd parallel
    CatchS(:,1) = modyrs';
    CatchS(:,2) = f_Sc';

    % 7. Annual predation mortality rate
    %mean (no spatial component)
    f_dmpred;
    %per day to per year
    f_ampred; %CHECK THIS!!!

    % 8. Ecosystem prod
    EcoProd(:,1) = modyrs';
    EcoProd(:,2) = eco_prod;


end