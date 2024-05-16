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

    %%
    load([fpath 'Means_Project_' vers '_' harv '_' cfile '.mat'])
    load([fpath 'Means_Project_' vers '_' harv '_prod_' cfile '.mat'])
    load([fpath 'Means_Project_' vers '_' harv '_rec_' cfile '.mat'])
    load([fpath 'Means_Project_' vers '_' harv '_mort_' cfile '.mat'])
    load([fpath 'Means_Project_' vers '_' harv '_catch_' cfile '.mat'])

    load([Cdir 'Data_grid_nemuro_',vers2,'.mat'],'GRD')

    %% Need to do area-weighted means & sums

    [nid,nt] = size(mMF);

    % Area in m2 instead of km2, matrix
    aream2 = GRD.AREA .* 1e6;
    aream2_mat = repmat(aream2,1,nt);

    % 1. Total biomass of sardine [monthly mean -> annual mean, gWW/m2]
    f_tbio = sum( (mSF + mMF) .* aream2_mat ,1,"omitnan");

    % 2. Growth/maturation from one size class to the next = rec
    %total rec into adult class MF [monthly sum -> annual sum, gWW]
    f_arec = sum(ysrMF,1,"omitnan");

    % 3. Energy available for growth scaled by the biomass = prod
    %mean prod of both size classes [g/m2]
    f_mprod = (ympSF + ympMF)/2; %mean of size classes
    f_tprod = (ympSF + ympMF);   %sum of size classes
    %mult by area, then take mean
    f_amprod = mean(f_mprod .* aream2_mat ,1,"omitnan");
    f_atprod = mean(f_tprod .* aream2_mat ,1,"omitnan");

    % 4. Population prod = biomass growing/maturing from the juvenile size
    % class to the adult size class at time=T divided by the
    % biomass of adults at time=(T-1)
    %f_pp = mean(ymppMF,1,"omitnan");

    % 5. Catch N of 42nd parallel [monthly sum -> annual sum, gWW]
    Nid = (GRD.LAT >= 42.0);
    f_Nc = sum(ytcMF(Nid,:),1,"omitnan");

    % 6. Catch S of 42nd parallel [monthly sum -> annual sum, gWW]                  
    Sid = (GRD.LAT < 42.0);
    f_Sc = sum(ytcMF(Sid,:),1,"omitnan");

    % 7. Annual predation mortality 
    %nat mort rate
    Nat_mrt = 0.1;                                  %per y
    f_mnat = Nat_mrt .* (mSF + mMF) .* aream2_mat;  %mean biomass
    %pred mort [monthly sum -> annual sum, gWW]
    f_mpred = (ysdSF + ysdMF);        
    %total mort, per day to per year
    f_ampred = sum( (f_mnat + f_mpred) ,1,"omitnan");

    %% 8. Ecosystem prod = Forage biomass / Zooplankton biomass
    load([Cdir 'feisty_',vers2,'_Tzoo_gWW_1980-2100.mat'],...
         'Zm_proj');

    [~,nmo] = size(Zm_proj);

    %% Year sums
    a = 1:12:nmo; % start of each yr
    b = 12:12:nmo; % end of each yr

    tZ = NaN*ones(nid,nt);

    for i = 1:(nt)
        tZ(:,i) = sum(Zm_proj(:,a(i):b(i)),2,"omitnan");
    end

    tZ00 = tZ .* aream2_mat;
    z_tbio = sum(tZ00,1,"omitnan");

    eco_prod = f_tbio ./ z_tbio;

    %% Save to make csv file
    modyrs = 2011:2100;
    save([fpath 'FutureSeas_Project_' vers '_' harv '_outputs.mat'],...
     'modyrs','f_tbio','f_arec','f_amprod','f_atprod','f_Nc','f_Sc',...
     'f_ampred','eco_prod') %'f_pp',

end
