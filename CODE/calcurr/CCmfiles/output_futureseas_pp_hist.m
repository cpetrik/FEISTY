% FEISTY output for Ensemble analysis

% 4. Population prod = biomass growing/maturing from the juvenile size
% class to the adult size class at time=T divided by the
% biomass of adults at time=(T-1)

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

    %% Hist
    load([fpath 'Means_Hist_' vers '_' harv '_pp_' cfile '.mat'],'ymppMF')

    [nid,nt] = size(ymppMF);
    
    load([Cdir 'Data_grid_nemuro_',vers2,'.mat'],'GRD')

    % Area in m2 instead of km2, matrix
    aream2 = GRD.AREA .* 1e6;
    aream2_mat = repmat(aream2,1,nt);

    %% Area-weighted mean

    % 4. Population prod = biomass growing/maturing from the juvenile size
    % class to the adult size class at time=T divided by the
    % biomass of adults at time=(T-1)
    f_pp = sum( ymppMF .* aream2_mat ,1,"omitnan") ./ sum(aream2_mat);
    
    %% Save to make csv file
    modyrs = 1980:2010;
    save([fpath 'FutureSeas_Hist_' vers '_' harv '_outputs.mat'],...
     'f_pp','-append')



    %% Proj
    load([fpath 'Means_Project_' vers '_' harv '_pp_' cfile '.mat'],'ymppMF')

    [nid,nt] = size(ymppMF);

    % Area in m2 instead of km2, matrix
    aream2 = GRD.AREA .* 1e6;
    aream2_mat = repmat(aream2,1,nt);

    %% Area-weighted mean

    % 4. Population prod = biomass growing/maturing from the juvenile size
    % class to the adult size class at time=T divided by the
    % biomass of adults at time=(T-1)
    f_pp = sum( ymppMF .* aream2_mat ,1,"omitnan") ./ sum(aream2_mat);
    
    %% Save to make csv file
    modyrs = 2011:2100;
    save([fpath 'FutureSeas_Project_' vers '_' harv '_outputs.mat'],...
     'f_pp','-append')

end
