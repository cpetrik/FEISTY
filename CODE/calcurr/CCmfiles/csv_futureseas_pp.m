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
    
    %% Create csv table
    
    %% Hist
    load([fpath 'FutureSeas_Hist_' vers '_' harv '_outputs.mat'],'f_pp','modyrs')

    hyr = length(modyrs);
    
    per = 'historical';
    period = repmat(per,hyr,1);
    
    ESM = repmat(vers2,hyr,1);

    mod = 'feisty';
    model = repmat(mod,hyr,1);

    modV = 'feisty-basic';
    modelV = repmat(modV,hyr,1);

    rund = '03/11/2024';
    rundate = repmat(rund,hyr,1);

    prov = 'colleen_petrik';
    provider = repmat(prov,hyr,1);

    Vnames = {'Year','Period','ESM','Var_of_interest','Units','Mean',...
        'Model','Model_Version','Run_Date','Comment','Provider'};

    % 4. Population prod - NOT CORRECT YET, rerunning
    %f_pp';
    pvar = 'population productivity';
    punit = '--';
    pcom = 'recruit biomass per adult biomass of all forage fish';
    pV = repmat(pvar,hyr,1);
    pU = repmat(punit,hyr,1);
    pC = repmat(pcom,hyr,1);

    popProd = table(modyrs',period,ESM,pV,pU,f_pp',model,modelV,rundate,pC,provider,...
        'VariableNames',Vnames);

    %% Future
    load([fpath 'FutureSeas_Project_' vers '_' harv '_outputs.mat'],'f_pp','modyrs')

    fyr = length(modyrs);
    
    per = 'future    ';
    period = repmat(per,fyr,1);
    
    ESM = repmat(vers2,fyr,1);

    mod = 'feisty';
    model = repmat(mod,fyr,1);

    modV = 'feisty-basic';
    modelV = repmat(modV,fyr,1);

    rund = '03/11/2024';
    rundate = repmat(rund,fyr,1);

    prov = 'colleen_petrik';
    provider = repmat(prov,fyr,1);


    % 4. Population prod - NOT CORRECT YET, rerunning
    %f_pp';
    pvar = 'population productivity';
    punit = '--';
    pcom = 'recruit biomass per adult biomass of all forage fish';
    pV = repmat(pvar,fyr,1);
    pU = repmat(punit,fyr,1);
    pC = repmat(pcom,fyr,1);

    popProd2 = table(modyrs',period,ESM,pV,pU,f_pp',model,modelV,rundate,pC,provider,...
        'VariableNames',Vnames);


    %% vertcat
    popProdB = [popProd;popProd2]; 
    
    %% write
    writetable(popProdB,[fpath 'FutureSeas_' vers '_' harv '_pp.csv'],'Delimiter',',')

end







