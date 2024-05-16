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
    
    %% Hist
    load([fpath 'FutureSeas_Hist_' vers '_' harv '_outputs.mat'])

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

    % 1. Total biomass of sardine
    bvar = 'total biomass';
    bunit = 'gWW';
    bcom = 'total forage fish biomass';
    bV = repmat(bvar,hyr,1);
    bU = repmat(bunit,hyr,1);
    bC = repmat(bcom,hyr,1);

    TotBiom = table(modyrs',period,ESM,bV,bU,f_tbio',model,modelV,rundate,bC,provider,...
        'VariableNames',Vnames);

    % 2. Growth/maturation from one size class to the next = rec
    %f_arec';
    gvar = 'recruitment';
    gunit = 'gWW';
    gcom = 'biomass transitioning from immature to mature stage of forage fish';
    gV = repmat(gvar,hyr,1);
    gU = repmat(gunit,hyr,1);
    gC = repmat(gcom,hyr,1);

    Rec = table(modyrs',period,ESM,gV,gU,f_arec',model,modelV,rundate,gC,provider,...
        'VariableNames',Vnames);

    % 3. Energy available for growth scaled by the biomass = prod
    %f_amprod';
    mvar = 'mean production';
    munit = 'gWW';
    mcom = 'energy available for growth or reproduction mean of forage fish size classes';
    mV = repmat(mvar,hyr,1);
    mU = repmat(munit,hyr,1);
    mC = repmat(mcom,hyr,1);

    mProd = table(modyrs',period,ESM,mV,mU,f_amprod',model,modelV,rundate,mC,provider,...
        'VariableNames',Vnames);

    tvar = 'total production';
    tunit = 'gWW';
    tcom = 'energy available for growth or reproduction sum of forage fish size classes';
    tV = repmat(tvar,hyr,1);
    tU = repmat(tunit,hyr,1);
    tC = repmat(tcom,hyr,1);

    tProd = table(modyrs',period,ESM,tV,tU,f_atprod',model,modelV,rundate,tC,provider,...
        'VariableNames',Vnames);

    % 4. Population prod - NOT CORRECT YET, rerunning
    %f_pp';
%     pvar = 'population productivity';
%     punit = '--';
%     pcom = 'recruit biomass per adult biomass of all forage fish';
%     pV = repmat(pvar,hyr,1);
%     pU = repmat(punit,hyr,1);
%     pC = repmat(pcom,hyr,1);
% 
%     popProd = table(modyrs',period,ESM,pV,pU,f_pp',model,modelV,rundate,pC,provider,...
%         'VariableNames',Vnames);

    % 5. Catch N of 42nd parallel
    %f_Nc';
    nvar = 'catch in North';
    nunit = 'gWW';
    ncom = 'latitude >= 42; all forage fish             ';
    nV = repmat(nvar,hyr,1);
    nU = repmat(nunit,hyr,1);
    nC = repmat(ncom,hyr,1);

    Ncatch = table(modyrs',period,ESM,nV,nU,f_Nc',model,modelV,rundate,nC,provider,...
        'VariableNames',Vnames);

    % 6. Catch S of 42nd parallel
    %f_Sc';
    svar = 'catch in South';
    sunit = 'gWW';
    scom = 'latitude < 42; all forage fish             ';
    sV = repmat(svar,hyr,1);
    sU = repmat(sunit,hyr,1);
    sC = repmat(scom,hyr,1);

    Scatch = table(modyrs',period,ESM,sV,sU,f_Sc',model,modelV,rundate,sC,provider,...
        'VariableNames',Vnames);

    % 7. Annual predation mortality rate
    %f_ampred'; 
    dvar = 'natural mortality';
    dunit = 'gWW';
    dcom = 'annual predation mortality rate of all forage fish';
    dV = repmat(dvar,hyr,1);
    dU = repmat(dunit,hyr,1);
    dC = repmat(dcom,hyr,1);

    Pred = table(modyrs',period,ESM,dV,dU,f_ampred',model,modelV,rundate,dC,provider,...
        'VariableNames',Vnames);

    % 8. Ecosystem prod
    %eco_prod;
    evar = 'ecosystem productivity';
    eunit = '--';
    ecom = 'all forage fish to their prey';
    eV = repmat(evar,hyr,1);
    eU = repmat(eunit,hyr,1);
    eC = repmat(ecom,hyr,1);

    ecoProd = table(modyrs',period,ESM,eV,eU,eco_prod',model,modelV,rundate,eC,provider,...
        'VariableNames',Vnames);
    

    %% Future
    load([fpath 'FutureSeas_Project_' vers '_' harv '_outputs.mat'])

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


    % 1. Total biomass of sardine
    bvar = 'total biomass';
    bunit = 'gWW';
    bcom = 'total forage fish biomass';
    bV = repmat(bvar,fyr,1);
    bU = repmat(bunit,fyr,1);
    bC = repmat(bcom,fyr,1);

    TotBiom2 = table(modyrs',period,ESM,bV,bU,f_tbio',model,modelV,rundate,bC,provider,...
        'VariableNames',Vnames);

    % 2. Growth/maturation from one size class to the next = rec
    %f_arec';
    gvar = 'recruitment';
    gunit = 'gWW';
    gcom = 'biomass transitioning from immature to mature stage of forage fish';
    gV = repmat(gvar,fyr,1);
    gU = repmat(gunit,fyr,1);
    gC = repmat(gcom,fyr,1);

    Rec2 = table(modyrs',period,ESM,gV,gU,f_arec',model,modelV,rundate,gC,provider,...
        'VariableNames',Vnames);

    % 3. Energy available for growth scaled by the biomass = prod
    %f_amprod';
    mvar = 'mean production';
    munit = 'gWW';
    mcom = 'energy available for growth or reproduction mean of forage fish size classes';
    mV = repmat(mvar,fyr,1);
    mU = repmat(munit,fyr,1);
    mC = repmat(mcom,fyr,1);

    mProd2 = table(modyrs',period,ESM,mV,mU,f_amprod',model,modelV,rundate,mC,provider,...
        'VariableNames',Vnames);

    tvar = 'total production';
    tunit = 'gWW';
    tcom = 'energy available for growth or reproduction sum of forage fish size classes';
    tV = repmat(tvar,fyr,1);
    tU = repmat(tunit,fyr,1);
    tC = repmat(tcom,fyr,1);

    tProd2 = table(modyrs',period,ESM,tV,tU,f_atprod',model,modelV,rundate,tC,provider,...
        'VariableNames',Vnames);

    % 4. Population prod - NOT CORRECT YET, rerunning
    %f_pp';
%     pvar = 'population productivity';
%     punit = '--';
%     pcom = 'recruit biomass per adult biomass of all forage fish';
%     pV = repmat(pvar,fyr,1);
%     pU = repmat(punit,fyr,1);
%     pC = repmat(pcom,fyr,1);
% 
%     popProd2 = table(modyrs',period,ESM,pV,pU,f_pp',model,modelV,rundate,pC,provider,...
%         'VariableNames',Vnames);

    % 5. Catch N of 42nd parallel
    %f_Nc';
    nvar = 'catch in North';
    nunit = 'gWW';
    ncom = 'latitude >= 42; all forage fish; 2010 effort';
    nV = repmat(nvar,fyr,1);
    nU = repmat(nunit,fyr,1);
    nC = repmat(ncom,fyr,1);

    Ncatch2 = table(modyrs',period,ESM,nV,nU,f_Nc',model,modelV,rundate,nC,provider,...
        'VariableNames',Vnames);

    % 6. Catch S of 42nd parallel
    %f_Sc';
    svar = 'catch in South';
    sunit = 'gWW';
    scom = 'latitude < 42; all forage fish; 2010 effort';
    sV = repmat(svar,fyr,1);
    sU = repmat(sunit,fyr,1);
    sC = repmat(scom,fyr,1);

    Scatch2 = table(modyrs',period,ESM,sV,sU,f_Sc',model,modelV,rundate,sC,provider,...
        'VariableNames',Vnames);

    % 7. Annual predation mortality rate
    %f_ampred'; 
    dvar = 'natural mortality';
    dunit = 'gWW';
    dcom = 'annual predation mortality rate of all forage fish';
    dV = repmat(dvar,fyr,1);
    dU = repmat(dunit,fyr,1);
    dC = repmat(dcom,fyr,1);

    Pred2 = table(modyrs',period,ESM,dV,dU,f_ampred',model,modelV,rundate,dC,provider,...
        'VariableNames',Vnames);

    % 8. Ecosystem prod
    %eco_prod;
    evar = 'ecosystem productivity';
    eunit = '--';
    ecom = 'all forage fish to their prey';
    eV = repmat(evar,fyr,1);
    eU = repmat(eunit,fyr,1);
    eC = repmat(ecom,fyr,1);

    ecoProd2 = table(modyrs',period,ESM,eV,eU,eco_prod',model,modelV,rundate,eC,provider,...
        'VariableNames',Vnames);

    %% vertcat
    TotBiomB = [TotBiom;TotBiom2];
    RecB     = [Rec;Rec2];
    mProdB   = [mProd;mProd2];
    tProdB   = [tProd;tProd2];
    %popProdB = [popProd;popProd2];
    NcatchB  = [Ncatch;Ncatch2];
    ScatchB  = [Scatch;Scatch2];
    PredB    = [Pred;Pred2];
    ecoProdB = [ecoProd;ecoProd2];

    %% write
    writetable(TotBiomB,[fpath 'FutureSeas_' vers '_' harv '_biom.csv'],'Delimiter',',')
    writetable(RecB,[fpath 'FutureSeas_' vers '_' harv '_rec.csv'],'Delimiter',',')
    writetable(mProdB,[fpath 'FutureSeas_' vers '_' harv '_mprod.csv'],'Delimiter',',')
    writetable(tProdB,[fpath 'FutureSeas_' vers '_' harv '_tprod.csv'],'Delimiter',',')
    %writetable(popProdB,[fpath 'FutureSeas_' vers '_' harv '_pp.csv'],'Delimiter',',')
    writetable(NcatchB,[fpath 'FutureSeas_' vers '_' harv '_Ncatch.csv'],'Delimiter',',')
    writetable(ScatchB,[fpath 'FutureSeas_' vers '_' harv '_Scatch.csv'],'Delimiter',',')
    writetable(PredB,[fpath 'FutureSeas_' vers '_' harv '_pred.csv'],'Delimiter',',')
    writetable(ecoProdB,[fpath 'FutureSeas_' vers '_' harv '_ecoprod.csv'],'Delimiter',',')


end







