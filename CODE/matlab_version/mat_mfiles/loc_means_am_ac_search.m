% FEISTY output at all locations
% Only last 12 months of 150 years saved 

clear all
close all

cfile = 'Dc_enc-b200_m4-b175-k086_c-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught','fmort'};
cols=cols';
spots=spots';
spots{16} = 'PUP';

aep = 5:5:60;
acp = 5:1:19;

%%
%needed for residence time: bio, nu, gamma, rep, rec, pred, nmort, fmort 
SF = NaN*ones(length(aep),length(acp),16); 
SP = SF; %bio
SD = SF;
MF = SF;
MP = SF;
MD = SF;
LP = SF;
LD = SF;
BI = SF;
nSF = SF; %nu
nSP = SF;
nSD = SF;
nMF = SF;
nMP = SF;
nMD = SF;
nLP = SF;
nLD = SF;
gSF = SF; %gamma
gSP = SF;
gSD = SF;
gMF = SF;
gMP = SF;
gMD = SF;
gLP = SF;
gLD = SF;
eMF = SF; %rep
eLP = SF;
eLD = SF;
rSF = SF; %rec
rSP = SF;
rSD = SF;
rMF = SF;
rMP = SF;
rMD = SF;
rLP = SF;
rLD = SF;
pSF = SF; %pred
pSP = SF;
pSD = SF;
pMF = SF;
pMP = SF;
pMD = SF;
pLP = SF;
pLD = SF;
mSF = SF; %nmort
mSP = SF;
mSD = SF;
mMF = SF;
mMP = SF;
mMD = SF;
mLP = SF;
mLD = SF;
fMF = SF; %fmort
fMP = SF;
fMD = SF;
fLP = SF;
fLD = SF;

%%
for k=1:length(acp)
    for j=1:length(aep)   
        h = acp(k);         % coeff on Cmax
        gam = aep(j);       % coeff on search area

    
    sfile = ['/Volumes/MIP/NC/Matlab_new_size/param_ensemble/',cfile,...
        '/Climatol_All_fish030_enc',num2str(gam),'_cmax',num2str(h),'_locs.mat'];
    load(sfile);
    
    % Last year
    [nt,nv,id] = size(S_Cobalt);
    time=1:nt;
    
    %biomass = 1
    SF(j,k,:)=squeeze(mean(S_Sml_f(:,1,:),1));
    SP(j,k,:)=squeeze(mean(S_Sml_p(:,1,:),1));
    SD(j,k,:)=squeeze(mean(S_Sml_d(:,1,:),1));
    MF(j,k,:)=squeeze(mean(S_Med_f(:,1,:),1));
    MP(j,k,:)=squeeze(mean(S_Med_p(:,1,:),1));
    MD(j,k,:)=squeeze(mean(S_Med_d(:,1,:),1));
    LP(j,k,:)=squeeze(mean(S_Lrg_p(:,1,:),1));
    LD(j,k,:)=squeeze(mean(S_Lrg_d(:,1,:),1));
    BI(j,k,:)=squeeze(mean(S_Cobalt(:,1,:),1));
    
    %nu = 15
    nSF(j,k,:)=squeeze(mean(S_Sml_f(:,15,:),1));
    nSP(j,k,:)=squeeze(mean(S_Sml_p(:,15,:),1));
    nSD(j,k,:)=squeeze(mean(S_Sml_d(:,15,:),1));
    nMF(j,k,:)=squeeze(mean(S_Med_f(:,15,:),1));
    nMP(j,k,:)=squeeze(mean(S_Med_p(:,15,:),1));
    nMD(j,k,:)=squeeze(mean(S_Med_d(:,15,:),1));
    nLP(j,k,:)=squeeze(mean(S_Lrg_p(:,15,:),1));
    nLD(j,k,:)=squeeze(mean(S_Lrg_d(:,15,:),1));
    
    %gamma = 16
    gSF(j,k,:)=squeeze(mean(S_Sml_f(:,16,:),1));
    gSP(j,k,:)=squeeze(mean(S_Sml_p(:,16,:),1));
    gSD(j,k,:)=squeeze(mean(S_Sml_d(:,16,:),1));
    gMF(j,k,:)=squeeze(mean(S_Med_f(:,16,:),1));
    gMP(j,k,:)=squeeze(mean(S_Med_p(:,16,:),1));
    gMD(j,k,:)=squeeze(mean(S_Med_d(:,16,:),1));
    gLP(j,k,:)=squeeze(mean(S_Lrg_p(:,16,:),1));
    gLD(j,k,:)=squeeze(mean(S_Lrg_d(:,16,:),1));
    
    %rep = 18 
    eMF(j,k,:)=squeeze(mean(S_Med_f(:,18,:),1));
    eLP(j,k,:)=squeeze(mean(S_Lrg_p(:,18,:),1));
    eLD(j,k,:)=squeeze(mean(S_Lrg_d(:,18,:),1));
    
    %rec = 19
    rSF(j,k,:)=squeeze(mean(S_Sml_f(:,19,:),1));
    rSP(j,k,:)=squeeze(mean(S_Sml_p(:,19,:),1));
    rSD(j,k,:)=squeeze(mean(S_Sml_d(:,19,:),1));
    rMF(j,k,:)=squeeze(mean(S_Med_f(:,19,:),1));
    rMP(j,k,:)=squeeze(mean(S_Med_p(:,19,:),1));
    rMD(j,k,:)=squeeze(mean(S_Med_d(:,19,:),1));
    rLP(j,k,:)=squeeze(mean(S_Lrg_p(:,19,:),1));
    rLD(j,k,:)=squeeze(mean(S_Lrg_d(:,19,:),1));
    
    %pred = 22 
    pSF(j,k,:)=squeeze(mean(S_Sml_f(:,22,:),1));
    pSP(j,k,:)=squeeze(mean(S_Sml_p(:,22,:),1));
    pSD(j,k,:)=squeeze(mean(S_Sml_d(:,22,:),1));
    pMF(j,k,:)=squeeze(mean(S_Med_f(:,22,:),1));
    pMP(j,k,:)=squeeze(mean(S_Med_p(:,22,:),1));
    pMD(j,k,:)=squeeze(mean(S_Med_d(:,22,:),1));
    pLP(j,k,:)=squeeze(mean(S_Lrg_p(:,22,:),1));
    pLD(j,k,:)=squeeze(mean(S_Lrg_d(:,22,:),1));
    
    %nmort = 23
    mSF(j,k,:)=squeeze(mean(S_Sml_f(:,23,:),1));
    mSP(j,k,:)=squeeze(mean(S_Sml_p(:,23,:),1));
    mSD(j,k,:)=squeeze(mean(S_Sml_d(:,23,:),1));
    mMF(j,k,:)=squeeze(mean(S_Med_f(:,23,:),1));
    mMP(j,k,:)=squeeze(mean(S_Med_p(:,23,:),1));
    mMD(j,k,:)=squeeze(mean(S_Med_d(:,23,:),1));
    mLP(j,k,:)=squeeze(mean(S_Lrg_p(:,23,:),1));
    mLD(j,k,:)=squeeze(mean(S_Lrg_d(:,23,:),1));
    
    %fmort = 26
    fMF(j,k,:)=squeeze(mean(S_Med_f(:,26,:),1));
    fMP(j,k,:)=squeeze(mean(S_Med_p(:,26,:),1));
    fMD(j,k,:)=squeeze(mean(S_Med_d(:,26,:),1));
    fLP(j,k,:)=squeeze(mean(S_Lrg_p(:,26,:),1));
    fLD(j,k,:)=squeeze(mean(S_Lrg_d(:,26,:),1));
    
    end
end

%%
nfile = ['/Volumes/MIP/NC/Matlab_new_size/param_ensemble/',cfile,'/'];
save([nfile 'Locs_Climatol_All_fish03_means_aenc_acmax_search.mat'])






