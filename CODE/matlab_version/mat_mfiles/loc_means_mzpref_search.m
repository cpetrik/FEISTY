% FEISTY output at all locations
% Only last 12 months of 150 years saved

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_A050_Sm025_nmort1_BE08_noCC_RE00100';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught','fmort'};
cols=cols';
spots=spots';
spots{16} = 'PUP';

mzp = 0.1:0.05:0.95;

%%
%needed for residence time: bio, nu, gamma, rep, rec, pred, nmort, fmort
SF = NaN*ones(length(mzp),16);
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
for k=1:length(mzp)
    mzpref = mzp(k);
    
    sfile = ['/Volumes/MIP/NC/Matlab_new_size/param_ensemble/',cfile,...
        '/Climatol_All_fish03_MZ0',num2str(100*mzpref),'_locs.mat'];
    load(sfile);
    
    % Last year
    [nt,nv,id] = size(S_Cobalt);
    time=1:nt;
    
    %biomass = 1
    SF(k,:)=squeeze(mean(S_Sml_f(:,1,:),1));
    SP(k,:)=squeeze(mean(S_Sml_p(:,1,:),1));
    SD(k,:)=squeeze(mean(S_Sml_d(:,1,:),1));
    MF(k,:)=squeeze(mean(S_Med_f(:,1,:),1));
    MP(k,:)=squeeze(mean(S_Med_p(:,1,:),1));
    MD(k,:)=squeeze(mean(S_Med_d(:,1,:),1));
    LP(k,:)=squeeze(mean(S_Lrg_p(:,1,:),1));
    LD(k,:)=squeeze(mean(S_Lrg_d(:,1,:),1));
    BI(k,:)=squeeze(mean(S_Cobalt(:,1,:),1));
    
    %nu = 15
    nSF(k,:)=squeeze(mean(S_Sml_f(:,15,:),1));
    nSP(k,:)=squeeze(mean(S_Sml_p(:,15,:),1));
    nSD(k,:)=squeeze(mean(S_Sml_d(:,15,:),1));
    nMF(k,:)=squeeze(mean(S_Med_f(:,15,:),1));
    nMP(k,:)=squeeze(mean(S_Med_p(:,15,:),1));
    nMD(k,:)=squeeze(mean(S_Med_d(:,15,:),1));
    nLP(k,:)=squeeze(mean(S_Lrg_p(:,15,:),1));
    nLD(k,:)=squeeze(mean(S_Lrg_d(:,15,:),1));
    
    %gamma = 16
    gSF(k,:)=squeeze(mean(S_Sml_f(:,16,:),1));
    gSP(k,:)=squeeze(mean(S_Sml_p(:,16,:),1));
    gSD(k,:)=squeeze(mean(S_Sml_d(:,16,:),1));
    gMF(k,:)=squeeze(mean(S_Med_f(:,16,:),1));
    gMP(k,:)=squeeze(mean(S_Med_p(:,16,:),1));
    gMD(k,:)=squeeze(mean(S_Med_d(:,16,:),1));
    gLP(k,:)=squeeze(mean(S_Lrg_p(:,16,:),1));
    gLD(k,:)=squeeze(mean(S_Lrg_d(:,16,:),1));
    
    %rep = 18
    eMF(k,:)=squeeze(mean(S_Med_f(:,18,:),1));
    eLP(k,:)=squeeze(mean(S_Lrg_p(:,18,:),1));
    eLD(k,:)=squeeze(mean(S_Lrg_d(:,18,:),1));
    
    %rec = 19
    rSF(k,:)=squeeze(mean(S_Sml_f(:,19,:),1));
    rSP(k,:)=squeeze(mean(S_Sml_p(:,19,:),1));
    rSD(k,:)=squeeze(mean(S_Sml_d(:,19,:),1));
    rMF(k,:)=squeeze(mean(S_Med_f(:,19,:),1));
    rMP(k,:)=squeeze(mean(S_Med_p(:,19,:),1));
    rMD(k,:)=squeeze(mean(S_Med_d(:,19,:),1));
    rLP(k,:)=squeeze(mean(S_Lrg_p(:,19,:),1));
    rLD(k,:)=squeeze(mean(S_Lrg_d(:,19,:),1));
    
    %pred = 22
    pSF(k,:)=squeeze(mean(S_Sml_f(:,22,:),1));
    pSP(k,:)=squeeze(mean(S_Sml_p(:,22,:),1));
    pSD(k,:)=squeeze(mean(S_Sml_d(:,22,:),1));
    pMF(k,:)=squeeze(mean(S_Med_f(:,22,:),1));
    pMP(k,:)=squeeze(mean(S_Med_p(:,22,:),1));
    pMD(k,:)=squeeze(mean(S_Med_d(:,22,:),1));
    pLP(k,:)=squeeze(mean(S_Lrg_p(:,22,:),1));
    pLD(k,:)=squeeze(mean(S_Lrg_d(:,22,:),1));
    
    %nmort = 23
    mSF(k,:)=squeeze(mean(S_Sml_f(:,23,:),1));
    mSP(k,:)=squeeze(mean(S_Sml_p(:,23,:),1));
    mSD(k,:)=squeeze(mean(S_Sml_d(:,23,:),1));
    mMF(k,:)=squeeze(mean(S_Med_f(:,23,:),1));
    mMP(k,:)=squeeze(mean(S_Med_p(:,23,:),1));
    mMD(k,:)=squeeze(mean(S_Med_d(:,23,:),1));
    mLP(k,:)=squeeze(mean(S_Lrg_p(:,23,:),1));
    mLD(k,:)=squeeze(mean(S_Lrg_d(:,23,:),1));
    
    %fmort = 26
    fMF(k,:)=squeeze(mean(S_Med_f(:,26,:),1));
    fMP(k,:)=squeeze(mean(S_Med_p(:,26,:),1));
    fMD(k,:)=squeeze(mean(S_Med_d(:,26,:),1));
    fLP(k,:)=squeeze(mean(S_Lrg_p(:,26,:),1));
    fLD(k,:)=squeeze(mean(S_Lrg_d(:,26,:),1));
    
end

%%
nfile = ['/Volumes/MIP/NC/Matlab_new_size/param_ensemble/',cfile,'/'];
save([nfile 'Locs_Climatol_All_fish03_means_mzpref_search.mat'])






