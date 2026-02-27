% Compare fish offline and online

clear
close all

cfile ='NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/',...
    cfile,'/Cobalt_Feisty/'];


%% Grid
vpath = '/Volumes/petrik-lab/Feisty/NC/MOM6-1D/BATS_vert/cobalt_only/';
load([vpath '20040101.ocean_grid_12mo_BATS.mat'],'zl','zl_long_name',...
    'zi')
dz = diff(zi);

%% Offline
cname = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
offp = ['/Volumes/petrik-lab/Feisty/NC/MOM6-1D/BATS_vert/offline_feisty/' cname '/'];

exper = 'BATS_10yr_COBALT2004_v4_HPcap_All_fish03';

load([offp exper '.mat'])

% All vars 75x120, "S_size_fn"

%% Online
onp = '/Volumes/petrik-lab/Feisty/NC/MOM6-1D/BATS_vert/cobalt_feisty/test_BATS_10yr/';

load([onp '20040101.ocean_feisty_biomass_means_10yr.mat'])

% "mFN2" 75 x 3660
% "tFN2" 1 x 3660

%% time integrated sums - offline
[nz,nnt] = size(S_Sml_f);
dz_mat = repmat(dz,1,nnt);

tfSF=sum(S_Sml_f.*dz_mat,1);
tfMF=sum(S_Med_f.*dz_mat,1);

tfSP=sum(S_Sml_p.*dz_mat,1);
tfMP=sum(S_Med_p.*dz_mat,1);
tfLP=sum(S_Lrg_p.*dz_mat,1);

tfSD=sum(S_Sml_d.*dz_mat,1);
tfMD=sum(S_Med_d,1);
tfLD=sum(S_Lrg_d,1);
tfB =sum(S_Bent_bio,1);


%% Mean biomass of years 4-8

afB  = mean(tfB(37:96));
afSF = mean(tfSF(37:96));
afSP = mean(tfSP(37:96));
afSD = mean(tfSD(37:96));
afMF = mean(tfMF(37:96));
afMP = mean(tfMP(37:96));
afMD = mean(tfMD(37:96));
afLP = mean(tfLP(37:96));
afLD = mean(tfLD(37:96));

anB  = mean(tB2(1099:2928));
anSF = mean(tSF2(1099:2928));
anSP = mean(tSP2(1099:2928));
anSD = mean(tSD2(1099:2928));
anMF = mean(tMF2(1099:2928));
anMP = mean(tMP2(1099:2928));
anMD = mean(tMD2(1099:2928));
anLP = mean(tLP2(1099:2928));
anLD = mean(tLD2(1099:2928));


anF = anSF + anMF;
anP = anSP + anMP + anLP;
anD = anSD + anMD + anLD;


afF = afSF + afMF;
afP = afSP + afMP + afLP;
afD = afSD + afMD + afLD;


%% Diffs 1yr ts
tdSF = anSF - afSF;
tdSP = anSP - afSP;
tdSD = anSD - afSD;
tdMF = anMF - afMF;
tdMP = anMP - afMP;
tdMD = anMD - afMD;
tdLP = anLP - afLP;
tdLD = anLD - afLD;

tdB = anB - afB;
tdF = anF - afF;
tdP = anP - afP;
tdD = anD - afD;

tpB = (anB - afB) ./ afB;
tpF = (anF - afF) ./ afF;
tpP = (anP - afP) ./ afP;
tpD = (anD - afD) ./ afD;

fbar(2) = tpF;
fbar(3) = tpP;
fbar(1) = tpD;

%% colors
cm10=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0.5 0.5 0.5; ...    %med grey
    0 0 0];...      %black


set(groot,'defaultAxesColorOrder',cm10);

%% Bar plot

figure(1)
barh(100*fbar,'FaceColor',[0 0.5 0.75])
set(gca,'YTickLabel',{'Demersal','Forage','Lg Pelagic'})
title('% Change at BATS')
xlim([-100 50])
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_mean_Pdiff_fish_BATS.png'])

