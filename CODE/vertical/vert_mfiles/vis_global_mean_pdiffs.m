% Plot percent differences between COBALT only and COBALT-FEISTY 
% COBALT all vars

clear
close all

%%
npath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FEISTYon/';
fpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/COBALTonly/';

load([fpath 'grid_OM4_05_COBALTv3.mat'],'wet',...
    'z_l_units','z_l_long_name','z_l','geolon','geolat','areacello')

dz = diff(z_l);

[ni,nj] = size(geolon);

%% ONLINE -----------------------------------------------------------
%npath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/';

cfile ='NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/',...
    cfile,'/Cobalt_Feisty/'];

exper = 'OM4_05_COBALTv3_FEISTYon_021326';

%% 
load([npath 'ocean_cobalt_nuts_month_z.199001-199412_means.mat'])

%% Chl from ug/kg to ug/m3
matarea = reshape(areacello,ni*nj,1);

sNO3 = reshape(sNO3,ni*nj,6);
sNH4 = reshape(sNH4,ni*nj,6);
sO2 = reshape(sO2,ni*nj,6);

NmNH4 = sum(sNH4(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');
NmNO3 = sum(sNO3(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');
NmO2 = sum(sO2(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');

clear tNH4 tNO3 vNH4 vNO3 mNH4 mNO3
clear tO2 tCHL vO2 vCHL mO2 mCHL mCHLs

%% Put everything in gC m^-^2 
load([npath 'ocean_cobalt_tracers_month_z.199001-199412_means.mat'])

sSP = reshape(sSP,ni*nj,6);
sLP = reshape(sLP,ni*nj,6);
sSZ = reshape(sSZ,ni*nj,6);
sDI = reshape(sDI,ni*nj,6);
sDE = reshape(sDE,ni*nj,6);

NmDE = sum(sDE(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');
NmDI = sum(sDI(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');
NmSP = sum(sSP(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');
NmLP = sum(sLP(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');
NmSZ = sum(sSZ(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');

clear vDE vDI mDE mDI
clear vSP vLP vSZ mSP mLP mSZ

%% Put everything in gC m^-^2 or gC m^-^2 d^-^1
load([npath 'ocean_cobalt_feisty_forcing_z.199001-199912_means.mat'])

sMZ = reshape(sMZ,ni*nj,6);
sLZ = reshape(sLZ,ni*nj,6);
sBDE = reshape(sDE,ni*nj,6);

NmBDE = sum(sBDE(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');
NmMZ = sum(sMZ(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');
NmLZ = sum(sLZ(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');

clear mDE mMZ mLZ 
clear vMZ vLZ 

%% OFFLINE -----------------------------------------------------------
%fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

mod = 'OM4_05_COBALTv3_FEISTYoff';

load([fpath 'ocean_cobalt_nuts_month_z.199001-199412_means.mat'])

%% Chl from ug/kg to ug/m3
matarea = reshape(areacello,ni*nj,1);

sNO3 = reshape(sNO3,ni*nj,6);
sNH4 = reshape(sNH4,ni*nj,6);
sO2 = reshape(sO2,ni*nj,6);

FmNH4 = sum(sNH4(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');
FmNO3 = sum(sNO3(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');
FmO2 = sum(sO2(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');

clear tNH4 tNO3 vNH4 vNO3 mNH4 mNO3
clear tO2 tCHL vO2 vCHL mO2 mCHL mCHLs

%% Put everything in gC m^-^2 
load([fpath 'ocean_cobalt_tracers_month_z.199001-199412_means.mat'])

sSP = reshape(sSP,ni*nj,6);
sLP = reshape(sLP,ni*nj,6);
sSZ = reshape(sSZ,ni*nj,6);
sDI = reshape(sDI,ni*nj,6);
sDE = reshape(sDE,ni*nj,6);

FmDE = sum(sDE(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');
FmDI = sum(sDI(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');
FmSP = sum(sSP(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');
FmLP = sum(sLP(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');
FmSZ = sum(sSZ(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');

clear vDE vDI mDE mDI
clear vSP vLP vSZ mSP mLP mSZ

%% Put everything in gC m^-^2 or gC m^-^2 d^-^1
load([fpath 'ocean_cobalt_feisty_forcing_z.199001-199912_means.mat'])

sMZ = reshape(sMZ,ni*nj,6);
sLZ = reshape(sLZ,ni*nj,6);
sBDE = reshape(sDE,ni*nj,6);

FmBDE = sum(sBDE(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');
FmMZ = sum(sMZ(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');
FmLZ = sum(sLZ(:,1) .* matarea,'omitnan') ./ sum(matarea,'omitnan');

clear mDE mMZ mLZ 
clear vMZ vLZ 

%% Percent Diffs Online - Offline / Offline
ltex = {'O2','POC','NO3','NH4','Diaz','SP','LP','SZ','MZ','LZ'};

mglob(1) = 100*(NmO2 - FmO2) ./ FmO2;
mglob(2) = 100* (NmDE - FmDE) ./ FmDE;
mglob(3) = 100*(NmNO3 - FmNO3) ./ FmNO3;
mglob(4) = 100*(NmNH4 - FmNH4) ./ FmNH4;
mglob(5) = 100* (NmDI - FmDI) ./ FmDI;
mglob(6) = 100* (NmSP - FmSP) ./ FmSP;
mglob(7) = 100* (NmLP - FmLP) ./ FmLP;
mglob(8) = 100* (NmSZ - FmSZ) ./ FmSZ;
mglob(9) = 100*(NmMZ - FmMZ) ./ FmMZ;
mglob(10) = 100*(NmLZ - FmLZ) ./ FmLZ;

dmBDE = 100*(NmBDE - FmBDE) ./ FmBDE;

save([fpath 'ocean_cobalt_feisty_all.199001-199412_global_Pdiff.mat'],...
    'mglob','dmBDE','ltex')
save([npath 'ocean_cobalt_feisty_all.199001-199412_global_Pdiff.mat'],...
    'mglob','dmBDE','ltex')

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
f2 = figure('Units','inches','Position',[1 3 7 8]);
barh(mglob,'FaceColor',[0 0.5 0.75])
set(gca,'YTickLabel',ltex)
title('% Change Globally')
xlim([-210 210])

print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_mean_Pdiff_all_global.png'])

%%
f3 = figure('Units','inches','Position',[1 3 7 8]);
barh(mglob(2:10),'FaceColor',[0 0.5 0.75])
set(gca,'YTickLabel',ltex(2:10))
title('% Change Globally')
xlim([-210 210])

print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_mean_Pdiff_all_global_noO2.png'])
