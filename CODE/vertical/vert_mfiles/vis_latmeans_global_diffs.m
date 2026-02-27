% nuts, LTL, and feisty forcing
% means by latitude

clear
close all

%% molN/kg to gC/m3 or gC/m3/d
S2D= 60*60*24;
D2Y = S2D*365;
%MZ and LZ molN/kg
N2Ckg= 1035 * (106/16) * 12.01;
%HPloss 'mol N kg-1 s-1'
N2CkgD = N2Ckg * S2D;
N2CkgY = N2Ckg * D2Y;
%fn_tot_btm = 'mol m-2 s-1'
N2CmD = (106/16) * 12.01 * S2D;
N2CmY = (106/16) * 12.01 * D2Y;

% molN/kg to gC/m3
NtoC= 1035 * (106/16) * 12.01;

%% ONLINE -----------------------------------------------------------
%npath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/';
npath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FEISTYon/';

cfile ='NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/',...
    cfile,'/Cobalt_Feisty/'];

exper = 'OM4_05_COBALTv3_FEISTYon_021326';

%% 
load([npath 'ocean_cobalt_nuts_month_z.199001-199412_means_lat.mat'])

% Chl from ug/kg to ug/m3

NvNH4 = vNH4;
NvNO3 = vNO3;
NvO2  = vO2;
NvCHL = vCHL*1035;

NmNH4 = mNH4;
NmNO3 = mNO3;
NmO2 = mO2;
NmCHL = mCHL*1035;
NmCHLs = mCHLs*1035;

clear tNH4 tNO3 vNH4 vNO3 mNH4 mNO3
clear tO2 tCHL vO2 vCHL mO2 mCHL mCHLs

%% Put everything in gC m^-^2 
load([npath 'ocean_cobalt_tracers_month_z.199001-199412_means_lat.mat'])

NvDE = NtoC * vDE;
NvDI = NtoC * vDI;
NvSP = NtoC * vSP;
NvLP = NtoC * vLP;
NvSZ = NtoC * vSZ;

NmDE = NtoC * mDE;
NmDI = NtoC * mDI;
NmSP = NtoC * mSP;
NmLP = NtoC * mLP;
NmSZ = NtoC * mSZ;

clear vDE vDI mDE mDI
clear vSP vLP vSZ mSP mLP mSZ

%% Put everything in gC m^-^2 or gC m^-^2 d^-^1
load([npath 'ocean_cobalt_feisty_forcing_z.199001-199412_means_lat.mat'])

NmBDE = N2CmD * mDE;
NmMZ = N2Ckg * mMZ;
NmLZ = N2Ckg * mLZ;

NvMZ = N2Ckg * vMZ;
NvLZ = N2Ckg * vLZ;

clear mDE mMZ mLZ 
clear vMZ vLZ 

%% OFFLINE -----------------------------------------------------------
%fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
fpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/COBALTonly/';

mod = 'OM4_05_COBALTv3_FEISTYoff';

%% Put everything in gC m^-^2 
load([fpath 'ocean_cobalt_nuts_month_z.199001-199412_means_lat.mat'])

FvNH4 = vNH4;
FvNO3 = vNO3;
FvO2  = vO2;
FvCHL = vCHL*1035;

FmNH4 = mNH4;
FmNO3 = mNO3;
FmO2  = mO2;
FmCHL = mCHL*1035;
FmCHLs = mCHLs*1035;

clear tNH4 tNO3 vNH4 vNO3 mNH4 mNO3
clear tO2 tCHL vO2 vCHL mO2 mCHL mCHLs

FmNO3(FmNO3<0) = 0;
FmO2(FmO2<0) = 0;
NmO2(NmO2<0) = 0;

%% Put everything in gC m^-^2 or gC m^-^2 d^-^1
load([fpath 'ocean_cobalt_feisty_forcing_z.199001-199412_means_lat.mat'])

FmBDE = N2CmD * mDE;
FmMZ = N2Ckg * mMZ;
FmLZ = N2Ckg * mLZ;

FvMZ = N2Ckg * vMZ;
FvLZ = N2Ckg * vLZ;

clear mDE mMZ mLZ 
clear vMZ vLZ 

%% Put everything in gC m^-^2 
load([fpath 'ocean_cobalt_tracers_month_z.199001-199412_means_lat.mat'])

FvDE = NtoC * vDE;
FvDI = NtoC * vDI;
FvSP = NtoC * vSP;
FvLP = NtoC * vLP;
FvSZ = NtoC * vSZ;

FmDE = NtoC * mDE;
FmDI = NtoC * mDI;
FmSP = NtoC * mSP;
FmLP = NtoC * mLP;
FmSZ = NtoC * mSZ;

clear vDE vDI mDE mDI
clear vSP vLP vSZ mSP mLP mSZ


%% Percent Diffs Online - Offline / Offline
dmNH4 = 100*(NmNH4 - FmNH4) ./ FmNH4;
dmNO3 = 100*(NmNO3 - FmNO3) ./ FmNO3;
dmO2 = 100*(NmO2 - FmO2) ./ FmO2;
dmCHL = 100*(NmCHL - FmCHL) ./ FmCHL;
dmCHLs = 100*(NmCHLs - FmCHLs) ./ FmCHLs;

dvNH4 = 100*(NvNH4 - FvNH4) ./ FvNH4;
dvNO3 = 100*(NvNO3 - FvNO3) ./ FvNO3;
dvO2 = 100*(NvO2 - FvO2) ./ FvO2;
dvCHL = 100*(NvCHL - FvCHL) ./ FvCHL;

dmDE = 100* (NmDE - FmDE) ./ FmDE;
dmDI = 100* (NmDI - FmDI) ./ FmDI;
dmSP = 100* (NmSP - FmSP) ./ FmSP;
dmLP = 100* (NmLP - FmLP) ./ FmLP;
dmSZ = 100* (NmSZ - FmSZ) ./ FmSZ;

dvDE = 100* (NvDE - FvDE) ./ FvDE;
dvDI = 100* (NvDI - FvDI) ./ FvDI;
dvSP = 100* (NvSP - FvSP) ./ FvSP;
dvLP = 100* (NvLP - FvLP) ./ FvLP;
dvSZ = 100* (NvSZ - FvSZ) ./ FvSZ;

dmBDE = 100*(NmBDE - FmBDE) ./ FmBDE;
dmMZ = 100*(NmMZ - FmMZ) ./ FmMZ;
dmLZ = 100*(NmLZ - FmLZ) ./ FmLZ;

dvMZ = 100*(NvMZ - FvMZ) ./ FvMZ;
dvLZ = 100*(NvLZ - FvLZ) ./ FvLZ;


%%
load([fpath 'grid_OM4_05_COBALTv3.mat'],'wet',...
    'z_l_units','z_l_long_name','z_l','geolon','geolat')

dz = diff(z_l);

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

%% Vert
figure(1)
subplot(2,3,1)
plot((dvNO3(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvNH4(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvDI(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvSP(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvLP(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on; 
plot((dvSZ(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvMZ(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on; 
plot((dvLZ(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvDE(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvO2(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on;
% legend({'NO3','NH4','Diaz','SP','LP','SZ','MZ','LZ','POC','O2'})
% legend('location','southeast')
title('% Change Equatorial')
ylabel('Depth (m)')

subplot(2,3,2)
plot((dvNO3(1:9,2)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvNH4(1:9,2)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvDI(1:9,2)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvSP(1:9,2)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvLP(1:9,2)),-1*z_l(1:9),'LineWidth',2); hold on; 
plot((dvSZ(1:9,2)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvMZ(1:9,2)),-1*z_l(1:9),'LineWidth',2); hold on; 
plot((dvLZ(1:9,2)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvDE(1:9,2)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvO2(1:9,2)),-1*z_l(1:9),'LineWidth',2); hold on;
% legend({'NO3','NH4','Diaz','SP','LP','SZ','MZ','LZ','POC','O2'})
% legend('location','southeast')
title('% Change Subtropical')
ylabel('Depth (m)')

subplot(2,3,4)
plot((dvNO3(1:9,3)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvNH4(1:9,3)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvDI(1:9,3)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvSP(1:9,3)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvLP(1:9,3)),-1*z_l(1:9),'LineWidth',2); hold on; 
plot((dvSZ(1:9,3)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvMZ(1:9,3)),-1*z_l(1:9),'LineWidth',2); hold on; 
plot((dvLZ(1:9,3)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvDE(1:9,3)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvO2(1:9,3)),-1*z_l(1:9),'LineWidth',2); hold on;
% legend({'NO3','NH4','Diaz','SP','LP','SZ','MZ','LZ','POC','O2'})
% legend('location','southeast')
title('% Change Temperate')
ylabel('Depth (m)')

subplot(2,3,5)
plot((dvNO3(1:9,4)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvNH4(1:9,4)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvDI(1:9,4)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvSP(1:9,4)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvLP(1:9,4)),-1*z_l(1:9),'LineWidth',2); hold on; 
plot((dvSZ(1:9,4)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvMZ(1:9,4)),-1*z_l(1:9),'LineWidth',2); hold on; 
plot((dvLZ(1:9,4)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvDE(1:9,4)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvO2(1:9,4)),-1*z_l(1:9),'LineWidth',2); hold on;
% legend({'NO3','NH4','Diaz','SP','LP','SZ','MZ','LZ','POC','O2'})
% legend('location','southeast')
title('% Change Polar')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_vert_Pdiff_all_lats.png'])


figure(2)
plot((dvNO3(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvNH4(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvDI(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvSP(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvLP(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on; 
plot((dvSZ(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvMZ(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on; 
plot((dvLZ(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvDE(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on;
plot((dvO2(1:9,1)),-1*z_l(1:9),'LineWidth',2); hold on;
legend({'NO3','NH4','Diaz','SP','LP','SZ','MZ','LZ','POC','O2'})
legend('location','southeast')
title('% Change Equatorial')
ylabel('Depth (m)')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_vert_Pdiff_all_eq.png'])


%% Bar plots
ltex = {'O2','POC','NO3','NH4','Diaz','SP','LP','SZ','MZ','LZ'};

mlat(1,:) = dmO2;
mlat(2,:) = dmDE;
mlat(3,:) = dmNO3;
mlat(4,:) = dmNH4;
mlat(5,:) = dmDI;
mlat(6,:) = dmSP;
mlat(7,:) = dmLP;
mlat(8,:) = dmSZ;
mlat(9,:) = dmMZ;
mlat(10,:) = dmLZ;

figure(3)
subplot(2,2,1)
barh(mlat(:,1),'FaceColor',[0 0.5 0.75])
set(gca,'YTickLabel',ltex)
title('% Change Equatorial')
xlim([-100 225])

subplot(2,2,2)
barh(mlat(:,2),'FaceColor',[0 0.5 0.75])
set(gca,'YTickLabel',ltex)
title('% Change Subtropics')
xlim([-100 225])

subplot(2,2,3)
barh(mlat(:,3),'FaceColor',[0 0.5 0.75])
set(gca,'YTickLabel',ltex)
title('% Change Temperate')
xlim([-100 225])

subplot(2,2,4)
barh(mlat(:,4),'FaceColor',[0 0.5 0.75])
set(gca,'YTickLabel',ltex)
title('% Change Polar')
xlim([-100 225])
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_mean_Pdiff_all_lats.png'])

