% Compare fish offline and online

clear
close all

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

%% monthly means - online
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
MNTH = repmat(MNTH,1,10);
aa = (cumsum(MNTH)+1);
a = [1,aa(1:end-1)]; %start of the month
b = cumsum(MNTH);    %end of the month

mnB = nan(nz,nnt);
mnSF = mnB;
mnSP = mnB;
mnSD = mnB;
mnMF = mnB;
mnMP = mnB;
mnMD = mnB;
mnLP = mnB;
mnLD = mnB;

tnB = nan(1,nnt);
tnSF = tnB;
tnSP = tnB;
tnSD = tnB;
tnMF = tnB;
tnMP = tnB;
tnMD = tnB;
tnLP = tnB;
tnLD = tnB;

for i = 1:length(MNTH)
    mnB(:,i)  = mean(mBE2(:,a(i):b(i)),2);
    mnSF(:,i) = mean(mSF2(:,a(i):b(i)),2);
    mnSP(:,i) = mean(mSP2(:,a(i):b(i)),2);
    mnSD(:,i) = mean(mSD2(:,a(i):b(i)),2);
    mnMF(:,i) = mean(mMF2(:,a(i):b(i)),2);
    mnMP(:,i) = mean(mMP2(:,a(i):b(i)),2);
    mnMD(:,i) = mean(mMD2(:,a(i):b(i)),2);
    mnLP(:,i) = mean(mLP2(:,a(i):b(i)),2);
    mnLD(:,i) = mean(mLD2(:,a(i):b(i)),2);

    tnB(:,i)  = mean(tB2(:,a(i):b(i)),2);
    tnSF(:,i) = mean(tSF2(:,a(i):b(i)),2);
    tnSP(:,i) = mean(tSP2(:,a(i):b(i)),2);
    tnSD(:,i) = mean(tSD2(:,a(i):b(i)),2);
    tnMF(:,i) = mean(tMF2(:,a(i):b(i)),2);
    tnMP(:,i) = mean(tMP2(:,a(i):b(i)),2);
    tnMD(:,i) = mean(tMD2(:,a(i):b(i)),2);
    tnLP(:,i) = mean(tLP2(:,a(i):b(i)),2);
    tnLD(:,i) = mean(tLD2(:,a(i):b(i)),2);
end

%%
% types combined  - offline
tfF = tfSF + tfMF;
tfP = tfSP + tfMP + tfLP;
tfD = tfSD + tfMD + tfLD;

mfF = S_Sml_f + S_Med_f;
mfP = S_Sml_p + S_Med_p + S_Lrg_p;
mfD = S_Sml_d + S_Med_d + S_Lrg_d;

% types combined - online
tnF = tnSF + tnMF;
tnP = tnSP + tnMP + tnLP;
tnD = tnSD + tnMD + tnLD;

mnF = mnSF + mnMF;
mnP = mnSP + mnMP + mnLP;
mnD = mnSD + mnMD + mnLD;

%% Diffs - 10yr time series
tdiffSF = tnSF - tfSF;
tdiffSP = tnSP - tfSP;
tdiffSD = tnSD - tfSD;
tdiffMF = tnMF - tfMF;
tdiffMP = tnMP - tfMP;
tdiffMD = tnMD - tfMD;
tdiffLP = tnLP - tfLP;
tdiffLD = tnLD - tfLD;

tdiffB = tnB - tfB;
tdiffF = tnF - tfF;
tdiffP = tnP - tfP;
tdiffD = tnD - tfD;

%% Diffs - 10yr vert 2D
mdiffSF = mnSF - S_Sml_f;
mdiffSP = mnSP - S_Sml_p;
mdiffSD = mnSD - S_Sml_d;
mdiffMF = mnMF - S_Med_f;
mdiffMP = mnMP - S_Med_p;
mdiffMD = mnMD - S_Med_d;
mdiffLP = mnLP - S_Lrg_p;
mdiffLD = mnLD - S_Lrg_d;

mdiffB = mnB - S_Bent_bio;
mdiffF = mnF - mfF;
mdiffP = mnP - mfP;
mdiffD = mnD - mfD;

%% Mean annual ts & vertical distrib of years 4-8
st = 1:12:nnt;  %start of the yr
en = 12:12:nnt; %end of the yr

afB = nan(12,1);
afSF = afB;
afSP = afB;
afSD = afB;
afMF = afB;
afMP = afB;
afMD = afB;
afLP = afB;
afLD = afB;

anB = nan(12,1);
anSF = anB;
anSP = anB;
anSD = anB;
anMF = anB;
anMP = anB;
anMD = anB;
anLP = anB;
anLD = anB;

for m = 1:12
    i = (36+m):12:96;
    afB(m)  = mean(tfB(:,i),2);
    afSF(m) = mean(tfSF(:,i),2);
    afSP(m) = mean(tfSP(:,i),2);
    afSD(m) = mean(tfSD(:,i),2);
    afMF(m) = mean(tfMF(:,i),2);
    afMP(m) = mean(tfMP(:,i),2);
    afMD(m) = mean(tfMD(:,i),2);
    afLP(m) = mean(tfLP(:,i),2);
    afLD(m) = mean(tfLD(:,i),2);

    anB(m)  = mean(tnB(:,i),2);
    anSF(m) = mean(tnSF(:,i),2);
    anSP(m) = mean(tnSP(:,i),2);
    anSD(m) = mean(tnSD(:,i),2);
    anMF(m) = mean(tnMF(:,i),2);
    anMP(m) = mean(tnMP(:,i),2);
    anMD(m) = mean(tnMD(:,i),2);
    anLP(m) = mean(tnLP(:,i),2);
    anLD(m) = mean(tnLD(:,i),2);
end

vfB  = mean(S_Bent_bio(:,st(4):en(8)),2);
vfSF = mean(S_Sml_f(:,st(4):en(8)),2);
vfSP = mean(S_Sml_p(:,st(4):en(8)),2);
vfSD = mean(S_Sml_d(:,st(4):en(8)),2);
vfMF = mean(S_Med_f(:,st(4):en(8)),2);
vfMP = mean(S_Med_p(:,st(4):en(8)),2);
vfMD = mean(S_Med_d(:,st(4):en(8)),2);
vfLP = mean(S_Lrg_p(:,st(4):en(8)),2);
vfLD = mean(S_Lrg_d(:,st(4):en(8)),2);

vnB  = mean(mnB(:,st(4):en(8)),2);
vnSF = mean(mnSF(:,st(4):en(8)),2);
vnSP = mean(mnSP(:,st(4):en(8)),2);
vnSD = mean(mnSD(:,st(4):en(8)),2);
vnMF = mean(mnMF(:,st(4):en(8)),2);
vnMP = mean(mnMP(:,st(4):en(8)),2);
vnMD = mean(mnMD(:,st(4):en(8)),2);
vnLP = mean(mnLP(:,st(4):en(8)),2);
vnLD = mean(mnLD(:,st(4):en(8)),2);

anF = anSF + anMF;
anP = anSP + anMP + anLP;
anD = anSD + anMD + anLD;

vnF = vnSF + vnMF;
vnP = vnSP + vnMP + vnLP;
vnD = vnSD + vnMD + vnLD;

afF = afSF + afMF;
afP = afSP + afMP + afLP;
afD = afSD + afMD + afLD;

vfF = vfSF + vfMF;
vfP = vfSP + vfMP + vfLP;
vfD = vfSD + vfMD + vfLD;

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

%% Diffs - 1yr vert 2D
vdSF = vnSF - vfSF;
vdSP = vnSP - vfSP;
vdSD = vnSD - vfSD;
vdMF = vnMF - vfMF;
vdMP = vnMP - vfMP;
vdMD = vnMD - vfMD;
vdLP = vnLP - vfLP;
vdLD = vnLD - vfLD;

vdB = vnB - vfB;
vdF = vnF - vfF;
vdP = vnP - vfP;
vdD = vnD - vfD;

vpB = (vnB - vfB) ./ vfB;
vpF = (vnF - vfF) ./ vfB;
vpP = (vnP - vfP) ./ vfB;
vpD = (vnD - vfD) ./ vfB;

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

%% Plots in time
y = 1:12;

% All size classes of all
figure(1)
%plot(y,(tdB),'Linewidth',1); hold on;
plot(y,(tdSF),'Linewidth',1); hold on;
plot(y,(tdMF),'Linewidth',1); hold on;
plot(y,(tdSP),'Linewidth',1); hold on;
plot(y,(tdMP),'Linewidth',1); hold on;
plot(y,(tdLP),'Linewidth',1); hold on;
plot(y,(tdSD),'Linewidth',1); hold on;
plot(y,(tdMD),'Linewidth',1); hold on;
plot(y,(tdLD),'Linewidth',1); hold on;
%legend('B','SF','MF','SP','MP','LP','SD','MD','LD')
legend('SF','MF','SP','MP','LP','SD','MD','LD')
legend('location','eastoutside')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
%ylabel(' Biomass (g m^-^2)')
ylabel('Difference in integrated Biomass (g m^-^2)')
stamp('')
%print('-dpng',[ppath exper '_ts_all_sizes.png'])

% Fn Types
figure(2)
%plot(y,(tdB),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,(tdF),'r','Linewidth',2); hold on;
plot(y,(tdP),'b','Linewidth',2); hold on;
plot(y,(tdD),'k','Linewidth',2); hold on;
%legend('B','F','P','D')
legend('F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-4 0.5])
xlabel('Time (y)')
ylabel('Difference in integrated Biomass (g m^-^2)')
stamp('')
% print('-dpng',[ppath exper '_ts_all_types.png'])

% Fn Types %diff
figure(3)
%plot(y,(tdB),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,100*(tpF),'r','Linewidth',2); hold on;
plot(y,100*(tpP),'b','Linewidth',2); hold on;
plot(y,100*(tpD),'k','Linewidth',2); hold on;
%legend('B','F','P','D')
legend('F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
%ylim([-4 0.5])
xlabel('Time (y)')
ylabel('% Difference in integrated Biomass')
stamp('')
% print('-dpng',[ppath exper '_ts_all_types.png'])

%% Plots with depth

% All size classes of all
figure(4)
%plot((vdB(1:40)),-1*zl(1:40),'Linewidth',1); hold on;
plot((vdSF(1:40)),-1*zl(1:40),'Linewidth',1); hold on;
plot((vdMF(1:40)),-1*zl(1:40),'Linewidth',1); hold on;
plot((vdSP(1:40)),-1*zl(1:40),'Linewidth',1); hold on;
plot((vdMP(1:40)),-1*zl(1:40),'Linewidth',1); hold on;
plot((vdLP(1:40)),-1*zl(1:40),'Linewidth',1); hold on;
plot((vdSD(1:40)),-1*zl(1:40),'Linewidth',1); hold on;
plot((vdMD(1:40)),-1*zl(1:40),'Linewidth',1); hold on;
plot((vdLD(1:40)),-1*zl(1:40),'Linewidth',1); hold on;
%legend('B','SF','MF','SP','MP','LP','SD','MD','LD')
legend('SF','MF','SP','MP','LP','SD','MD','LD')
legend('location','eastoutside')
ylabel('Depth (m)')
xlabel('Difference in integrated Biomass (g m^-^2)')
stamp('')
%print('-dpng',[ppath exper '_ts_all_sizes.png'])

% Fn Types
figure(5)
%plot((vdB),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot((vdF(1:40)),-1*zl(1:40),'r','Linewidth',2); hold on;
plot((vdP(1:40)),-1*zl(1:40),'b','Linewidth',2); hold on;
plot((vdD(1:40)),-1*zl(1:40),'k','Linewidth',2); hold on;
%legend('B','F','P','D')
legend('F','P','D')
legend('location','eastoutside')
ylabel('Depth (m)')
xlabel('Difference in integrated Biomass (g m^-^2)')
stamp('')
% print('-dpng',[ppath exper '_ts_all_types.png'])

% Fn Types %diff
figure(6)
%plot((tdB),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(100*(vpF(1:40)),-1*zl(1:40),'r','Linewidth',2); hold on;
plot(100*(vpP(1:40)),-1*zl(1:40),'b','Linewidth',2); hold on;
plot(100*(vpD(1:40)),-1*zl(1:40),'k','Linewidth',2); hold on;
%legend('B','F','P','D')
legend('F','P','D')
legend('location','eastoutside')
ylabel('Depth (m)')
xlabel('Difference in integrated Biomass (g m^-^2)')
stamp('')
% print('-dpng',[ppath exper '_ts_all_types.png'])

