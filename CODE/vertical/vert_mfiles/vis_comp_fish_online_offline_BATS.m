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

tnB = nan(nz,nnt);
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

afB = nan(nz,12);
afSF = afB;
afSP = afB;
afSD = afB;
afMF = afB;
afMP = afB;
afMD = afB;
afLP = afB;
afLD = afB;

anB = nan(nz,12);
anSF = anB;
anSP = anB;
anSD = anB;
anMF = anB;
anMP = anB;
anMD = anB;
anLP = anB;
anLD = anB;

for i = 4:8
    afB(:,i)  = mean(tfB(:,st(i):en(i)),2);
    afSF(:,i) = mean(tfSF(:,st(i):en(i)),2);
    afSP(:,i) = mean(tfSP(:,st(i):en(i)),2);
    afSD(:,i) = mean(tfSD(:,st(i):en(i)),2);
    afMF(:,i) = mean(tfMF(:,st(i):en(i)),2);
    afMP(:,i) = mean(tfMP(:,st(i):en(i)),2);
    afMD(:,i) = mean(tfMD(:,st(i):en(i)),2);
    afLP(:,i) = mean(tfLP(:,st(i):en(i)),2);
    afLD(:,i) = mean(tfLD(:,st(i):en(i)),2);

    anB(:,i)  = mean(tnB(:,st(i):en(i)),2);
    anSF(:,i) = mean(tnSF(:,st(i):en(i)),2);
    anSP(:,i) = mean(tnSP(:,st(i):en(i)),2);
    anSD(:,i) = mean(tnSD(:,st(i):en(i)),2);
    anMF(:,i) = mean(tnMF(:,st(i):en(i)),2);
    anMP(:,i) = mean(tnMP(:,st(i):en(i)),2);
    anMD(:,i) = mean(tnMD(:,st(i):en(i)),2);
    anLP(:,i) = mean(tnLP(:,st(i):en(i)),2);
    anLD(:,i) = mean(tnLD(:,st(i):en(i)),2);
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

%% Diffs - 1yr time series

%% Diffs - 1 yr vert distrib

%%
zl_ts = repmat(zl,1,length(tdays));
[zl2,tts] = meshgrid(tdays,-1*zl);

%flipud seafloor - offline
mMD2 = flipud(mMD2);
mLD2 = flipud(mLD2);
mBE2 = flipud(mBE2);
Bts = flipud(tts);


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
nmo = length(P);
nyr = nmo/12;
y = (1:nmo)/12;

% All size classes of all
figure(1)
plot(y,log10(B),'Linewidth',1); hold on;
plot(y,log10(tSF),'Linewidth',1); hold on;
plot(y,log10(tMF),'Linewidth',1); hold on;
plot(y,log10(tSP),'Linewidth',1); hold on;
plot(y,log10(tMP),'Linewidth',1); hold on;
plot(y,log10(tLP),'Linewidth',1); hold on;
plot(y,log10(tSD),'Linewidth',1); hold on;
plot(y,log10(tMD),'Linewidth',1); hold on;
plot(y,log10(tLD),'Linewidth',1); hold on;
legend('B','SF','MF','SP','MP','LP','SD','MD','LD')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-5 2])
xlabel('Time (y)')
%ylabel('log10 Biomass (g m^-^2)')
ylabel('log10 integrated Biomass (g m^-^2)')
title('Spinup')
stamp('')
print('-dpng',[ppath exper '_ts_all_sizes.png'])

% Fn Types
figure(2)
plot(y,log10(B),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,log10(F),'r','Linewidth',2); hold on;
plot(y,log10(P),'b','Linewidth',2); hold on;
plot(y,log10(D),'k','Linewidth',2); hold on;
legend('B','F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-5 2])
xlabel('Time (y)')
%ylabel('log10 Biomass (g m^-^2)')
ylabel('log10 integrated Biomass (g m^-^2)')
title('Spinup')
stamp('')
print('-dpng',[ppath exper '_ts_all_types.png'])

%% Plots
figure(3)
subplot(3,3,1)
pcolor(tts,zl2,log10(S_Sml_f));
shading flat;
colorbar
clim([-5 2])
title('SF')

subplot(3,3,2)
pcolor(tts,zl2,log10(S_Sml_p));
shading flat;
colorbar
clim([-5 2])
title('SP')

subplot(3,3,3)
pcolor(tts,zl2,log10(S_Sml_d));
shading flat;
colorbar
clim([-5 2])
title('SD')

subplot(3,3,4)
pcolor(tts,zl2,log10(S_Med_f));
shading flat;
colorbar
title('MF')
clim([-5 2])
ylabel('Depth (m)')

subplot(3,3,5)
pcolor(tts,zl2,log10(S_Med_p));
shading flat;
colorbar
clim([-5 2])
title('MP')

subplot(3,3,6)
pcolor(tts,zl2,log10(S_Med_d(2:76,:)));
shading flat;
colorbar
clim([-5 2])
title('MD')

subplot(3,3,7)
pcolor(tts,zl2,log10(S_Bent_bio(2:76,:)));
shading flat;
colorbar
clim([-5 2])
title('BE')

subplot(3,3,8)
pcolor(tts,zl2,log10(S_Lrg_p));
shading flat;
colorbar
title('LP')
clim([-5 2])
xlabel('Time (y)')

subplot(3,3,9)
pcolor(tts,zl2,log10(S_Lrg_d(2:76,:)));
shading flat;
colorbar
clim([-5 2])
title('LD')
print('-dpng',[ppath exper '_depth_ts_fntypes.png'])


%% All zoom

figure(6)
subplot(3,3,1)
pcolor(tts(1:45,:),zl2(1:45,:),log10(S_Sml_f(1:45,:)));
shading flat;
colorbar
clim([-5 0])
title('SF')

subplot(3,3,2)
pcolor(tts(1:45,:),zl2(1:45,:),log10(S_Sml_p(1:45,:)));
shading flat;
colorbar
clim([-5 0])
title('SP')

subplot(3,3,3)
pcolor(tts(1:45,:),zl2(1:45,:),log10(S_Sml_d(1:45,:)));
shading flat;
colorbar
clim([-5 0])
title('SD')

subplot(3,3,4)
pcolor(tts(1:45,:),zl2(1:45,:),log10(S_Med_f(1:45,:)));
shading flat;
colorbar
title('MF')
clim([-5 0])
ylabel('Depth (m)')

subplot(3,3,5)
pcolor(tts(1:45,:),zl2(1:45,:),log10(S_Med_p(1:45,:)));
shading flat;
colorbar
clim([-5 0])
title('MP')

subplot(3,3,6)
pcolor(tts(73:75,:),zl2(73:75,:),log10(S_Med_d(74:76,:)));
shading flat;
colorbar
clim([-5 0])
title('MD')

subplot(3,3,7)
pcolor(tts(73:75,:),zl2(73:75,:),log10(S_Bent_bio(74:76,:)));
shading flat;
colorbar
clim([-3 1])
title('BE')

subplot(3,3,8)
pcolor(tts(1:45,:),zl2(1:45,:),log10(S_Lrg_p(1:45,:)));
shading flat;
colorbar
title('LP')
clim([-3 1])
xlabel('Time (y)')

subplot(3,3,9)
pcolor(tts(73:75,:),zl2(73:75,:),log10(S_Lrg_d(74:76,:)));
shading flat;
colorbar
clim([-3 1])
title('LD')
print('-dpng',[ppath exper '_depth_ts_fntypes_zoom.png'])


