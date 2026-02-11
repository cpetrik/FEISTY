% Look at COBALT-FEISTY nonvert 1D at BATS

clear
close all

vpath = '/Volumes/petrik-lab/Feisty/NC/MOM6-1D/BATS_vert/cobalt_only/';

cname = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath = ['/Volumes/petrik-lab/Feisty/NC/MOM6-1D/BATS_vert/offline_feisty/' cname '/'];

exper = 'BATS_spinup_COBALT2004_v4_HPcap_All_fish03';

ppath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/MAPP-METF/NCAR3/online_FEISTY/GFDL_MOM6_1D/Vertical/Offline/';

%%
load([fpath exper '.mat'])

[nz,nt] = size(S_Sml_f);

load([vpath '20040101.ocean_grid_12mo_BATS.mat'],'zl','zl_long_name',...
    'zi')

dz = diff(zi);
dz_mat = repmat(dz,1,nt);

%% Save last month for initializing hindcast runs
Sml_f.bio = S_Sml_f(:,end);
Sml_p.bio = S_Sml_p(:,end);
Sml_d.bio = S_Sml_d(:,end);
Med_f.bio = S_Med_f(:,end);
Med_p.bio = S_Med_p(:,end);
Med_d.bio = S_Med_d(:,end);
Lrg_p.bio = S_Lrg_p(:,end);
Lrg_d.bio = S_Lrg_d(:,end);
BENT.bio  = S_Bent_bio(75,end);

save([fpath 'Last_mo_' exper '_' cname '.mat'],'Sml_f','Sml_p','Sml_d',...
    'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')

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


%% time totals
% tSF=sum(S_Sml_f,1);
% tMF=sum(S_Med_f,1);
% 
% tSP=sum(S_Sml_p,1);
% tMP=sum(S_Med_p,1);
% tLP=sum(S_Lrg_p,1);
% 
% tSD=sum(S_Sml_d,1);
% tMD=sum(S_Med_d,1);
% tLD=sum(S_Lrg_d,1);

%integrated
tSF=sum(S_Sml_f.*dz_mat,1);
tMF=sum(S_Med_f.*dz_mat,1);

tSP=sum(S_Sml_p.*dz_mat,1);
tMP=sum(S_Med_p.*dz_mat,1);
tLP=sum(S_Lrg_p.*dz_mat,1);

tSD=sum(S_Sml_d.*dz_mat,1);
tMD=sum(S_Med_d,1);
tLD=sum(S_Lrg_d,1);

%%
F = tSF + tMF;
P = tSP + tMP + tLP;
D = tSD + tMD + tLD;

B = sum(S_Bent_bio,1);

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

%%
% grid depth and time
%zl_ts = repmat(zl,1,nmo);
[tts,zl2] = meshgrid(y,-1*zl);

S_Med_d(76,:) = zeros;
S_Bent_bio(76,:) = zeros;
S_Lrg_d(76,:) = zeros;

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
xlabel('Time (d)')

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
xlabel('Time (d)')

subplot(3,3,9)
pcolor(tts(73:75,:),zl2(73:75,:),log10(S_Lrg_d(74:76,:)));
shading flat;
colorbar
clim([-3 1])
title('LD')
print('-dpng',[ppath exper '_depth_ts_fntypes_zoom.png'])

