% Visualize output of FEISTY
% CESM FOSI
% Time series plots and maps

clear all
close all

%% Fish data
%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
mod = 'v15_All_fish03';
%mod = 'All_fish03';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end
load([fpath 'Time_Means_FOSI_' mod '_' cfile '.mat']);
load([fpath 'Space_Means_FOSI_' mod '_' cfile '.mat']);
load([fpath 'Annual_Means_FOSI_' mod '_' cfile '.mat'],'mz_mtf');
% load([fpath 'Time_Means_FOSI_' cfile '.mat']);
% load([fpath 'Space_Means_FOSI_' cfile '.mat']);
% load([fpath 'Annual_Means_FOSI_' cfile '.mat'],'mz_mtf');

% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);

[ni,nj]=size(TLONG);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

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

cmBP50=cbrewer('seq','BuPu',50,'PCHIP');

%% Plots in time
t = 1:length(sp_tmean); %time;
y = t/12;

% All size classes of all
figure(1)
plot(y,log10(sf_tmean),'Linewidth',1); hold on;
plot(y,log10(mf_tmean),'Linewidth',1); hold on;
plot(y,log10(sp_tmean),'Linewidth',1); hold on;
plot(y,log10(mp_tmean),'Linewidth',1); hold on;
plot(y,log10(lp_tmean),'Linewidth',1); hold on;
plot(y,log10(sd_tmean),'Linewidth',1); hold on;
plot(y,log10(md_tmean),'Linewidth',1); hold on;
plot(y,log10(ld_tmean),'Linewidth',1); hold on;
plot(y,log10(b_tmean),'Linewidth',1); hold on;
legend('SF','MF','SP','MP','LP','SD','MD','LD','B')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-3 1])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('FOSI')
stamp(mod)
print('-dpng',[ppath 'FOSI_',mod,'_all_sizes.png'])

%% Types together
F = sf_tmean+mf_tmean;
P = sp_tmean+mp_tmean+lp_tmean;
D = sd_tmean+md_tmean+ld_tmean;
B = b_tmean;

figure(2)
plot(y,log10(B),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,log10(F),'r','Linewidth',2); hold on;
plot(y,log10(P),'b','Linewidth',2); hold on;
plot(y,log10(D),'k','Linewidth',2); hold on;
legend('B','F','P','D')
legend('location','west')
xlim([y(1) y(end)])
%ylim([0.05 0.5])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('FOSI')
stamp(mod)
print('-dpng',[ppath 'FOSI_',mod,'_all_types.png'])

%% Plots in space

Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zb=NaN*ones(ni,nj);

Zsf(GRD.ID)=sf_sbio;
Zsp(GRD.ID)=sp_sbio;
Zsd(GRD.ID)=sd_sbio;
Zmf(GRD.ID)=mf_sbio;
Zmp(GRD.ID)=mp_sbio;
Zmd(GRD.ID)=md_sbio;
Zlp(GRD.ID)=lp_sbio;
Zld(GRD.ID)=ld_sbio;
Zb(GRD.ID)=b_sbio;

All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
AllF = Zsf+Zmf;
AllP = Zsp+Zmp+Zlp;
AllD = Zsd+Zmd+Zld;
AllS = Zsp+Zsf+Zsd;
AllM = Zmp+Zmf+Zmd;
AllL = Zlp+Zld;
FracPD = AllP ./ (AllP+AllD);
FracPF = AllP ./ (AllP+AllF);
FracLM = AllL ./ (AllL+AllM);

%% save outputs for comparison
save([fpath 'Plot_Means_FOSI_' mod '_' cfile '.mat'],'F','P','D','B',...
    'AllF','AllP','AllD','AllS','AllM','AllL');

%% bent
figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(Zb))
colormap(cmBP50)                %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 2]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('FOSI log10 mean benthic biomass (g m^-^2)')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_global_BENT.png'])

%% All 4 on subplots
figure(4)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(AllF))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 mean All F (g m^-^2)')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(AllD))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All D (g m^-^2)')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(AllP))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All P (g m^-^2)')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(All))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All fishes (g m^-^2)')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_global_All_subplot.png'])

%% Ratios on subplots red-white-blue
% 3 figure subplot P:D, P:F, M:L
figure(5)
subplot('Position',[0 0.53 0.5 0.5])
%P:D
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,FracPD)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Fraction Large Pelagics vs. Demersals')

%P:F
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,FracPF)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Fraction Large Pelagics vs. Forage Fishes')

%L:M
subplot('Position',[0.25 0.0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,FracLM)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Fraction Large vs. Medium')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_global_ratios_subplot.png'])

%% MZ loss plots
[nx,nt] = size(mz_mtf);
%Mean fraction
Cmz_smfrac = mz_smfrac;
%Total times it happens over time
Cmz_ttover = mz_ttf/nx;
%Total times it happens in each year in space?
Cmz_stover5 = mz_stf/nt;

%%
%happens whole year
test2=floor(Cmz_stover5);
%histogram(test2)
sum(test2)/length(test2) % = 2.1675

%happens >=50% of year
test4=round(Cmz_stover5);
%histogram(test4)
sum(test4)/length(test4) % = 2.4561

%% Plot in time
figure(6)
plot(y, Cmz_ttover,'k','LineWidth',2); hold on;
xlabel('Years')
ylabel('Fraction of grid cells over-consumed')
print('-dpng',[ppath 'FOSI_',mod '_timeseries_zoop_overcon.png'])

%% Plots in space
CFmz=NaN*ones(ni,nj);
COmz=NaN*ones(ni,nj);
COmz5=NaN*ones(ni,nj);
CFmz(GRD.ID)=Cmz_smfrac;
COmz(GRD.ID)=Cmz_stover5;

% save([bpath 'FOSI_' mod '_ts_map_zoop_overcon.mat'],...
%     'Cmz_ttover','CFmz','COmz','-append');

cmBP=cbrewer('seq','BuPu',10,'PCHIP');
cmBP2 = cmBP;
cmBP2(11,:) = [0 0 0];

%% 2 plots of both maps
figure(7)
%1 - m frac
subplot('Position',[0.01 0.5 0.9 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,CFmz)
colormap(cmBP2)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1.1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'Mean fraction MZ hploss consumed','HorizontalAlignment','center')

%2 - m over
subplot('Position',[0.01 0.01 0.9 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,COmz)
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'Mean times MZ overconsumed','HorizontalAlignment','center')
print('-dpng',[ppath 'FOSI_',mod '_global_zoop_overcon.png'])
