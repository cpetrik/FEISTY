% Visualize output of POEM
% Spinup with IPSL downscaled model
% 200 years
% Saved as mat files

clear 
close all

%% Fish data
pp = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CC/';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';
vers = 'HAD';
mod = 'hadley';
harv = 'All_fishobs';

fpath=['/Volumes/petrik-lab/Feisty/NC/NEMURO/',cfile,'/',vers,'/'];

ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end
load([fpath 'Means_Spinup_' vers '_' harv '_' cfile '.mat']);

% Map data
cpath = ['/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/' vers 'down/'];
load([cpath 'feisty_' mod '_gridspec.mat'])%,'LON','LAT');
load([cpath 'Data_grid_nemuro_' mod '.mat']);

%%
F = sf_tmean+mf_tmean;
P = sp_tmean+mp_tmean+lp_tmean;
D = sd_tmean+md_tmean+ld_tmean;
B = b_tmean;

%%
[ni,nj]=size(LON);
plotminlat=30; %Set these bounds for your data
plotmaxlat=48;
plotminlon=-134;
plotmaxlon=-115;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

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
    
cm21=[1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    0 1 1;...     %c
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.75 0.75 0.75;... %lt grey
    0.5 0.5 0.5;...    %med grey
    49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    1 1 0;...      %yellow
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255]; %peach

set(groot,'defaultAxesColorOrder',cm10);

cmBP50=cbrewer('seq','BuPu',50,'PCHIP');

%% Plots in time
y = 1:length(sp_tmean); %time;

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
ylim([-5 2])
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
title(['Spinup ' harv])
stamp(cfile)
print('-dpng',[ppath 'Spinup_' vers '_' harv '_all_sizes.png'])

%%
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
ylabel('log10 Biomass (g m^-^2)')
title(['Spinup ' harv])
print('-dpng',[ppath 'Spinup_' vers '_' harv '_all_types.png'])


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

Zsf(GRD.ID)=sf_mean;
Zsp(GRD.ID)=sp_mean;
Zsd(GRD.ID)=sd_mean;
Zmf(GRD.ID)=mf_mean;
Zmp(GRD.ID)=mp_mean;
Zmd(GRD.ID)=md_mean;
Zlp(GRD.ID)=lp_mean;
Zld(GRD.ID)=ld_mean;
Zb(GRD.ID)=b_mean;

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

%% bent
figure(3)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','on','FLineWidth',1)
surfm(LAT,LON,log10(Zb))
colormap(cmBP50)
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean benthic biomass (g m^-^2)')
stamp(vers)
print('-dpng',[ppath 'Spinup_' vers '_' harv '_BENT.png'])

%% All 4 on subplots
figure(4)
% all F
subplot('Position',[0 0.5 0.5 0.5])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','on','FLineWidth',1)
surfm(LAT,LON,log10(AllF))
colormap(cmBP50)
clim([-1 3]);
colorbar('Position',[0.475 0.25 0.035 0.5],'orientation','vertical')
set(gcf,'renderer','painters')
text(0,0.93,'\bf log_1_0 mean All F (g m^-^2)','HorizontalAlignment','center')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','on','FLineWidth',1)
surfm(LAT,LON,log10(AllD))
colormap(cmBP50)
clim([-1 3]);
set(gcf,'renderer','painters')
text(0,0.93,'\bf log_1_0 mean All D (g m^-^2)','HorizontalAlignment','center')
%title('log10 mean All D (g m^-^2)')

% All P
subplot('Position',[0.5 0.5 0.5 0.5])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','on','FLineWidth',1)
surfm(LAT,LON,log10(AllP))
colormap(cmBP50)
clim([-1 3]);
set(gcf,'renderer','painters')
text(0,0.93,'\bf log_1_0 mean All P (g m^-^2)','HorizontalAlignment','center')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','on','FLineWidth',1)
surfm(LAT,LON,log10(All))
colormap(cmBP50)
clim([-1 3]);
set(gcf,'renderer','painters')
text(0,0.93,'\bf log_1_0 mean All fishes (g m^-^2)','HorizontalAlignment','center')
stamp(vers)
print('-dpng',[ppath 'Spinup_' vers '_' harv '_All_subplot.png'])

%% Ratios on subplots red-white-blue
% 3 figure subplot P:D, P:F, M:L
figure(5)
subplot('Position',[0 0.5 0.5 0.5])
%P:D
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','on','FLineWidth',1)
surfm(LAT,LON,FracPD)
cmocean('balance')
clim([0 1]);
colorbar('Position',[0.475 0.55 0.035 0.4],'orientation','vertical')
set(gcf,'renderer','painters')
text(0,0.93,'\bf Fraction Large Pelagics vs. Demersals','HorizontalAlignment','center')

%P:F
subplot('Position',[0.5 0.5 0.5 0.5])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','on','FLineWidth',1)
surfm(LAT,LON,FracPF)
cmocean('balance')
clim([0 1]);
set(gcf,'renderer','painters')
text(0,0.93,'\bf Fraction Large Pelagics vs. Forage Fishes','HorizontalAlignment','center')

%L:M
subplot('Position',[0.25 0.0 0.5 0.5])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','on','FLineWidth',1)
surfm(LAT,LON,FracLM)
cmocean('balance')
clim([0 1]);
set(gcf,'renderer','painters')
text(0,0.93,'\bf Fraction Large vs. Medium','HorizontalAlignment','center')
stamp(vers)
print('-dpng',[ppath 'Spinup_' vers '_' harv '_ratios_subplot.png'])

