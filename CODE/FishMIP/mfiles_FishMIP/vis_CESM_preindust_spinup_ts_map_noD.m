% Visualize output of FEISTY
% Preindustrial 1800-2100 with spinup biomass
% Time series plots and maps

clear all
close all

% Fish data
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP_CESM/';
cfile = 'NoDc_enc70-b200_m4-b175-k086_c20-b250_noD_J100_A050_Sm025_nmort1_BE00_noCC_RE00100';
fpath=['/Volumes/FEISTY/NC/FishMIP/CESM1-BEC/' cfile '/'];
harv = 'pristine';
tharv = 'F=0';
ppath = [pp cfile];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'Means_last_mo_spinup_' cfile '.mat']);

%% Map data
cpath = '/Volumes/FEISTY/Fish-MIP/CESM/';
load([cpath 'gridspec_cesm.mat']);
load([cpath 'Data_grid_cesm.mat']);
[ni,nj]=size(LON);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
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

%% Plots in time
t = 1:length(sp_tmean); %time;
y = 1850 + (t-1)/12;

% All size classes of all
figure(5)
plot(y,(sf_tmean),'Linewidth',1); hold on;
plot(y,(mf_tmean),'Linewidth',1); hold on;
plot(y,(sp_tmean),'Linewidth',1); hold on;
plot(y,(mp_tmean),'Linewidth',1); hold on;
plot(y,(lp_tmean),'Linewidth',1); hold on;
legend('SF','MF','SP','MP','LP')
legend('location','eastoutside')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title(['Spinup ' harv])
stamp(cfile)
print('-dpng',[ppath '/Spinup_pre_',harv,'_all_sizes.png'])

%%
F = sf_tmean+mf_tmean;
P = sp_tmean+mp_tmean+lp_tmean;

figure(7)
plot(y,(F),'r','Linewidth',2); hold on;
plot(y,(P),'b','Linewidth',2); hold on;
legend('F','P')
legend('location','east')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title(['Spinup ' tharv])
stamp(harv)
print('-dpng',[ppath '/Spinup_pre_',harv,'_all_types.png'])


%% Plots in space

Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);

Zsf(GRD.ID)=Sml_f.bio;
Zsp(GRD.ID)=Sml_p.bio;
Zmf(GRD.ID)=Med_f.bio;
Zmp(GRD.ID)=Med_p.bio;
Zlp(GRD.ID)=Lrg_p.bio;

%% Diff maps of all fish
All = Zsp+Zsf+Zmp+Zmf+Zlp;
AllF = Zsf+Zmf;
AllP = Zsp+Zmp+Zlp;
AllS = Zsp+Zsf;
AllM = Zmp+Zmf;
FracPF = AllP ./ (AllP+AllF);
FracLM = AllP ./ (AllP+AllM);

%% ALL
figure(3)
% 2ndary prod
subplot('Position',[0.01 0.52 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(AllF))
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colormap('jet')
colorbar
set(gcf,'renderer','painters')
text(0,1.75,'All F','HorizontalAlignment','center')
text(-2.75,1.75,'A')

% All
subplot('Position',[0.01 0 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(All))
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colormap('jet')
colorbar
set(gcf,'renderer','painters')
text(0,1.75,'All fish','HorizontalAlignment','center')
text(-2.75,1.75,'C')

% P
subplot('Position',[0.51 0.52 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(AllP))
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colormap('jet')
colorbar
set(gcf,'renderer','painters')
text(0,1.75,'All P','HorizontalAlignment','center')
text(-2.75,1.75,'B')

% D
subplot('Position',[0.51 0 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,FracPF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.75,'Frac P vs F','HorizontalAlignment','center')
text(-2.75,1.75,'D')
%stamp(cfile)
print('-dpng',[ppath '/Spinup_pre_',harv,'_global_types.png'])

