% Visualize output of FEISTY forced with CESM
% NPP control exper 1850-2100 initialized with spinup biomass
% Time series plots and maps

clear all
close all

%% Fish data
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP_CESM/';
cfile = 'NoDc_enc70-b200_m4-b175-k086_c20-b250_noD_J100_A050_Sm025_nmort1_BE00_noCC_RE00100';
harv = 'pristine';
tharv = 'F=0';
fpath=['/Volumes/FEISTY/NC/FishMIP/CESM1-BEC/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'Means_NPP_cont_' cfile '.mat']);

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
%y = 1850 + (t)/12;


%% All size classes of all
figure(5)
plot(y,log10(sf_tmean),'Linewidth',1); hold on;
plot(y,log10(mf_tmean),'Linewidth',1); hold on;
plot(y,log10(sp_tmean),'Linewidth',1); hold on;
plot(y,log10(mp_tmean),'Linewidth',1); hold on;
plot(y,log10(lp_tmean),'Linewidth',1); hold on;
legend('SF','MF','SP','MP','LP')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-2 2])
xlabel('Year')
ylabel('log_1_0 Biomass (g m^-^2)')
title(['NPP control ' harv])
stamp(cfile)
print('-dpng',[ppath 'NPP_cont_',harv,'_ts_all_sizes.png'])

%%
F = sf_tmean+mf_tmean;
P = sp_tmean+mp_tmean+lp_tmean;

figure(7)
plot(y,log10(F),'r','Linewidth',2); hold on;
plot(y,log10(P),'b','Linewidth',2); hold on;
legend('F','P')
legend('location','northeast')
xlim([y(1) y(end)])
ylim([-1 1])
xlabel('Year')
ylabel('log_1_0 Biomass (g m^-^2)')
title(['NPP control ' harv])
stamp(cfile)
print('-dpng',[ppath 'NPP_cont_',harv,'_ts_all_types.png'])

%% Save for plotting together
NF = F;
NP = P;
Ny = y;
save([fpath 'Ts_Means_NPP_cont_' cfile '.mat'],'NF','NP','Ny');


 
