% Visualize output of FEISTY
% Preindustrial 1849-1949 to spinup biomass
% Time series plots and maps

clear all
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];

mod = 'fishobs_ctrlclim_15arcmin';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP/Phase3a/';
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end
load([fpath 'Spinup_ctrlclim_All_fishobs_end_cycle_10.mat']);

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

%% means
sf_tmean = mean(ySF,'omitnan');
sp_tmean = mean(ySP,'omitnan');
sd_tmean = mean(ySD,'omitnan');
mf_tmean = mean(yMF,'omitnan');
mp_tmean = mean(yMP,'omitnan');
md_tmean = mean(yMD,'omitnan');
lp_tmean = mean(yLP,'omitnan');
ld_tmean = mean(yLD,'omitnan');
b_tmean = mean(yB,'omitnan');

%% Plots in time
y = 1:10;

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
%ylim([-5 2])
xlabel('# Cycles')
ylabel('log_1_0 Biomass (g m^-^2)')
title('Spinup')
stamp(mod)
print('-dpng',[ppath 'test_cycles_Spinup_empHP_',mod,'_all_sizes.png'])

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
legend('location','east')
xlim([y(1) y(end)])
ylim([-5 2])
xlabel('# Cycles')
ylabel('log_1_0 Biomass (g m^-^2)')
title('Spinup')
stamp(mod)
print('-dpng',[ppath 'test_cycles_Spinup_empHP_',mod,'_all_types.png'])
