% Visualize output of FEISTY forced with CESM
% Historic 1850-2005 initialized with spinup biomass
% Time series plots and maps

clear all
close all

%% Fish data
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP_CESM/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'pristine';
tharv = 'F=0';
fpath=['/Volumes/GFDL/NC/FishMIP/CESM1-BEC/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'Means_Historic_' harv '_' cfile '.mat']);

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

% % Large Pelagics
% figure(1)
% subplot(4,1,1)
% plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
% plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
% plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
% xlim([y(1) y(end)])
% title('Historic 3km Large Pelagics')
% legend('Larvae','Juveniles','Adults')
% legend('location','southeast')
% stamp(cfile)
% 
% subplot(4,1,2)
% plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
% xlim([y(1) y(end)])
% title('Larvae')
% ylabel('log10 Biomass (g m^-^2)')
% 
% subplot(4,1,3)
% plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
% xlim([y(1) y(end)])
% title('Juveniles')
% 
% subplot(4,1,4)
% plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
% xlim([y(1) y(end)])
% title('Adults')
% xlabel('Time (mo)')
% print('-dpng',[ppath 'Historic_',harv,'_P_time.png'])
% 
% % Forage fishes
% figure(2)
% subplot(3,1,1)
% plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
% plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
% xlim([y(1) y(end)])
% title('Historic 3km Forage Fishes')
% ylabel('log10 Biomass (g m^-^2)')
% legend('Immature','Adults')
% legend('location','southeast')
% stamp(cfile)
% 
% subplot(3,1,2)
% plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
% xlim([y(1) y(end)])
% title('Immature')
% ylabel('log10 Biomass (g m^-^2)')
% 
% subplot(3,1,3)
% plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
% xlim([y(1) y(end)])
% title('Adults')
% xlabel('Time (mo)')
% ylabel('log10 Biomass (g m^-^2)')
% print('-dpng',[ppath 'Historic_',harv,'_F_time.png'])
% 
% % Demersals
% figure(3)
% subplot(4,1,1)
% plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
% plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
% plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
% xlim([y(1) y(end)])
% title('Historic 3km Demersal Piscivores')
% legend('Larvae','Juveniles','Adults')
% legend('location','southeast')
% stamp(cfile)
% 
% subplot(4,1,2)
% plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
% xlim([y(1) y(end)])
% title('Larvae')
% ylabel('log10 Biomass (g m^-^2)')
% 
% subplot(4,1,3)
% plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
% xlim([y(1) y(end)])
% title('Juveniles')
% 
% subplot(4,1,4)
% plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
% xlim([y(1) y(end)])
% title('Adults')
% xlabel('Time (mo)')
% print('-dpng',[ppath 'Historic_',harv,'_D_time.png'])
% 
% % Benthic inverts
% figure(4)
% subplot(2,1,1)
% plot(y,log10(b_tmean),'b','Linewidth',1); hold on;
% xlim([y(1) y(end)])
% title('Historic 3km Benthic Inverts')
% xlabel('Time (mo)')
% ylabel('log10 Biomass (g m^-^2)')
% 
% subplot(2,1,2)
% plot(y,(b_tmean),'b','Linewidth',1); hold on;
% xlim([y(1) y(end)])
% ylabel('Biomass (g m^-^2)')
% print('-dpng',[ppath 'Historic_',harv,'_B_time.png'])

%% All size classes of all
figure(5)
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
ylim([-2 2])
xlabel('Year')
ylabel('log_1_0 Biomass (g m^-^2)')
title(['Historic ' harv])
stamp(cfile)
print('-dpng',[ppath 'Historic_',harv,'_ts_all_sizes.png'])

%%
F = sf_tmean+mf_tmean;
P = sp_tmean+mp_tmean+lp_tmean;
D = sd_tmean+md_tmean+ld_tmean;
B = b_tmean;

figure(7)
plot(y,log10(B),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,log10(F),'r','Linewidth',2); hold on;
plot(y,log10(P),'b','Linewidth',2); hold on;
plot(y,log10(D),'k','Linewidth',2); hold on;
legend('B','F','P','D')
legend('location','northeast')
xlim([y(1) y(end)])
ylim([-1 1])
xlabel('Year')
ylabel('log_1_0 Biomass (g m^-^2)')
title(['Historic ' harv])
stamp(cfile)
print('-dpng',[ppath 'Historic_',harv,'_ts_all_types.png'])

%%
figure(8)
subplot(4,1,1)
plot(y,(F),'r','Linewidth',2); hold on;
xlim([y(1) y(end)])
title(['Historic ' harv ' biomass (g m^-^2)'])
ylabel('Forage fish')

subplot(4,1,3)
plot(y,(D),'k','Linewidth',2); hold on;
xlim([y(1) y(end)])
ylabel('Demersal fish')

subplot(4,1,4)
plot(y,(B),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
xlim([y(1) y(end)])
xlabel('Year')
ylabel('Benthic invertebrates')

subplot(4,1,2)
plot(y,(P),'b','Linewidth',2); hold on;
xlim([y(1) y(end)])
ylabel('Large pelagic fish')
stamp(cfile)
print('-dpng',[ppath 'Historic_',harv,'_all_types_subplot.png'])

%%
cm4=[0 0 1;...    %b
    1 0 0;...     %r
    0.5 0.5 0.5; ...    %med grey
    0 0 0];...      %black

set(groot,'defaultAxesColorOrder',cm4);

figure(90)
yyaxis left
plot(y,(P),'Linewidth',2);
ylabel('P Biomass (g m^-^2)')
xlim([y(1) y(end)])
yyaxis right
plot(y,(F),'Linewidth',2);
ylabel('F Biomass (g m^-^2)')
xlim([y(1) y(end)])
xlabel('Time (y)')
title(['Historic ' harv])
stamp(cfile)
%print('-dpng',[ppath 'Historic_',harv,'_FvP.png'])
 
