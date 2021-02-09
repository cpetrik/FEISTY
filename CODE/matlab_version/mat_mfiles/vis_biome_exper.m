% Visualize output of Spinup 
% 150 years
% Saved as mat files

clear all
close all

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
exper = 'Biome_exper_control_';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/' exper];

ppath = [pp cfile '/Biome_exper/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

load([fpath 'Means_fished_' cfile '.mat']);

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
y = (time/12);
F = sf_tmean+mf_tmean;
P = sp_tmean+mp_tmean+lp_tmean;
D = sd_tmean+md_tmean+ld_tmean;
B = b_tmean;

% All size classes of all
figure(1)
plot(y,log10(b_tmean),'Linewidth',1); hold on;
plot(y,log10(sf_tmean),'Linewidth',1); hold on;
plot(y,log10(mf_tmean),'Linewidth',1); hold on;
plot(y,log10(sp_tmean),'Linewidth',1); hold on;
plot(y,log10(mp_tmean),'Linewidth',1); hold on;
plot(y,log10(lp_tmean),'Linewidth',1); hold on;
plot(y,log10(sd_tmean),'Linewidth',1); hold on;
plot(y,log10(md_tmean),'Linewidth',1); hold on;
plot(y,log10(ld_tmean),'Linewidth',1); hold on;
legend('B','SF','MF','SP','MP','LP','SD','MD','LD')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-3 1])
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
title('CORE')
stamp(cfile)
%print('-dpng',[ppath 'CORE_fished_all_sizes.png'])

figure(2)
plot(y,log10(B),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,log10(F),'r','Linewidth',2); hold on;
plot(y,log10(P),'b','Linewidth',2); hold on;
plot(y,log10(D),'k','Linewidth',2); hold on;
legend('B','F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-1 1])
xlabel('Time (y)')
ylabel('log10 Biomass (g m^-^2)')
title(['CORE'])
%print('-dpng',[ppath 'CORE_fished_all_types.png'])

%% 
st=1:11:88;
en=11:11:88;

figure(3)
for i=1:8
plot(y,log10(b_mean(st(i):en(i),:))); hold on;
end
xlabel('Time (y)')
ylabel('log10 Biomass (g m^-^2)')
title('Benthos')
legend('SLC','SECCS','SECSS','SCoast','NLC','NECCS','NECSS','NCoast')
legend('location','eastoutside')

figure(4)
for i=1:8
plot(y,log10(mf_mean(st(i):en(i),:))); hold on;
end
xlabel('Time (y)')
ylabel('log10 Biomass (g m^-^2)')
title('MF')
legend('SLC','SECCS','SECSS','SCoast','NLC','NECCS','NECSS','NCoast')
legend('location','eastoutside')

figure(5)
for i=1:8
plot(y,log10(lp_mean(st(i):en(i),:))); hold on;
end
xlabel('Time (y)')
ylabel('log10 Biomass (g m^-^2)')
title('LP')
legend('SLC','SECCS','SECSS','SCoast','NLC','NECCS','NECSS','NCoast')
legend('location','eastoutside')

figure(6)
for i=1:8
plot(y,log10(ld_mean(st(i):en(i),:))); hold on;
end
xlabel('Time (y)')
ylabel('log10 Biomass (g m^-^2)')
title('LD')
legend('SLC','SECCS','SECSS','SCoast','NLC','NECCS','NECSS','NCoast')
legend('location','eastoutside')

%%
plot(y,log10(ld_mean(1:11,:)));

%%
plot(time/12,log10(MF.bio(1:20,:)));