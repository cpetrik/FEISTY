% Visualize output of POEM with different initial biomasses
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
locs = abbrev';

%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
close all

% colors
load('MyColormaps.mat')
cm9=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...  %b
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

set(groot,'defaultAxesColorOrder',cm9);

%% Loop over IC biomasses
icb=[1e-10;1;10;1e3];
Fts = NaN*ones(length(icb), 12*200);
Pts = Fts;
Dts = Fts;
Fts_locs = NaN*ones(length(locs),12*150,length(icb));
Pts_locs = Fts_locs;
Dts_locs = Fts_locs;
Bts_locs = Fts_locs;
for n=1:length(icb)
    
    load([fpath 'Means_Climatol_' harv '_ICbiom',num2str(icb(n)),'_' cfile '.mat']);
    
    %% Plots in time
    y = time;
    nt = length(time);
    F = sf_tmean(1:nt)+mf_tmean;
    P = sp_tmean+mp_tmean+lp_tmean;
    D = sd_tmean+md_tmean+ld_tmean;
    
    if (n<4)
        Fts(n,1:1800) = F;
        Pts(n,1:1800) = P;
        Dts(n,1:1800) = D;
    else
        Fts(n,:) = F;
        Pts(n,:) = P;
        Dts(n,:) = D;
    end
    
    Flocs = sf_ts_locs+mf_ts_locs;
    Plocs = sp_ts_locs+mp_ts_locs+lp_ts_locs;
    Dlocs = sd_ts_locs+md_ts_locs+ld_ts_locs;
    
    Fts_locs(:,:,n) = Flocs(:,1:1800);
    Pts_locs(:,:,n) = Plocs(:,1:1800);
    Dts_locs(:,:,n) = Dlocs(:,1:1800);
    Bts_locs(:,:,n) = b_ts_locs(:,1:1800);
    
    %% All size classes of all
    
    f4=figure(4);
    subplot(2,2,n)
    plot(y,log10(sf_tmean(1:nt)),'Linewidth',1); hold on;
    plot(y,log10(mf_tmean),'Linewidth',1); hold on;
    plot(y,log10(sp_tmean),'Linewidth',1); hold on;
    plot(y,log10(mp_tmean),'Linewidth',1); hold on;
    plot(y,log10(lp_tmean),'Linewidth',1); hold on;
    plot(y,log10(sd_tmean),'Linewidth',1); hold on;
    plot(y,log10(md_tmean),'Linewidth',1); hold on;
    plot(y,log10(ld_tmean),'Linewidth',1); hold on;
    legend('SF','MF','SP','MP','LP','SD','MD','LD')
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    ylim([-2.5 1.5])
    xlabel('Time (mo)')
    ylabel('log10 Biomass (g m^-^2)')
    title(['ICbiom = ' num2str(icb(n))])
    
    f5=figure(5);
    subplot(2,2,n)
    plot(y,log10(F),'r','Linewidth',2); hold on;
    plot(y,log10(P),'b','Linewidth',2); hold on;
    plot(y,log10(D),'k','Linewidth',2); hold on;
    legend('F','P','D')
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    ylim([-2.5 1.5])
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    title(['ICbiom = ' num2str(icb(n))])
    
end
print(f4,'-dpng',[ppath 'Climatol_' harv '_all_sizes_ICbiom.png'])
print(f5,'-dpng',[ppath 'Climatol_' harv '_all_types_ICbiom.png'])

%%
cm4=[0 0.7 0;...   %g
    1 0 0;...     %r
    0 0 0.75;...  %b
    0 0 0];...      %black
    set(groot,'defaultAxesColorOrder',cm4);

figure(1)
plot(y,log10(Fts),'Linewidth',2); hold on;
plot(y,repmat(log10(Fts(3,1800)),length(y),1),'--b','Linewidth',2); hold on;
legend(num2str(icb))
xlim([y(1) y(end)])
ylim([-1 1])
xlabel('Time (y)')
ylabel('log10 Biomass (g m^-^2)')
title('F')
print('-dpng',[ppath 'Climatol_' harv '_F_ICbiom.png'])

figure(2)
plot(y,log10(Pts),'Linewidth',2); hold on;
plot(y,repmat(log10(Pts(3,1800)),length(y),1),'--b','Linewidth',2); hold on;
legend(num2str(icb))
xlim([y(1) y(end)])
ylim([-1 1])
xlabel('Time (y)')
ylabel('log10 Biomass (g m^-^2)')
title('P')
print('-dpng',[ppath 'Climatol_' harv '_P_ICbiom.png'])
%%
figure(3)
plot(y,log10(Dts),'Linewidth',2); hold on;
legend(num2str(icb))
xlim([y(1) y(end)])
ylim([-1 1])
xlabel('Time (y)')
ylabel('log10 Biomass (g m^-^2)')
title('D')
print('-dpng',[ppath 'Climatol_' harv '_D_ICbiom.png'])

%% Loop over locations
% for i=1:length(locs)
%     
%     f6=figure(6);
%     subplot(4,4,i)
%     plot(y,log10(squeeze(Fts_locs(i,:,:)))); hold on;
%     %legend(num2str(icb))
%     xlim([y(1) y(end)])
%     ylim([-2 2])
%     if (i>12)
%         xlabel('Time (y)')
%     end
%     if (i==5)
%         ylabel('log10 Biomass (g m^-^2)')
%     end
%     title([locs{i} ' F'])
%     
%     f7=figure(7);
%     subplot(4,4,i)
%     plot(y,log10(squeeze(Pts_locs(i,:,:)))); hold on;
%     %legend(num2str(icb))
%     xlim([y(1) y(end)])
%     ylim([-2 2])
%     if (i>12)
%         xlabel('Time (y)')
%     end
%     if (i==5)
%         ylabel('log10 Biomass (g m^-^2)')
%     end
%     title([locs{i} ' P'])
%     
%     f8=figure(8);
%     subplot(4,4,i)
%     plot(y,log10(squeeze(Dts_locs(i,:,:)))); hold on;
%     %legend(num2str(icb))
%     xlim([y(1) y(end)])
%     ylim([-2 2])
%     if (i>12)
%         xlabel('Time (y)')
%     end
%     if (i==5)
%         ylabel('log10 Biomass (g m^-^2)')
%     end
%     title([locs{i} ' D'])
%     
% end
% print(f6,'-dpng',[ppath 'Climatol_' harv '_Flocs_ICbiom.png'])
% print(f7,'-dpng',[ppath 'Climatol_' harv '_Plocs_ICbiom.png'])
% print(f8,'-dpng',[ppath 'Climatol_' harv '_Dlocs_ICbiom.png'])
    



