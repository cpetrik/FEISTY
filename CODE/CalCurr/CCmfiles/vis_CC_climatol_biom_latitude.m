% Visualize output of POEM
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files
% Plot biomass vs. latitude

clear all
close all
warning off

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end

load([fpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

% plot info
[ni,nj]=size(lon);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

geolat_t=lat;
geolon_t=lon;

cmS=cbrewer('div','Spectral',11);
cmS2=cmS(6:end,:);
cmS2=flipud(cmS2);

cmRB=cbrewer('div','RdYlBu',40);
cmRB=flipud(cmRB);
cmPO=cbrewer('div','PuOr',40);
cmPO=flipud(cmPO);

% save('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat',...
%     'mycmap','-append')
load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
cmPR = mycmap;

%% All fish
lat_vec = lat(ID);
dep_vec = depth(ID);
shelf = (dep_vec<200);

AllF = sf_mean+mf_mean;
AllP = sp_mean+mp_mean+lp_mean;
AllD = sd_mean+md_mean+ld_mean;
All = AllF+AllP+AllD;
AllS = sp_mean+sf_mean+sd_mean;
AllM = mp_mean+mf_mean+md_mean;
AllL = lp_mean+ld_mean;
FracPD = AllP ./ (AllP+AllD);
FracPF = AllP ./ (AllP+AllF);
FracLM = AllL ./ (AllM+AllL);

%% All Biom

[N,edges,bin] = histcounts(lat_vec,35);

lat_bin = nan*ones(length(lat_vec),1);
for n = 1:36
    lat_bin(bin==n) = edges(n);
end

mlatF = nan*ones(length(edges),1);
mlatP = mlatF;
mlatD = mlatF;
tlatF = mlatF;
tlatP = mlatF;
tlatD = mlatF;
mlatS = mlatF;
mlatM = mlatF;
mlatL = mlatF;
tlatS = mlatF;
tlatM = mlatF;
tlatL = mlatF;
for n = 1:35
    mlatF(n) = nanmean(AllF(bin==n));
    tlatF(n) = nansum(AllF(bin==n));
    mlatP(n) = nanmean(AllP(bin==n));
    tlatP(n) = nansum(AllP(bin==n));
    mlatD(n) = nanmean(AllD(bin==n));
    tlatD(n) = nansum(AllD(bin==n));
    
    mlatS(n) = nanmean(AllS(bin==n));
    tlatS(n) = nansum(AllS(bin==n));
    mlatM(n) = nanmean(AllM(bin==n));
    tlatM(n) = nansum(AllM(bin==n));
    mlatL(n) = nanmean(AllL(bin==n));
    tlatL(n) = nansum(AllL(bin==n));
end

%% Boxplot
figure(1)
subplot(2,2,1)
boxplot(log10(AllF),lat_bin,'Symbol','k.','OutlierSize',1,'PlotStyle','compact')
ylim([-4 2])
title('F')
subplot(2,2,2)
boxplot(log10(AllP),lat_bin,'Symbol','k.','OutlierSize',1,'PlotStyle','compact')
ylim([-4 2])
title('P')
subplot(2,2,3)
boxplot(log10(AllD),lat_bin,'Symbol','k.','OutlierSize',1,'PlotStyle','compact')
ylim([-3 2])
title('D')
subplot(2,2,4)
boxplot(log10(All),lat_bin,'Symbol','k.','OutlierSize',1,'PlotStyle','compact')
ylim([-3 2])
title('All')
print('-dpng',[ppath 'Climatol_' harv '_biom_latitude_boxplot.png'])

%% Boxplot fractions
figure(2)
subplot(2,2,1)
boxplot(FracPD,lat_bin,'Symbol','k.','OutlierSize',1,'PlotStyle','compact')
ylim([0 1])
title('P / (P+D)')
set(gca,'XTick',edges(1):2:edges(end),'XTickLabel',[])
subplot(2,2,2)
boxplot(FracPF,lat_bin,'Symbol','k.','OutlierSize',1,'PlotStyle','compact')
ylim([0 1])
title('P / (P+F)')
subplot(2,2,3)
boxplot(FracLM,lat_bin,'Symbol','k.','OutlierSize',1,'PlotStyle','compact')
ylim([0 1])
title('L / (L+M)')
print('-dpng',[ppath 'Climatol_' harv '_fraction_latitude_boxplot.png'])

%% Plot means
figure(3)
subplot(2,2,1)
plot(log10(mlatF),-edges,'k')
ylim([-90 90])
xlim([-4 2])
ylabel('Latitude')
title('mean F')
subplot(2,2,2)
plot(log10(mlatP),-edges,'k')
ylim([-90 90])
xlim([-4 2])
title('mean P')
subplot(2,2,3)
plot(log10(mlatD),-edges,'k')
ylim([-90 90])
xlim([-4 2])
ylabel('Latitude')
title('mean D')
subplot(2,2,4)
plot(log10(mlatD+mlatP+mlatF),-edges,'k')
ylim([-90 90])
xlim([-4 2])
title('mean All')


%% Plot together
figure(4)
plot(log10(mlatD+mlatP+mlatF),-edges,'k','LineWidth',2); hold on
plot(log10(mlatF),-edges,'color',[0.98 0.19 0],'LineWidth',2); hold on
plot(log10(mlatP),-edges,'color',[0 0 0.65],'LineWidth',2); hold on
plot(log10(mlatD),-edges,'color',[0.1 0.65 0.10],'LineWidth',2); hold on
ylabel('Latitude')
xlabel('log_1_0 Mean biomass (g m^-^2)')
ylim([-90 90])
xlim([-4 2])
legend('Total','F','P','D')
legend('location','west')
print('-dpng',[ppath 'Climatol_' harv '_mean_biom_latitude_plot.png'])


%% Save for using R instead
bTab = table(edges',mlatF,mlatP,mlatD);

coast = int32(shelf);

lat_bin = lat_bin-1;
Tab = table(lat_vec,lat_bin,coast,AllF,AllP,AllD,All,AllS,AllM,AllL,...
    FracPD,FracPF,FracLM);

writetable(Tab,[fpath 'Lat_bio_Climatol_' harv '_' cfile '.csv'],...
    'Delimiter',',')
save([fpath 'Lat_bio_Climatol_' harv '_' cfile '.mat'],'N','edges','bin',...
    'lat_vec','lat_bin','AllF','AllP','AllD','All','AllS','AllM','AllL',...
    'FracPD','FracPF','FracLM');





