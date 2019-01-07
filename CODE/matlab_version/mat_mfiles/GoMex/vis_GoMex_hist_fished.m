% Visualize output of POEM
% Historic time period (1861-2005) in Gulf of Mexico
% 145 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/GoMex/'];

load([fpath 'Means_Historic_',harv,'_' cfile '.mat']);
load([fpath 'GoMex_ts_means_bio_prod_fish_Historic_' harv '_' cfile '.mat']);

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

%% Pick which time period mean
% 1951-2000
sp_smean=sp_mean50;
sf_smean=sf_mean50;
sd_smean=sd_mean50;
mp_smean=mp_mean50;
mf_smean=mf_mean50;
md_smean=md_mean50;
lp_smean=lp_mean50;
ld_smean=ld_mean50;
b_smean=b_mean50;

%% colors
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

%% Plots in time
y = 1860+(1/12):(1/12):2005;
F = sf_tmean+mf_tmean;
P = sp_tmean+mp_tmean+lp_tmean;
D = sd_tmean+md_tmean+ld_tmean;

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
legend('SF','MF','SP','MP','LP','SD','MD','LD')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-3 1])
xlabel('Year')
ylabel('log_1_0 Biomass (g m^-^2)')
title('Historic fished')
stamp(cfile)
print('-dpng',[ppath 'Hist_',harv,'_GoMex_all_sizes.png'])

figure(2)
plot(y,log10(F),'r','Linewidth',2); hold on;
plot(y,log10(P),'b','Linewidth',2); hold on;
plot(y,log10(D),'k','Linewidth',2); hold on;
legend('F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-0.4 0.6])
xlabel('Year')
ylabel('log_1_0 Biomass (g m^-^2)')
title(['Historic fished'])
print('-dpng',[ppath 'Hist_',harv,'_GoMex_all_types.png'])

%% Time series
all_bio = sp_tmean+sf_tmean+sd_tmean+mp_tmean+mf_tmean+md_tmean+lp_tmean+ld_tmean;

figure(70)
plot(y,log10(all_bio),'k','LineWidth',2)
ylim([0 1])
xlim([1860 2005])
xlabel('Year')
ylabel('All fish mean biomass (log_1_0 g/m^2)')
title('Historic fished')
print('-dpng',[ppath 'Hist_',harv,'_GoMex_ts_mbio.png'])

%% Plots in space
[ni,nj]=size(geolon_t);

Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zb=NaN*ones(ni,nj);

Zsf(grid(:,1))=sf_smean;
Zsp(grid(:,1))=sp_smean;
Zsd(grid(:,1))=sd_smean;
Zmf(grid(:,1))=mf_smean;
Zmp(grid(:,1))=mp_smean;
Zmd(grid(:,1))=md_smean;
Zlp(grid(:,1))=lp_smean;
Zld(grid(:,1))=ld_smean;
Zb(grid(:,1))=b_smean;

% plot info
plotminlat=15; %Set these bounds for your data
plotmaxlat=35;
plotminlon=-100;
plotmaxlon=-75;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%% bent
figure(4)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Historic fished 1951-2000 log_1_0 mean benthic biomass (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_',harv,'_GoMex_BENT.png'])

%% Diff maps of all fish
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

%% ALL
figure(14)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 2]);
hcb = colorbar('h');
ylim(hcb,[-1 2])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic fished 1951-2000 log_1_0 mean biomass All Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_',harv,'_GoMex_All.png'])

% all F
figure(15)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic fished 1951-2000 log_1_0 mean biomass All F (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_',harv,'_GoMex_AllF.png'])

% all D
figure(16)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic fished 1951-2000 log_1_0 mean biomass All D (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_',harv,'_GoMex_AllD.png'])

% All P
figure(17)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic fished 1951-2000 log_1_0 mean biomass All P (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_',harv,'_GoMex_AllP.png'])

%% All 4 on subplots
figure(18)
% all F
subplot('Position',[0 0.5 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('log_1_0 mean All F (g m^-^2)')

% all D
subplot('Position',[0 0 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.45 0.25 0.03 0.5],'orientation','vertical')
set(gcf,'renderer','painters')
title('log_1_0 mean All D (g m^-^2)')

% All P
subplot('Position',[0.5 0.5 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('log_1_0 mean All P (g m^-^2)')

% All
subplot('Position',[0.5 0 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('log_1_0 mean All fishes (g m^-^2)')
%stamp(cfile)
print('-dpng',[ppath 'Hist_',harv,'_GoMex_All_subplot.png'])

%% All 4 on subplots
figure(38)
% all F
subplot('Position',[0 0.5 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('log_1_0 mean All F (g m^-^2)')

% all D
subplot('Position',[0 0 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.45 0.25 0.03 0.5],'orientation','vertical')
set(gcf,'renderer','painters')
title('log_1_0 mean All D (g m^-^2)')

% All P
subplot('Position',[0.5 0.5 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('log_1_0 mean All P (g m^-^2)')

% Bent
subplot('Position',[0.5 0 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('log_1_0 mean B (g m^-^2)')
%stamp(cfile)
print('-dpng',[ppath 'Hist_',harv,'_GoMex_Bent_subplot.png'])

%% Ratios on subplots red-white-blue
% 3 figure subplot P:D, P:F, M:L
figure(19)
subplot('Position',[0 0.5 0.48 0.45])
%P:D
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Fraction Large Pelagics vs. Demersals')

%P:F
subplot('Position',[0.5 0.5 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracPF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Fraction Large Pelagics vs. Forage Fishes')

%L:M
subplot('Position',[0.25 0.0 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracLM)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.7 0.02 0.035 0.425],'orientation','vertical')
set(gcf,'renderer','painters')
title('Fraction Large vs. Medium')
%stamp(cfile)
print('-dpng',[ppath 'Hist_',harv,'_GoMex_ratios_subplot.png'])


%% Save
save([fpath 'Map_vals_Historic_',harv,'_' cfile '.mat'],'b_tmean',...
    'sf_tmean','mf_tmean','sp_tmean','mp_tmean','lp_tmean',...
    'sd_tmean','md_tmean','ld_tmean','F','P','D',...
    'AllF','AllP','AllD','AllS','AllM','AllL',...
    'FracPD','FracPF','FracLM','Zb')



