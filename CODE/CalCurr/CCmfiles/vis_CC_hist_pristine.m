% Visualize output of POEM
% Historic time period (1861-2005) at all locations
% 145 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];

load([fpath 'Means_Historic_pristine_' cfile '.mat']);

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

%% Pick which time period mean
% 1990-1994
sp_smean=sp_mean5;
sf_smean=sf_mean5;
sd_smean=sd_mean5;
mp_smean=mp_mean5;
mf_smean=mf_mean5;
md_smean=md_mean5;
lp_smean=lp_mean5;
ld_smean=ld_mean5;
b_smean=b_mean5;

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
ylim([-5 2])
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
title('Historic Pristine')
stamp(cfile)
print('-dpng',[ppath 'Hist_pristine_all_sizes.png'])

figure(2)
plot(y,log10(F),'r','Linewidth',2); hold on;
plot(y,log10(P),'b','Linewidth',2); hold on;
plot(y,log10(D),'k','Linewidth',2); hold on;
legend('F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-5 2])
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
title(['Historic Pristine'])
print('-dpng',[ppath 'Hist_pristine_all_types.png'])

% %% Recruitment
% st=1:12:length(y);
% en=12:12:length(y);
% MPy = NaN*ones((length(time)/12),1);
% SFy = MPy;
% MDy = MPy;
% 
% for n=1:(length(time)/12)
%     MPy(n) = nansum(mp_trec(st(n):en(n)));
%     SFy(n) = nansum(sf_trec(st(n):en(n)));
%     MDy(n) = nansum(md_trec(st(n):en(n)));
% end
% 
% %%
% figure(3)
% subplot(3,1,2)
% plot(1861:2005,log10(MPy),'b','Linewidth',2); hold on;
% xlim([1861 2005])
% ylabel('log10 annual Recruitment (g m^-^2)')
% title('Large pelagics')
% stamp(cfile)
% 
% subplot(3,1,1)
% plot(1861:2005,log10(SFy),'r','Linewidth',2); hold on;
% xlim([1861 2005])
% title('Forage fishes')
% 
% subplot(3,1,3)
% plot(1861:2005,log10(MDy),'k','Linewidth',2); hold on;
% xlim([1861 2005])
% title('Demersals')
% print('-dpng',[ppath 'Hist_pristine_recruitment.png'])

%% Time series
all_bio = sp_tmean+sf_tmean+sd_tmean+mp_tmean+mf_tmean+md_tmean+lp_tmean+ld_tmean;

figure(70)
plot(y,log10(all_bio),'k','LineWidth',2)
ylim([-2 2])
xlim([1860 2005])
xlabel('Year')
ylabel('All fish mean biomass (g/m^2)')
title('Historic pristine')
print('-dpng',[ppath 'Hist_pristine_ts_mbio.png'])

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
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%% bent
figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2.5 0.5]);
hcb = colorbar('h');
ylim(hcb,[-2.5 0.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic pristine 1956-2005 log10 mean benthic biomass (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_pristine_global_BENT.png'])

%
% mgZb = (Zb/9)*1e3;
% figure(5)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(mgZb))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-0.8 2.3]);
% hcb = colorbar('h');
% ylim(hcb,[-0.8 2.3])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Historic pristine 1956-2005 log10 mean benthic biomass (mg C m^-^2)')
% stamp(cfile)
% print('-dpng',[ppath 'Hist_pristine_global_BENT_mgC.png'])

% % sp
% figure(6)
% surf(geolon_t,geolat_t,log10(Zsp)); view(2); hold on;
% shading flat
% title('Historic pristine 1956-2005 log10 mean Larval P biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Hist_pristine_global_SP.png'])
%
% % sf
% figure(7)
% surf(geolon_t,geolat_t,log10(Zsf)); view(2); hold on;
% shading flat
% title('Historic pristine 1956-2005 log10 mean Larval F biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Hist_pristine_global_SF.png'])
%
% % sd
% figure(8)
% surf(geolon_t,geolat_t,log10(Zsd)); view(2); hold on;
% shading flat
% title('Historic pristine 1956-2005 log10 mean Larval D biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Hist_pristine_global_SD.png'])
%
% % mp
% figure(9)
% surf(geolon_t,geolat_t,log10(Zmp)); view(2); hold on;
% shading flat
% title('Historic pristine 1956-2005 log10 mean Juvenile P biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Hist_pristine_global_MP.png'])
%
% % mf
% figure(10)
% surf(geolon_t,geolat_t,log10(Zmf)); view(2); hold on;
% shading flat
% title('Historic pristine 1956-2005 log10 mean Adult F biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Hist_pristine_global_MF.png'])
%
% % md
% figure(11)
% surf(geolon_t,geolat_t,log10(Zmd)); view(2); hold on;
% shading flat
% title('Historic pristine 1956-2005 log10 mean Juvenile D biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Hist_pristine_global_MD.png'])
%
% % lp
% figure(12)
% surf(geolon_t,geolat_t,log10(Zlp)); view(2); hold on;
% shading flat
% title('Historic pristine 1956-2005 log10 mean Adult P biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Hist_pristine_global_LP.png'])
%
% % ld
% figure(13)
% surf(geolon_t,geolat_t,log10(Zld)); view(2); hold on;
% shading flat
% title('Historic pristine 1956-2005 log10 mean Adult D biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Hist_pristine_global_LD.png'])

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
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 2]);
hcb = colorbar('h');
ylim(hcb,[-1 2])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic pristine 1956-2005 log10 mean biomass All Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_pristine_global_All.png'])

% all F
figure(15)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic pristine 1956-2005 log10 mean biomass All F (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_pristine_global_AllF.png'])

% all D
figure(16)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic pristine 1956-2005 log10 mean biomass All D (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_pristine_global_AllD.png'])

% All P
figure(17)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic pristine 1956-2005 log10 mean biomass All P (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_pristine_global_AllP.png'])

%% All 4 on subplots
figure(18)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 mean All F (g m^-^2)')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('log10 mean All D (g m^-^2)')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('log10 mean All P (g m^-^2)')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('log10 mean All fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_pristine_global_All_subplot.png'])

%% Ratios on subplots red-white-blue
% 3 figure subplot P:D, P:F, M:L
figure(19)
subplot('Position',[0 0.53 0.5 0.5])
%P:D
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Fraction Large Pelagics vs. Demersals')

%P:F
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracPF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Fraction Large Pelagics vs. Forage Fishes')

%L:M
subplot('Position',[0.25 0.0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracLM)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Fraction Large vs. Medium')
stamp(cfile)
print('-dpng',[ppath 'Hist_pristine_global_ratios_subplot.png'])

