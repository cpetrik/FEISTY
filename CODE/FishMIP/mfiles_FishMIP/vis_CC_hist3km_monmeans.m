% Visualize output of POEM
% Historic 1988-2010 with 3km model
% Saved as mat files

clear all
close all

% Fish data
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
%harv = 'pristine';
%tharv = 'F=0';
harv = 'All_fish05';
tharv = 'All fish F=0.5';
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/CalCurr/'];
ppath = [pp cfile '/CC/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'Means_Historic3km_' harv '_' cfile '.mat']);

% Map data
cpath = '/Volumes/GFDL/NEMURO/3km/';
load([cpath 'gridspec_3km.mat'],'LON','LAT');
load([cpath 'Data_grid_3km_hist.mat']);
[ni,nj]=size(LON);
plotminlat=32; %Set these bounds for your data
plotmaxlat=44;
plotminlon=-129;
plotmaxlon=-116;
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
y = 1988 + (t-1)/12;

% All size classes of all
figure(4)
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
xlabel('Time (y)')
ylabel('log10 Biomass (g m^-^2)')
title(['Historic 3km ' tharv])
stamp(harv)
% print('-dpng',[ppath 'Historic3km_',harv,'_log10_all_sizes.png'])

figure(5)
F = sf_tmean+mf_tmean;
P = sp_tmean+mp_tmean+lp_tmean;
D = sd_tmean+md_tmean+ld_tmean;
B = b_tmean;

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
title(['Historic 3km ' tharv])
stamp(harv)
% print('-dpng',[ppath 'Historic3km_',harv,'_log10_all_types.png'])


%% Plots in space
[ni,nj]=size(LON);

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

ocean=NaN*ones(ni,nj);
ocean(GRD.ID)=ones(size(sf_mean));

%% ocean cells
% figure(55)
% surf(LON,LAT,ocean); view(2); hold on;
% shading flat
% title('Water cells')
% colormap('jet')
% colorbar('h')
% caxis([1 2])
% stamp(harv)
% print('-dpng',[ppath 'Ocean_cells.png'])

% bent
% figure(50)
% axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','on','FLineWidth',1)
% surfm(LAT,LON,log10(Zb))
% colormap('jet')
% % load coast;                     %decent looking coastlines
% % h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-2.5 0.5]);
% hcb = colorbar('h');
% ylim(hcb,[-2.5 0.5])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('log10 mean benthic biomass (g m^-^2)')
% stamp(harv)
% print('-dpng',[ppath 'Historic3km_',harv,'_CC_BENT.png'])

%
% mgZb = (Zb/9)*1e3;
% figure(51)
% axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...     'Grid','on','FLineWidth',1)
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(LAT,LON,log10(mgZb))
% colormap('jet')
% caxis([-0.8 2.3]);
% hcb = colorbar('h');
% ylim(hcb,[-0.8 2.3])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('log10 mean benthic biomass (mg C m^-^2)')
% stamp(harv)
% print('-dpng',[ppath 'Historic3km_',harv,'_CC_BENT_mgC.png'])
%
% % sp
% figure(11)
% surf(LON,LAT,log10(Zsp)); view(2); hold on;
% shading flat
% title('log10 mean Larval P biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(harv)
% print('-dpng',[ppath 'Historic3km_',harv,'_CC_SP.png'])
% 
% % sf
% figure(12)
% surf(LON,LAT,log10(Zsf)); view(2); hold on;
% shading flat
% title('log10 mean Larval F biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(harv)
% print('-dpng',[ppath 'Historic3km_',harv,'_CC_SF.png'])
% 
% % sd
% figure(13)
% surf(LON,LAT,log10(Zsd)); view(2); hold on;
% shading flat
% title('log10 mean Larval D biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(harv)
% print('-dpng',[ppath 'Historic3km_',harv,'_CC_SD.png'])
% 
% % mp
% figure(14)
% surf(LON,LAT,log10(Zmp)); view(2); hold on;
% shading flat
% title('log10 mean Juvenile P biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(harv)
% print('-dpng',[ppath 'Historic3km_',harv,'_CC_MP.png'])
% 
% % mf
% figure(15)
% surf(LON,LAT,log10(Zmf)); view(2); hold on;
% shading flat
% title('log10 mean Adult F biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(harv)
% print('-dpng',[ppath 'Historic3km_',harv,'_CC_MF.png'])
% 
% % md
% figure(16)
% surf(LON,LAT,log10(Zmd)); view(2); hold on;
% shading flat
% title('log10 mean Juvenile D biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(harv)
% print('-dpng',[ppath 'Historic3km_',harv,'_CC_MD.png'])
% 
% % lp
% figure(17)
% surf(LON,LAT,log10(Zlp)); view(2); hold on;
% shading flat
% title('log10 mean Adult P biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(harv)
% print('-dpng',[ppath 'Historic3km_',harv,'_CC_LP.png'])
% 
% % ld
% figure(18)
% surf(LON,LAT,log10(Zld)); view(2); hold on;
% shading flat
% title('log10 mean Adult D biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(harv)
% print('-dpng',[ppath 'Historic3km_',harv,'_CC_LD.png'])

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
% figure(21)
% axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','on','FLineWidth',1)
% surfm(LAT,LON,log10(All))
% colormap('jet')
% caxis([-1 2]);
% hcb = colorbar('h');
% ylim(hcb,[-1 2])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('log10 mean biomass All Fishes (g m^-^2)')
% stamp(harv)
% print('-dpng',[ppath 'Historic3km_',harv,'_CC_All.png'])
% 
% % all F
% figure(22)
% axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','on','FLineWidth',1)
% surfm(LAT,LON,log10(AllF))
% colormap('jet')
% caxis([-1 1.5]);
% hcb = colorbar('h');
% set(gcf,'renderer','painters')
% title('log10 mean biomass All F (g m^-^2)')
% stamp(harv)
% print('-dpng',[ppath 'Historic3km_',harv,'_CC_AllF.png'])
% 
% % all D
% figure(23)
% axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','on','FLineWidth',1)
% surfm(LAT,LON,log10(AllD))
% colormap('jet')
% caxis([-1 1.5]);
% hcb = colorbar('h');
% set(gcf,'renderer','painters')
% title('log10 mean biomass All D (g m^-^2)')
% stamp(harv)
% print('-dpng',[ppath 'Historic3km_',harv,'_CC_AllD.png'])
% 
% % All P
% figure(24)
% axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','on','FLineWidth',1)
% surfm(LAT,LON,log10(AllP))
% colormap('jet')
% caxis([-1 1.5]);
% hcb = colorbar('h');
% set(gcf,'renderer','painters')
% title('log10 mean biomass All P (g m^-^2)')
% stamp(harv)
% print('-dpng',[ppath 'Historic3km_',harv,'_CC_AllP.png'])

%% All 4 on subplots
figure(18)
% all F
subplot('Position',[0.025 0.5 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','on','FLineWidth',1)
surfm(LAT,LON,log10(AllF))
colormap('jet')
caxis([-1 2]);
set(gcf,'renderer','painters')
text(0,0.835,'\bf log_1_0 mean All F (g m^-^2)','HorizontalAlignment','center')

% all D
subplot('Position',[0.025 0 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','on','FLineWidth',1)
surfm(LAT,LON,log10(AllD))
colormap('jet')
caxis([-1 2]);
set(gcf,'renderer','painters')
text(0,0.835,'\bf log_1_0 mean All D (g m^-^2)','HorizontalAlignment','center')
%title('log10 mean All D (g m^-^2)')

% All P
subplot('Position',[0.525 0.5 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','on','FLineWidth',1)
surfm(LAT,LON,log10(AllP))
colormap('jet')
caxis([-1 2]);
set(gcf,'renderer','painters')
text(0,0.835,'\bf log_1_0 mean All P (g m^-^2)','HorizontalAlignment','center')

% All
subplot('Position',[0.525 0 0.49 0.49])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','on','FLineWidth',1)
surfm(LAT,LON,log10(All))
colormap('jet')
caxis([-1 2]);
colorbar('Position',[0.5 0.25 0.035 0.5],'orientation','vertical')
set(gcf,'renderer','painters')
text(0,0.835,'\bf log_1_0 mean All fishes (g m^-^2)','HorizontalAlignment','center')
stamp(harv)
print('-dpng',[ppath 'Historic3km_',harv,'_CC_All_subplot.png'])

%% Ratios on subplots red-white-blue
% 3 figure subplot P:D, P:F, M:L
% figure(19)
% subplot('Position',[0 0.5 0.49 0.49])
% %P:D
% axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','on','FLineWidth',1)
% surfm(LAT,LON,FracPD)
% cmocean('balance')
% caxis([0 1]);
% set(gcf,'renderer','painters')
% text(0,0.835,'\bf Fraction Large Pelagics vs. Demersals','HorizontalAlignment','center')
% 
% %P:F
% subplot('Position',[0.5 0.5 0.49 0.49])
% axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','on','FLineWidth',1)
% surfm(LAT,LON,FracPF)
% cmocean('balance')
% caxis([0 1]);
% colorbar('Position',[0.475 0.55 0.035 0.4],'orientation','vertical')
% set(gcf,'renderer','painters')
% text(0,0.835,'\bf Fraction Large Pelagics vs. Forage Fishes','HorizontalAlignment','center')
% 
% %L:M
% subplot('Position',[0.25 0.0 0.49 0.49])
% axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','on','FLineWidth',1)
% surfm(LAT,LON,FracLM)
% cmocean('balance')
% caxis([0 1]);
% set(gcf,'renderer','painters')
% text(0,0.835,'\bf Fraction Large vs. Medium','HorizontalAlignment','center')
% stamp(harv)
% print('-dpng',[ppath 'Historic3km_',harv,'_CC_ratios_subplot.png'])
% 
