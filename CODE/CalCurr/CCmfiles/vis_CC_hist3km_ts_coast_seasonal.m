% Visualize output of POEM
% Historic 1988-2010 with 3km model
% Time series plots and maps

clear all
close all

% Fish data
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
%harv = 'pristine';
%tharv = 'F=0';
harv = 'All_fish01';
tharv = 'All fish F=0.1';
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

% Large Pelagics
% figure(1)
% subplot(4,1,1)
% plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
% plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
% plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
% xlim([y(1) y(end)])
% title('Historic 3km Large Pelagics')
% legend('Larvae','Juveniles','Adults')
% legend('location','southeast')
% stamp(harv)
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
% print('-dpng',[ppath 'Historic3km_',harv,'_P_time.png'])

% Forage fishes
% figure(2)
% subplot(3,1,1)
% plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
% plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
% xlim([y(1) y(end)])
% title('Historic 3km Forage Fishes')
% ylabel('log10 Biomass (g m^-^2)')
% legend('Immature','Adults')
% legend('location','southeast')
% stamp(harv)
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
% print('-dpng',[ppath 'Historic3km_',harv,'_F_time.png'])

% Demersals
% figure(3)
% subplot(4,1,1)
% plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
% plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
% plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
% xlim([y(1) y(end)])
% title('Historic 3km Demersal Piscivores')
% legend('Larvae','Juveniles','Adults')
% legend('location','southeast')
% stamp(harv)
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
% print('-dpng',[ppath 'Historic3km_',harv,'_D_time.png'])

% Benthic inverts
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
% print('-dpng',[ppath 'Historic3km_',harv,'_B_time.png'])

%% All size classes of all
figure(5)
plot(y,(sf_tmean_coast),'Linewidth',1); hold on;
plot(y,(mf_tmean_coast),'Linewidth',1); hold on;
plot(y,(sp_tmean_coast),'Linewidth',1); hold on;
plot(y,(mp_tmean_coast),'Linewidth',1); hold on;
plot(y,(lp_tmean_coast),'Linewidth',1); hold on;
plot(y,(sd_tmean_coast),'Linewidth',1); hold on;
plot(y,(md_tmean_coast),'Linewidth',1); hold on;
plot(y,(ld_tmean_coast),'Linewidth',1); hold on;
plot(y,(b_tmean_coast),'Linewidth',1); hold on;
legend('SF','MF','SP','MP','LP','SD','MD','LD','B')
legend('location','eastoutside')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title(['Historic 3km ' tharv ' coastal (<200m)'])
stamp(harv)
print('-dpng',[ppath 'Historic3km_',harv,'_all_sizes_coast.png'])

figure(6)
subplot(2,1,2)
plot(y,(sf_tmean_coast),'Linewidth',1,'color',cm10(1,:)); hold on;
plot(y,(mf_tmean_coast),'Linewidth',1,'color',cm10(2,:)); hold on;
plot(y,(sp_tmean_coast),'Linewidth',1,'color',cm10(3,:)); hold on;
plot(y,(sd_tmean_coast),'Linewidth',1,'color',cm10(6,:)); hold on;
plot(y,(md_tmean_coast),'Linewidth',1,'color',cm10(7,:)); hold on;
plot(y,(ld_tmean_coast),'Linewidth',1,'color',cm10(8,:)); hold on;
plot(y,(b_tmean_coast),'Linewidth',1,'color',cm10(9,:)); hold on;
legend('SF','MF','SP','SD','MD','LD','B')
legend('location','eastoutside')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')

subplot(2,1,1)
plot(y,(mp_tmean_coast),'r','Linewidth',1); hold on;
plot(y,(lp_tmean_coast),'color',[0.5 0 0],'Linewidth',1); hold on;
xlim([y(1) y(end)])
legend('MP','LP')
legend('location','eastoutside')
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title(['Historic 3km ' tharv ' coastal (<200m)'])
stamp(harv)
print('-dpng',[ppath 'Historic3km_',harv,'_all_sizes_subplot_coast.png'])

F = sf_tmean_coast+mf_tmean_coast;
P = sp_tmean_coast+mp_tmean_coast+lp_tmean_coast;
D = sd_tmean_coast+md_tmean_coast+ld_tmean_coast;
B = b_tmean_coast;

figure(7)
plot(y,(B),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,(F),'r','Linewidth',2); hold on;
plot(y,(P),'b','Linewidth',2); hold on;
plot(y,(D),'k','Linewidth',2); hold on;
legend('B','F','P','D')
legend('location','east')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title(['Historic 3km ' tharv ' coastal (<200m)'])
stamp(harv)
print('-dpng',[ppath 'Historic3km_',harv,'_all_types_coast.png'])

figure(8)
subplot(2,1,2)
plot(y,(B),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,(F),'r','Linewidth',2); hold on;
plot(y,(D),'k','Linewidth',2); hold on;
legend('B','F','D')
legend('location','northwest')
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')

subplot(2,1,1)
plot(y,(P),'b','Linewidth',2); hold on;
legend('P')
legend('location','northwest')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title(['Historic 3km ' tharv ' coastal (<200m)'])
stamp(harv)
print('-dpng',[ppath 'Historic3km_',harv,'_all_types_subplot_coast.png'])

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
title(['Historic 3km ' tharv ' coastal (<200m)'])
stamp(harv)
print('-dpng',[ppath 'Historic3km_',harv,'_FvP_coastal.png'])

%% Seasonal climatology (means)
seas(1,:) = 1:4:size(sB,2); 
seas(2,:) = 2:4:size(sB,2);
seas(3,:) = 3:4:size(sB,2); 
seas(4,:) = 4:4:size(sB,2);
csB = NaN*ones(length(sB),4);
csSF = csB;
csSP = csB;
csSD = csB;
csMF = csB;
csMP = csB;
csMD = csB;
csLP = csB;
csLD = csB;
for i = 1:4
    csB(:,i) = mean(sB(:,seas(i,:)),2);
    csSF(:,i) = mean(sSF(:,seas(i,:)),2);
    csSP(:,i) = mean(sSP(:,seas(i,:)),2);
    csSD(:,i) = mean(sSD(:,seas(i,:)),2);
    csMF(:,i) = mean(sMF(:,seas(i,:)),2);
    csMP(:,i) = mean(sMP(:,seas(i,:)),2);
    csMD(:,i) = mean(sMD(:,seas(i,:)),2);
    csLP(:,i) = mean(sLP(:,seas(i,:)),2);
    csLD(:,i) = mean(sLD(:,seas(i,:)),2);
end
 
%% Plots in space
[ni,nj]=size(LON);

ngrid = NaN*ones(ni,nj); 
nyr = 4;%23;
stext = {'JFM','AMJ','JAS','OND'};

Zsf=NaN*ones(ni*nj,nyr);
Zsp=NaN*ones(ni*nj,nyr);
Zsd=NaN*ones(ni*nj,nyr);
Zmf=NaN*ones(ni*nj,nyr);
Zmp=NaN*ones(ni*nj,nyr);
Zmd=NaN*ones(ni*nj,nyr);
Zlp=NaN*ones(ni*nj,nyr);
Zld=NaN*ones(ni*nj,nyr);
Zb=NaN*ones(ni*nj,nyr);

Zsf(GRD.ID,:)=csSF;
Zsp(GRD.ID,:)=csSP;
Zsd(GRD.ID,:)=csSD;
Zmf(GRD.ID,:)=csMF;
Zmp(GRD.ID,:)=csMP;
Zmd(GRD.ID,:)=csMD;
Zlp(GRD.ID,:)=csLP;
Zld(GRD.ID,:)=csLD;
Zb(GRD.ID,:)=csB;

Zsf = reshape(Zsf,ni,nj,nyr);
Zsp = reshape(Zsp,ni,nj,nyr);
Zsd = reshape(Zsd,ni,nj,nyr);
Zmf = reshape(Zmf,ni,nj,nyr);
Zmp = reshape(Zmp,ni,nj,nyr);
Zmd = reshape(Zmd,ni,nj,nyr);
Zlp = reshape(Zlp,ni,nj,nyr);
Zld = reshape(Zld,ni,nj,nyr);
Zb = reshape(Zb,ni,nj,nyr);

%% bent
figure(9)
for n=1:nyr
    subplot(2,2,n)
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','on','FLineWidth',1)
    surfm(LAT,LON,log10(squeeze(Zb(:,:,n))))
    colormap('jet')
    % load coast;                     %decent looking coastlines
    % h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-2 1]);
    set(gcf,'renderer','painters')
    if (n==1)
        title({'Benthic inverts (log_1_0 g m^-^2)'; stext{n}})
    else
        title(stext{n})
    end
end
colorbar('Position',[0.5 0.35 0.03 0.3],'orientation','vertical')
stamp(harv)
print('-dpng',[ppath 'Historic3km_',harv,'_CC_yrs_BENT_seasonal.png'])

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
figure(10)
for n=1:nyr
    subplot(2,2,n)
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','on','FLineWidth',1)
    surfm(LAT,LON,log10(squeeze(All(:,:,n))))
    colormap('jet')
    % load coast;                     %decent looking coastlines
    % h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([0 2]);
    set(gcf,'renderer','painters')
    if (n==1)
        title({'All fishes (log_1_0 g m^-^2)'; stext{n}})
    else
        title(stext{n})
    end
end
colorbar('Position',[0.5 0.35 0.03 0.3],'orientation','vertical')
stamp(harv)
print('-dpng',[ppath 'Historic3km_',harv,'_CC_yrs_All_seasonal.png'])

% all F
figure(11)
for n=1:nyr
    subplot(2,2,n)
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','on','FLineWidth',1)
    surfm(LAT,LON,log10(squeeze(AllF(:,:,n))))
    colormap('jet')
    % load coast;                     %decent looking coastlines
    % h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-1 1]);
    set(gcf,'renderer','painters')
    if (n==1)
        title({'Forage fish (log_1_0 g m^-^2)'; stext{n}})
    else
        title(stext{n})
    end
end
colorbar('Position',[0.5 0.35 0.03 0.3],'orientation','vertical')
stamp(harv)
print('-dpng',[ppath 'Historic3km_',harv,'_CC_yrs_AllF_seasonal.png'])

% all D
figure(12)
for n=1:nyr
    subplot(2,2,n)
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','on','FLineWidth',1)
    surfm(LAT,LON,log10(squeeze(AllD(:,:,n))))
    colormap('jet')
    % load coast;                     %decent looking coastlines
    % h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-1 1]);
    set(gcf,'renderer','painters')
    if (n==1)
        title({'Demersals (log_1_0 g m^-^2)'; stext{n}})
    else
        title(stext{n})
    end
end
colorbar('Position',[0.5 0.35 0.03 0.3],'orientation','vertical')
stamp(harv)
print('-dpng',[ppath 'Historic3km_',harv,'_CC_yrs_AllD_seasonal.png'])

% All P
figure(13)
for n=1:nyr
    subplot(2,2,n)
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','on','FLineWidth',1)
    surfm(LAT,LON,log10(squeeze(AllP(:,:,n))))
    colormap('jet')
    % load coast;                     %decent looking coastlines
    % h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([0 2]);
    set(gcf,'renderer','painters')
    if (n==1)
        title({'Large pelagics (log_1_0 g m^-^2)'; stext{n}})
    else
        title(stext{n})
    end
end
colorbar('Position',[0.5 0.35 0.03 0.3],'orientation','vertical')
stamp(harv)
print('-dpng',[ppath 'Historic3km_',harv,'_CC_yrs_AllP_seasonal.png'])


