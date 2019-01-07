% Visualize output of POEM
% Historic 1988-2010 with 10km model
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
load([fpath 'Means_Historic10km_' harv '_' cfile '.mat']);

% Map data
cpath = '/Volumes/GFDL/NEMURO/10km/';
load([cpath 'gridspec_10km.mat'],'LON','LAT');
load([cpath 'Data_grid_10km_hist.mat']);
[ni,nj]=size(LON);
plotminlat=30; %Set these bounds for your data
plotmaxlat=48;
plotminlon=-134;
plotmaxlon=-115;
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
figure(1)
subplot(4,1,1)
plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Historic 10km Large Pelagics')
legend('Larvae','Juveniles','Adults')
legend('location','southeast')
stamp(cfile)

subplot(4,1,2)
plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Larvae')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,3)
plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Juveniles')

subplot(4,1,4)
plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Adults')
xlabel('Time (mo)')
print('-dpng',[ppath 'Historic10km_',harv,'_P_time.png'])

% Forage fishes
figure(2)
subplot(3,1,1)
plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Historic 10km Forage Fishes')
ylabel('log10 Biomass (g m^-^2)')
legend('Immature','Adults')
legend('location','southeast')
stamp(cfile)

subplot(3,1,2)
plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Immature')
ylabel('log10 Biomass (g m^-^2)')

subplot(3,1,3)
plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Adults')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Historic10km_',harv,'_F_time.png'])

% Demersals
figure(3)
subplot(4,1,1)
plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Historic 10km Demersal Piscivores')
legend('Larvae','Juveniles','Adults')
legend('location','southeast')
stamp(cfile)

subplot(4,1,2)
plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Larvae')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,3)
plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Juveniles')

subplot(4,1,4)
plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Adults')
xlabel('Time (mo)')
print('-dpng',[ppath 'Historic10km_',harv,'_D_time.png'])

% Benthic inverts
figure(4)
subplot(2,1,1)
plot(y,log10(b_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Historic 10km Benthic Inverts')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')

subplot(2,1,2)
plot(y,(b_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
ylabel('Biomass (g m^-^2)')
print('-dpng',[ppath 'Historic10km_',harv,'_B_time.png'])

%% All size classes of all
figure(5)
plot(y,(sf_tmean),'Linewidth',1); hold on;
plot(y,(mf_tmean),'Linewidth',1); hold on;
plot(y,(sp_tmean),'Linewidth',1); hold on;
plot(y,(mp_tmean),'Linewidth',1); hold on;
plot(y,(lp_tmean),'Linewidth',1); hold on;
plot(y,(sd_tmean),'Linewidth',1); hold on;
plot(y,(md_tmean),'Linewidth',1); hold on;
plot(y,(ld_tmean),'Linewidth',1); hold on;
plot(y,(b_tmean),'Linewidth',1); hold on;
legend('SF','MF','SP','MP','LP','SD','MD','LD','B')
legend('location','eastoutside')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title(['Historic 10km ' tharv])
stamp(cfile)
print('-dpng',[ppath 'Historic10km_',harv,'_all_sizes.png'])

figure(6)
subplot(2,1,2)
plot(y,(sf_tmean),'Linewidth',1,'color',cm10(1,:)); hold on;
plot(y,(mf_tmean),'Linewidth',1,'color',cm10(2,:)); hold on;
plot(y,(sp_tmean),'Linewidth',1,'color',cm10(3,:)); hold on;
plot(y,(sd_tmean),'Linewidth',1,'color',cm10(6,:)); hold on;
plot(y,(md_tmean),'Linewidth',1,'color',cm10(7,:)); hold on;
plot(y,(ld_tmean),'Linewidth',1,'color',cm10(8,:)); hold on;
plot(y,(b_tmean),'Linewidth',1,'color',cm10(9,:)); hold on;
legend('SF','MF','SP','SD','MD','LD','B')
legend('location','eastoutside')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')

subplot(2,1,1)
plot(y,(mp_tmean),'r','Linewidth',1); hold on;
plot(y,(lp_tmean),'color',[0.5 0 0],'Linewidth',1); hold on;
xlim([y(1) y(end)])
legend('MP','LP')
legend('location','eastoutside')
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title(['Historic 10km ' tharv])
stamp(cfile)
print('-dpng',[ppath 'Historic10km_',harv,'_all_sizes_subplot.png'])

F = sf_tmean+mf_tmean;
P = sp_tmean+mp_tmean+lp_tmean;
D = sd_tmean+md_tmean+ld_tmean;
B = b_tmean;

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
title(['Historic 10km ' tharv])
stamp(cfile)
print('-dpng',[ppath 'Historic10km_',harv,'_all_types.png'])

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
title(['Historic 10km ' tharv])
stamp(cfile)
print('-dpng',[ppath 'Historic10km_',harv,'_all_types_subplot.png'])


%% Plots in space
[ni,nj]=size(LON);

ngrid = NaN*ones(ni,nj); 

Zsf=NaN*ones(ni*nj,23);
Zsp=NaN*ones(ni*nj,23);
Zsd=NaN*ones(ni*nj,23);
Zmf=NaN*ones(ni*nj,23);
Zmp=NaN*ones(ni*nj,23);
Zmd=NaN*ones(ni*nj,23);
Zlp=NaN*ones(ni*nj,23);
Zld=NaN*ones(ni*nj,23);
Zb=NaN*ones(ni*nj,23);

Zsf(GRD.ID,:)=mSF;
Zsp(GRD.ID,:)=mSP;
Zsd(GRD.ID,:)=mSD;
Zmf(GRD.ID,:)=mMF;
Zmp(GRD.ID,:)=mMP;
Zmd(GRD.ID,:)=mMD;
Zlp(GRD.ID,:)=mLP;
Zld(GRD.ID,:)=mLD;
Zb(GRD.ID,:)=mB;

Zsf = reshape(Zsf,ni,nj,23);
Zsp = reshape(Zsp,ni,nj,23);
Zsd = reshape(Zsd,ni,nj,23);
Zmf = reshape(Zmf,ni,nj,23);
Zmp = reshape(Zmp,ni,nj,23);
Zmd = reshape(Zmd,ni,nj,23);
Zlp = reshape(Zlp,ni,nj,23);
Zld = reshape(Zld,ni,nj,23);
Zb = reshape(Zb,ni,nj,23);

%% bent
figure(9)
for n=1:23
    subplot(5,5,n)
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','on','FLineWidth',1)
    surfm(LAT,LON,log10(squeeze(Zb(:,:,n))))
    colormap('jet')
    % load coast;                     %decent looking coastlines
    % h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-2 1]);
    set(gcf,'renderer','painters')
    if (n==3)
        title({'Benthic inverts (log_1_0 g m^-^2)'; num2str(1987+n)})
    else
        title(num2str(1987+n))
    end
end
colorbar('Position',[0.625 0.15 0.25 0.03],'orientation','horizontal')
stamp(cfile)
print('-dpng',[ppath 'Historic10km_',harv,'_CC_yrs_BENT.png'])

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
for n=1:23
    subplot(5,5,n)
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','on','FLineWidth',1)
    surfm(LAT,LON,log10(squeeze(All(:,:,n))))
    colormap('jet')
    % load coast;                     %decent looking coastlines
    % h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([0 2]);
    set(gcf,'renderer','painters')
    if (n==3)
        title({'All fishes (log_1_0 g m^-^2)'; num2str(1987+n)})
    else
        title(num2str(1987+n))
    end
end
colorbar('Position',[0.625 0.15 0.25 0.03],'orientation','horizontal')
stamp(cfile)
print('-dpng',[ppath 'Historic10km_',harv,'_CC_yrs_All.png'])

%% all F
figure(11)
for n=1:23
    subplot(5,5,n)
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','on','FLineWidth',1)
    surfm(LAT,LON,log10(squeeze(AllF(:,:,n))))
    colormap('jet')
    % load coast;                     %decent looking coastlines
    % h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-1 1]);
    set(gcf,'renderer','painters')
    if (n==3)
        title({'Forage fish (log_1_0 g m^-^2)'; num2str(1987+n)})
    else
        title(num2str(1987+n))
    end
end
colorbar('Position',[0.625 0.15 0.25 0.03],'orientation','horizontal')
stamp(cfile)
print('-dpng',[ppath 'Historic10km_',harv,'_CC_yrs_AllF.png'])

% all D
figure(12)
for n=1:23
    subplot(5,5,n)
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','on','FLineWidth',1)
    surfm(LAT,LON,log10(squeeze(AllD(:,:,n))))
    colormap('jet')
    % load coast;                     %decent looking coastlines
    % h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-1 1]);
    set(gcf,'renderer','painters')
    if (n==3)
        title({'Demersals (log_1_0 g m^-^2)'; num2str(1987+n)})
    else
        title(num2str(1987+n))
    end
end
colorbar('Position',[0.625 0.15 0.25 0.03],'orientation','horizontal')
stamp(cfile)
print('-dpng',[ppath 'Historic10km_',harv,'_CC_yrs_AllD.png'])

% All P
figure(13)
for n=1:23
    subplot(5,5,n)
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','on','FLineWidth',1)
    surfm(LAT,LON,log10(squeeze(AllP(:,:,n))))
    colormap('jet')
    % load coast;                     %decent looking coastlines
    % h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([0 2]);
    set(gcf,'renderer','painters')
    if (n==3)
        title({'Large pelagics (log_1_0 g m^-^2)'; num2str(1987+n)})
    else
        title(num2str(1987+n))
    end
end
colorbar('Position',[0.625 0.15 0.25 0.03],'orientation','horizontal')
stamp(cfile)
print('-dpng',[ppath 'Historic10km_',harv,'_CC_yrs_AllP.png'])


%% Ratios on subplots red-white-blue
% 3 figure subplot P:D, P:F, M:L
figure(14)
for n=1:23
    subplot(5,5,n)
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','on','FLineWidth',1)
    surfm(LAT,LON,squeeze(FracPD(:,:,n)))
    cmocean('balance')
    % load coast;                     %decent looking coastlines
    % h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([0 1]);
    set(gcf,'renderer','painters')
    if (n==3)
        title({'Fraction Large Pelagics vs. Demersals'; num2str(1987+n)})
    else
        title(num2str(1987+n))
    end
end
colorbar('Position',[0.625 0.15 0.25 0.03],'orientation','horizontal')
stamp(cfile)
print('-dpng',[ppath 'Historic10km_',harv,'_CC_yrs_PD.png'])

%P:F
figure(15)
for n=1:23
    subplot(5,5,n)
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','on','FLineWidth',1)
    surfm(LAT,LON,squeeze(FracPF(:,:,n)))
    cmocean('balance')
    % load coast;                     %decent looking coastlines
    % h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([0 1]);
    set(gcf,'renderer','painters')
    if (n==3)
        title({'Fraction Large Pelagics vs. Forage Fish'; num2str(1987+n)})
    else
        title(num2str(1987+n))
    end
end
colorbar('Position',[0.625 0.15 0.25 0.03],'orientation','horizontal')
stamp(cfile)
print('-dpng',[ppath 'Historic10km_',harv,'_CC_yrs_PF.png'])

%L:M
figure(16)
for n=1:23
    subplot(5,5,n)
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','on','FLineWidth',1)
    surfm(LAT,LON,squeeze(FracLM(:,:,n)))
    cmocean('balance')
    % load coast;                     %decent looking coastlines
    % h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([0 1]);
    set(gcf,'renderer','painters')
    if (n==3)
        title({'Fraction Large vs. Medium'; num2str(1987+n)})
    else
        title(num2str(1987+n))
    end
end
colorbar('Position',[0.625 0.15 0.25 0.03],'orientation','horizontal')
stamp(cfile)
print('-dpng',[ppath 'Historic10km_',harv,'_CC_yrs_LM.png'])

%% Debug
% Zsf2=NaN*ones(ni,nj);
% Zsp2=NaN*ones(ni,nj);
% Zsd2=NaN*ones(ni,nj);
% Zmf2=NaN*ones(ni,nj);
% Zmp2=NaN*ones(ni,nj);
% Zmd2=NaN*ones(ni,nj);
% Zlp2=NaN*ones(ni,nj);
% Zld2=NaN*ones(ni,nj);
% Zb2=NaN*ones(ni,nj);
% 
% Zsf2(GRD.ID)=sf_mean;
% Zsp2(GRD.ID)=sp_mean;
% Zsd2(GRD.ID)=sd_mean;
% Zmf2(GRD.ID)=mf_mean;
% Zmp2(GRD.ID)=mp_mean;
% Zmd2(GRD.ID)=md_mean;
% Zlp2(GRD.ID)=lp_mean;
% Zld2(GRD.ID)=ld_mean;
% Zb2(GRD.ID)=b_mean;
% 
% AllP2 = Zsp2+Zmp2+Zlp2;
% All2 = Zsp2+Zsf2+Zsd2+Zmp2+Zmf2+Zmd2+Zlp2+Zld2;
% AllF2 = Zsf2+Zmf2;
% AllD2 = Zsd2+Zmd2+Zld2;
% 
% testP=AllP(:,:,20);
% testF=AllF(:,:,20);
% testD=AllD(:,:,20);
% testA=All(:,:,20);
