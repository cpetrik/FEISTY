% Visualize output of FEISTY
% 10 cycles of 1961-1980 ctrlclim to spinup biomass
% Time series plots and maps

clear
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];

% mods = {'All_fish_obs','All_fish_obs_v1.2','All_fish_obs_v2',...
%     'All_fish_obs_v3','All_fish_obs_v3.2',...e lan   sharee 
%     'All_fishobs_assessment','All_fishobs_effective','All_fishobs_nominal'};

mods = {'FFmsy_creep','FFmsy_nominal','FFmsymax_creep','FFmsymax_nominal',...
    'FFmsymin_creep','FFmsymin_nominal'};

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP/Phase3a/';
ppath = [pp cfile '/OneDeg/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%% Map data
%cpath = '/Volumes/MIP/Fish-MIP/Phase3/OneDeg/';
cpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';
load([cpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat']);
load([cpath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat'], 'GRD');

[ni,nj]=size(LON);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
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

set(groot,'defaultAxesColorOrder',cm10);

load coastlines;                     %decent looking coastlines

%%
for i=1:length(mods)
    %harv = harvs{i};
    mod = ['All_fishobs_' mods{i}];
    %mod = [harv '_ctrlclim_onedeg'];
    
    %mod = mods{i};

    load([fpath 'Means_PI_ctrlclim_',mod,'_' cfile '.mat'],'time',...
        'sf_tmean','sp_tmean','sd_tmean',...
        'mf_tmean','mp_tmean','md_tmean',...
        'lp_tmean','ld_tmean','b_tmean',...
        'mf_tmcatch','mp_tmcatch','md_tmcatch',...
        'lp_tmcatch','ld_tmcatch',...
        'sf_smean','sp_smean','sd_smean',...
        'mf_smean','mp_smean','md_smean',...
        'lp_smean','ld_smean','b_smean',...
        'mf_mcatch','mp_mcatch','md_mcatch',...
        'lp_mcatch','ld_mcatch');
    %TOO MUCH to load everything, CHANGE to tmean and smeans


    %% Plots in time
    %t = 1:length(sp_tmean); %time;
    y = 1841 + (time)/12;

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
    ylim([-2 1])
    xlabel('Time (y)')
    ylabel('log_1_0 Biomass (g m^-^2)')
    title('PI')
    stamp(mod)
    print('-dpng',[ppath 'PI_empHP_',mod,'_all_sizes.png'])

    %%
    figure(2)
    plot(y,log10(mf_tmcatch),'Linewidth',1); hold on;
    plot(y,log10(mp_tmcatch),'Linewidth',1); hold on;
    plot(y,log10(lp_tmcatch),'Linewidth',1); hold on;
    plot(y,log10(md_tmcatch),'Linewidth',1); hold on;
    plot(y,log10(ld_tmcatch),'Linewidth',1); hold on;
    legend('MF','MP','LP','MD','LD')
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    ylim([-7.5 -3.5])
    xlabel('Time (y)')
    ylabel('log_1_0 Yield (g m^-^2)')
    title('PI')
    stamp(mod)
    print('-dpng',[ppath 'PI_empHP_',mod,'_yield_all_sizes.png'])

    %% Types together
    F = sf_tmean+mf_tmean;
    P = sp_tmean+mp_tmean+lp_tmean;
    D = sd_tmean+md_tmean+ld_tmean;
    B = b_tmean;

    yF = mf_tmcatch;
    yP = mp_tmcatch+lp_tmcatch;
    yD = md_tmcatch+ld_tmcatch;

    %%
    figure(3)
    plot(y,log10(B),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(F),'r','Linewidth',2); hold on;
    plot(y,log10(P),'b','Linewidth',2); hold on;
    plot(y,log10(D),'k','Linewidth',2); hold on;
    legend('B','F','P','D')
    legend('location','east')
    xlim([y(1) y(end)])
    ylim([-1 1])
    xlabel('Time (y)')
    ylabel('log_1_0 Biomass (g m^-^2)')
    title('PI')
    stamp(mod)
    print('-dpng',[ppath 'PI_empHP_',mod,'_all_types.png'])

    figure(4)
    plot(y,log10(yF),'r','Linewidth',2); hold on;
    plot(y,log10(yP),'b','Linewidth',2); hold on;
    plot(y,log10(yD),'k','Linewidth',2); hold on;
    legend('F','P','D')
    legend('location','east')
    xlim([y(1) y(end)])
    ylim([-6.5 -3.5])
    xlabel('Time (y)')
    ylabel('log_1_0 Yield (g m^-^2)')
    title('PI')
    stamp(mod)
    print('-dpng',[ppath 'PI_empHP_',mod,'_yield_all_types.png'])

    %% Plots in space
    Zsf=NaN*ones(ni,nj);
    Zsp=NaN*ones(ni,nj);
    Zsd=NaN*ones(ni,nj);
    Zmf=NaN*ones(ni,nj);
    Zmp=NaN*ones(ni,nj);
    Zmd=NaN*ones(ni,nj);
    Zlp=NaN*ones(ni,nj);
    Zld=NaN*ones(ni,nj);
    Zb=NaN*ones(ni,nj);

    Zsf(GRD.ID)=sf_smean;
    Zsp(GRD.ID)=sp_smean;
    Zsd(GRD.ID)=sd_smean;
    Zmf(GRD.ID)=mf_smean;
    Zmp(GRD.ID)=mp_smean;
    Zmd(GRD.ID)=md_smean;
    Zlp(GRD.ID)=lp_smean;
    Zld(GRD.ID)=ld_smean;
    Zb(GRD.ID)=b_smean;

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

    %%
    Cmf=NaN*ones(ni,nj);
    Cmp=NaN*ones(ni,nj);
    Cmd=NaN*ones(ni,nj);
    Clp=NaN*ones(ni,nj);
    Cld=NaN*ones(ni,nj);

    Cmf(GRD.ID)=mf_mcatch;
    Cmp(GRD.ID)=mp_mcatch;
    Cmd(GRD.ID)=md_mcatch;
    Clp(GRD.ID)=lp_mcatch;
    Cld(GRD.ID)=ld_mcatch;

    CAllF = Cmf;
    CAllP = Cmp+Clp;
    CAllD = Cmd+Cld;
    CAll = CAllF + CAllP + CAllD;

    %% bent
    figure(5)
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(Zb))
    cmocean('matter')
    h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-1 2]);
    hcb = colorbar('h');
    set(gcf,'renderer','painters')
    title('PI 1949 log10 mean benthic biomass (g m^-^2)')
    stamp(mod)
    print('-dpng',[ppath 'PI_empHP_',mod,'_global_BENT.png'])

    %% All 4 on subplots
    figure(6)
    % all F
    subplot('Position',[0 0.51 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(AllF))
    cmocean('matter')
    h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-2 2]);
    colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title('log10 mean All F (g m^-^2)')

    % all D
    subplot('Position',[0 0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(AllD))
    cmocean('matter')
    h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-2 2]);
    set(gcf,'renderer','painters')
    title('log10 mean All D (g m^-^2)')

    % All P
    subplot('Position',[0.5 0.51 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(AllP))
    cmocean('matter')
    h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-2 2]);
    set(gcf,'renderer','painters')
    title('log10 mean All P (g m^-^2)')

    % All
    subplot('Position',[0.5 0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(All))
    cmocean('matter')
    h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-2 2]);
    set(gcf,'renderer','painters')
    title('log10 mean All fishes (g m^-^2)')
    stamp(mod)
    print('-dpng',[ppath 'PI_empHP_',mod,'_global_All_subplot.png'])

    %% Ratios on subplots red-white-blue
    % 3 figure subplot P:D, P:F, M:L
%     figure(7)
%     subplot('Position',[0 0.53 0.5 0.5])
%     %P:D
%     axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%         'Grid','off','FLineWidth',1)
%     surfm(LAT,LON,FracPD)
%     cmocean('balance')
%     h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
%     caxis([0 1]);
%     set(gcf,'renderer','painters')
%     title('Fraction Large Pelagics vs. Demersals')
% 
%     %P:F
%     subplot('Position',[0.5 0.53 0.5 0.5])
%     axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%         'Grid','off','FLineWidth',1)
%     surfm(LAT,LON,FracPF)
%     cmocean('balance')
%     h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
%     caxis([0 1]);
%     set(gcf,'renderer','painters')
%     title('Fraction Large Pelagics vs. Forage Fishes')
% 
%     %L:M
%     subplot('Position',[0.25 0.0 0.5 0.5])
%     axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%         'Grid','off','FLineWidth',1)
%     surfm(LAT,LON,FracLM)
%     cmocean('balance')
%     h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
%     caxis([0 1]);
%     colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
%     set(gcf,'renderer','painters')
%     title('Fraction Large vs. Medium')
%     stamp(mod)
%     print('-dpng',[ppath 'PI_empHP_',mod,'_global_ratios_subplot.png'])

    %% All 4 Catch on subplots
    figure(8)
    % all F
    subplot('Position',[0 0.51 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(CAllF))
    cmocean('matter')
    h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-6 -2])
    colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title('log10 mean All F (g m^-^2)')

    % all D
    subplot('Position',[0 0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(CAllD))
    cmocean('matter')
    h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-6 -2])
    set(gcf,'renderer','painters')
    title('log10 mean All D (g m^-^2)')

    % All P
    subplot('Position',[0.5 0.51 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(CAllP))
    cmocean('matter')
    h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-6 -2])
    set(gcf,'renderer','painters')
    title('log10 mean All P (g m^-^2)')

    % All
    subplot('Position',[0.5 0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(CAll))
    cmocean('matter')
    h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-6 -2])
    set(gcf,'renderer','painters')
    title('log10 mean All fishes (g m^-^2)')
    stamp(mod)
    print('-dpng',[ppath 'PI_empHP_',mod,'_global_yield_All_subplot.png'])

end

