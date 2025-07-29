% Visualize output of FEISTY
% SSP 534 2040-2300
% Time series plots and maps

clear
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

esms = {'IPSL','UKESM','CESM2-WACCM','CESM2-WACCM'};
e2 = {'ipsl','ukesm','cesm2','cesm2'};

pp = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/FishMIP/wg2300/';
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

cmBP50=cbrewer('seq','BuPu',50,'PCHIP');

for m=1 %1:length(esms)

    close all
    
    mod = esms{m};
    mod2 = e2{m};

    if m==3
        exper = [mod '_ssp534_zooc_pristine'];
    elseif m==4
        exper = [mod '_ssp534_zmeso_pristine'];
    elseif m==2
        exper = [mod '_ssp534_pristine'];
    else
        exper = [mod '_ssp534-over_pristine'];
    end
    
    fpath=['/Volumes/petrik-lab/Feisty/NC/WG2300/',cfile,'/',mod,'/'];

    load([fpath 'Means_' exper '_' cfile '.mat']);

    % Map data
    if m==2
        cpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/';
    else
        cpath = ['/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/',mod,'/'];
    end
    load([cpath 'gridspec_',mod2,'_cmip6_2300.mat']);
    load([cpath 'Data_grid_',mod2,'_cmip6_2300.mat']);

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
    
    %% Plots in time
    t = 1:length(sp_tmean); %time;
    y = 2040 + (t-1)/12;

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
    ylim([-3 1])
    xlabel('Time (y)')
    ylabel('log_1_0 Biomass (g m^-^2)')
    title([mod ' SSP534'])
    stamp(exper)
    print('-dpng',[ppath exper,'_ts_all_sizes.png'])

    %% Types together
    F = sf_tmean+mf_tmean;
    P = sp_tmean+mp_tmean+lp_tmean;
    D = sd_tmean+md_tmean+ld_tmean;
    B = b_tmean;

    figure(2)
    plot(y,log10(B),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(F),'r','Linewidth',2); hold on;
    plot(y,log10(P),'b','Linewidth',2); hold on;
    plot(y,log10(D),'k','Linewidth',2); hold on;
    legend('B','F','P','D')
    legend('location','northeast')
    xlim([y(1) y(end)])
    ylim([-2 2])
    xlabel('Time (y)')
    ylabel('log_1_0 Biomass (g m^-^2)')
    title([mod ' SSP534'])
    stamp(exper)
    print('-dpng',[ppath exper,'_ts_all_types.png'])

    %% Plots in space

    % SPACE MEANS HAVE NANS BUT NOT TIME MEANS?! 

    %2100 means
    Zsf=NaN*ones(ni,nj);
    Zsp=NaN*ones(ni,nj);
    Zsd=NaN*ones(ni,nj);
    Zmf=NaN*ones(ni,nj);
    Zmp=NaN*ones(ni,nj);
    Zmd=NaN*ones(ni,nj);
    Zlp=NaN*ones(ni,nj);
    Zld=NaN*ones(ni,nj);
    Zb=NaN*ones(ni,nj);

    Zsf(GRD.ID)=sf_mean1;
    Zsp(GRD.ID)=sp_mean1;
    Zsd(GRD.ID)=sd_mean1;
    Zmf(GRD.ID)=mf_mean1;
    Zmp(GRD.ID)=mp_mean1;
    Zmd(GRD.ID)=md_mean1;
    Zlp(GRD.ID)=lp_mean1;
    Zld(GRD.ID)=ld_mean1;
    Zb(GRD.ID)=b_mean1;

    All1 = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
    AllF1 = Zsf+Zmf;
    AllP1 = Zsp+Zmp+Zlp;
    AllD1 = Zsd+Zmd+Zld;
    AllS1 = Zsp+Zsf+Zsd;
    AllM1 = Zmp+Zmf+Zmd;
    AllL1 = Zlp+Zld;
    FracPD1 = AllP1 ./ (AllP1+AllD1);
    FracPF1 = AllP1 ./ (AllP1+AllF1);
    FracLM1 = AllL1 ./ (AllL1+AllM1);

    %2300 means
    Esf=NaN*ones(ni,nj);
    Esp=NaN*ones(ni,nj);
    Esd=NaN*ones(ni,nj);
    Emf=NaN*ones(ni,nj);
    Emp=NaN*ones(ni,nj);
    Emd=NaN*ones(ni,nj);
    Elp=NaN*ones(ni,nj);
    Eld=NaN*ones(ni,nj);
    Eb=NaN*ones(ni,nj);

    Esf(GRD.ID)=sf_mean3;
    Esp(GRD.ID)=sp_mean3;
    Esd(GRD.ID)=sd_mean3;
    Emf(GRD.ID)=mf_mean3;
    Emp(GRD.ID)=mp_mean3;
    Emd(GRD.ID)=md_mean3;
    Elp(GRD.ID)=lp_mean3;
    Eld(GRD.ID)=ld_mean3;
    Eb(GRD.ID)=b_mean3;

    All3 = Esp+Esf+Esd+Emp+Emf+Emd+Elp+Eld;
    AllF3 = Esf+Emf;
    AllP3 = Esp+Emp+Elp;
    AllD3 = Esd+Emd+Eld;
    AllS3 = Esp+Esf+Esd;
    AllM3 = Emp+Emf+Emd;
    AllL3 = Elp+Eld;
    FracPD3 = AllP3 ./ (AllP3+AllD3);
    FracPF3 = AllP3 ./ (AllP3+AllF3);
    FracLM3 = AllL3 ./ (AllL3+AllM3);

    %% bent
    figure(3)
    subplot('Position',[0 0.51 0.9 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(Zb))
    colormap(cmBP50)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-1 2]);
    title('SSP534 2100 log10 mean benthic biomass (g m^-^2)')
    
    subplot('Position',[0 0 0.9 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(Eb))
    colormap(cmBP50)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-1 2]);
    hcb = colorbar('h');
    set(gcf,'renderer','painters')
    title('SSP534 2300 log10 mean benthic biomass (g m^-^2)')
    stamp(exper)
    print('-dpng',[ppath exper,'_global_BENT_2100_2300.png'])

    %% All 4 on subplots
    figure(4)
    % all F
    subplot('Position',[0 0.51 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(AllF1))
    colormap(cmBP50)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-2 2]);
    colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title('log10 mean All F (g m^-^2) 2100')

    % all D
    subplot('Position',[0 0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(AllD1))
    colormap(cmBP50)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-2 2]);
    set(gcf,'renderer','painters')
    title('log10 mean All D (g m^-^2)')

    % All P
    subplot('Position',[0.5 0.51 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(AllP1))
    colormap(cmBP50)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-2 2]);
    set(gcf,'renderer','painters')
    title('log10 mean All P (g m^-^2)')

    % All
    subplot('Position',[0.5 0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(All1))
    colormap(cmBP50)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-2 2]);
    set(gcf,'renderer','painters')
    title('log10 mean All fishes (g m^-^2)')
    stamp(exper)
    print('-dpng',[ppath exper,'_global_All_subplot_2100.png'])

    %% All 4 on subplots
    figure(14)
    % all F
    subplot('Position',[0 0.51 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(AllF3))
    colormap(cmBP50)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-2 2]);
    colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title('log10 mean All F (g m^-^2) 2300')

    % all D
    subplot('Position',[0 0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(AllD3))
    colormap(cmBP50)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-2 2]);
    set(gcf,'renderer','painters')
    title('log10 mean All D (g m^-^2)')

    % All P
    subplot('Position',[0.5 0.51 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(AllP3))
    colormap(cmBP50)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-2 2]);
    set(gcf,'renderer','painters')
    title('log10 mean All P (g m^-^2)')

    % All
    subplot('Position',[0.5 0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(All3))
    colormap(cmBP50)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-2 2]);
    set(gcf,'renderer','painters')
    title('log10 mean All fishes (g m^-^2)')
    stamp(exper)
    print('-dpng',[ppath exper,'_global_All_subplot_2300.png'])

    %% Ratios on subplots red-white-blue
    % 3 figure subplot P:D, P:F, M:L
    figure(5)
    subplot('Position',[0 0.53 0.5 0.5])
    %P:D
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,FracPD1)
    cmocean('balance')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([0 1]);
    set(gcf,'renderer','painters')
    title('Fraction Large Pelagics vs. Demersals 2100')

    %P:F
    subplot('Position',[0.5 0.53 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,FracPF1)
    cmocean('balance')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([0 1]);
    set(gcf,'renderer','painters')
    title('Fraction Large Pelagics vs. Forage Fishes')

    %L:M
    subplot('Position',[0.25 0.0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,FracLM1)
    cmocean('balance')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([0 1]);
    colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title('Fraction Large vs. Medium')
    stamp(exper)
    print('-dpng',[ppath exper,'_global_ratios_subplot_2100.png'])


    %% Ratios on subplots red-white-blue
    % 3 figure subplot P:D, P:F, M:L
    figure(15)
    subplot('Position',[0 0.53 0.5 0.5])
    %P:D
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,FracPD3)
    cmocean('balance')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([0 1]);
    set(gcf,'renderer','painters')
    title('Fraction Large Pelagics vs. Demersals 2300')

    %P:F
    subplot('Position',[0.5 0.53 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,FracPF3)
    cmocean('balance')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([0 1]);
    set(gcf,'renderer','painters')
    title('Fraction Large Pelagics vs. Forage Fishes')

    %L:M
    subplot('Position',[0.25 0.0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,FracLM3)
    cmocean('balance')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([0 1]);
    colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title('Fraction Large vs. Medium')
    stamp(exper)
    print('-dpng',[ppath exper,'_global_ratios_subplot_2300.png'])


end