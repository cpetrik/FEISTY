% Visualize output of Spinup
% 50 years using no move IC
% Saved as mat files

clear
close all

%%
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/CORE/'];
ppath = [pp cfile '/CORE/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%%
load('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');
load('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

[ni,nj]=size(geolon_t);
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90;
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
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

%%
behave = {'enc','nu','mort','prey','preyconc'};

for b = 5 %1:length(behave)

    %%
    exper = ['Spinup1988_move_',behave{b},'_v27_dt12h_All_fish03_'];
    if b==4
        exper2 = 'Spinup1988_move_ingest_v27_dt12h_All_fish03_';
    else
        exper2 = exper;
    end

    load([fpath 'Means_' exper2 cfile '.mat']);

    close all

    %% Plots in time
    y = time/12;
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
    ylim([-3 2])
    xlabel('Time (mo)')
    ylabel('log10 Biomass (g m^-^2)')
    title('Spinup')
    stamp(exper2)
    print('-dpng',[ppath exper2 'all_sizes.png'])

    figure(2)
    plot(y,log10(B),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(F),'r','Linewidth',2); hold on;
    plot(y,log10(P),'b','Linewidth',2); hold on;
    plot(y,log10(D),'k','Linewidth',2); hold on;
    legend('B','F','P','D')
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    ylim([-2 2])
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    title(['Spinup'])
    print('-dpng',[ppath exper2 'all_types.png'])

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

    Zsf(GRD.ID)=sf_mean;
    Zsp(GRD.ID)=sp_mean;
    Zsd(GRD.ID)=sd_mean;
    Zmf(GRD.ID)=mf_mean;
    Zmp(GRD.ID)=mp_mean;
    Zmd(GRD.ID)=md_mean;
    Zlp(GRD.ID)=lp_mean;
    Zld(GRD.ID)=ld_mean;
    Zb(GRD.ID)=b_mean;

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

    %% bent
    figure(3)
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,log10(Zb))
    cmocean('dense')
    load coastlines;
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-2 2]);
    hcb = colorbar('h');
    set(gcf,'renderer','painters')
    title('log10 mean benthic biomass (g m^-^2)')
    stamp(exper2)
    print('-dpng',[ppath exper2 'global_BENT.png'])

    %% ALL
    figure(4)
    % all F
    subplot('Position',[0 0.51 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,log10(AllF))
    cmocean('dense')
    load coastlines;
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-2 2]);
    colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title('log10 mean All F (g m^-^2)')

    % all D
    subplot('Position',[0 0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,log10(AllD))
    cmocean('dense')
    load coastlines;
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-2 2]);
    set(gcf,'renderer','painters')
    title('log10 mean All D (g m^-^2)')

    % All P
    subplot('Position',[0.5 0.51 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,log10(AllP))
    cmocean('dense')
    load coastlines;
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-2 2]);
    set(gcf,'renderer','painters')
    title('log10 mean All P (g m^-^2)')

    % All
    subplot('Position',[0.5 0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,log10(All))
    cmocean('dense')
    load coastlines;
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-2 2]);
    set(gcf,'renderer','painters')
    title('log10 mean All fishes (g m^-^2)')
    stamp(exper2)
    print('-dpng',[ppath exper2 'All_subplot.png'])

    %% Ratios on subplots 3 figure subplot P:D, P:F, M:L
    figure(5)
    subplot('Position',[0 0.53 0.5 0.5])
    %P:D
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,FracPD)
    cmocean('balance')
    load coastlines;
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([0 1]);
    set(gcf,'renderer','painters')
    title('Fraction Large Pelagics vs. Demersals')

    %P:F
    subplot('Position',[0.5 0.53 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,FracPF)
    cmocean('balance')
    load coastlines;
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([0 1]);
    set(gcf,'renderer','painters')
    title('Fraction Large Pelagics vs. Forage Fishes')

    %L:M
    subplot('Position',[0.25 0.0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,FracLM)
    cmocean('balance')
    load coastlines;
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([0 1]);
    colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title('Fraction Large vs. Medium')
    stamp(exper2)
    print('-dpng',[ppath exper2 'global_ratios_subplot.png'])

end

