% Movement maps with jet colorbar for comp to Watson

clear 
close all

%%
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/CORE/'];
ppath = [pp cfile '/CORE/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%%
exper = 'CORE_Hindcast1988_no_move_';
load([fpath 'Means_' exper cfile '.mat']);

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

Zsf(GRD.ID)=sf_mean20;
Zsp(GRD.ID)=sp_mean20;
Zsd(GRD.ID)=sd_mean20;
Zmf(GRD.ID)=mf_mean20;
Zmp(GRD.ID)=mp_mean20;
Zmd(GRD.ID)=md_mean20;
Zlp(GRD.ID)=lp_mean20;
Zld(GRD.ID)=ld_mean20;
Zb(GRD.ID)=b_mean20;

All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
AllF = Zsf+Zmf;
AllP = Zsp+Zmp+Zlp;
AllD = Zsd+Zmd+Zld;

%% Only P
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP))
colormap('jet')
colorbar
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All P (g m^-^2)')

print('-dpng',[ppath exper 'global_LgPel_jet.png'])

%%
mods = {'enc','nu','mort','ingest','preyconc'};

for m =1:length(mods)

    close all
    exper = ['CORE_Hindcast_move_',mods{m},'_v27_dt12h_All_fish03_spin20IC'];

    load([fpath 'Means_' exper '_' cfile '.mat']);


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

    Zsf(GRD.ID)=sf_mean20;
    Zsp(GRD.ID)=sp_mean20;
    Zsd(GRD.ID)=sd_mean20;
    Zmf(GRD.ID)=mf_mean20;
    Zmp(GRD.ID)=mp_mean20;
    Zmd(GRD.ID)=md_mean20;
    Zlp(GRD.ID)=lp_mean20;
    Zld(GRD.ID)=ld_mean20;
    Zb(GRD.ID)=b_mean20;

    All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
    AllF = Zsf+Zmf;
    AllP = Zsp+Zmp+Zlp;
    AllD = Zsd+Zmd+Zld;
   

    %% All P
    clf
    figure(3)
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,log10(AllP))
    colormap('jet')
    colorbar
    load coastlines;
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-2 2]);
    set(gcf,'renderer','painters')
    title('log10 mean All P (g m^-^2)')

    print('-dpng',[ppath exper 'global_LgPel_jet.png'])


end
