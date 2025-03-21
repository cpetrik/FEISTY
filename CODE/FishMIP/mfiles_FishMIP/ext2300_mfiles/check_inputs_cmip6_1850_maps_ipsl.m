% Compare daily interp forcings across ESMs

clear
close all

%% Paths

esms = {'IPSL','UKESM'};
e2 = {'ipsl','ukesm'};

ppath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/FishMIP/wg2300/';

for m=3%1:length(esms)

    close all

    mod = esms{m};
    mod2 = e2{m};
    exper = [mod '_spinup_pristine'];

    % Map data
    if m==3
        cpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/';
        hpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
    else
        cpath = ['/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/',mod,'/'];
        hpath = ['/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/',mod,'/hist/'];
    end
    load([cpath 'gridspec_',mod2,'_cmip6_2300.mat']);
    load([cpath 'Data_grid_',mod2,'_cmip6_2300.mat']);
    load([hpath 'Data_',mod2,'_hist_daily_1850.mat'])

    [ni,nj]=size(LON);

    plotminlat=-90; %Set these bounds for your data
    plotmaxlat=90;
    plotminlon=-180;
    plotmaxlon=180;
    clatlim=[plotminlat plotmaxlat];
    clonlim=[plotminlon plotmaxlon];

    load coastlines

    %%
    ESM.dZm = 10 .^ (-2.925 + 1.964.*log10(ESM.Zm+eps) + 1.958e-2.*ESM.Tp);

    c_Tp = mean(ESM.Tp,2);
    c_Tb = mean(ESM.Tb,2);
    c_D = mean(ESM.det,2);
    c_Z = mean(ESM.Zm,2);
    c_Zl = mean(ESM.dZm,2);

    [mi,mj]=size(LON);
    cTp=NaN*ones(mi,mj);
    cTb=NaN*ones(mi,mj);
    cD=NaN*ones(mi,mj);
    cZ=NaN*ones(mi,mj);
    cZl=NaN*ones(mi,mj);

    cTp(GRD.ID)=c_Tp;
    cTb(GRD.ID)=c_Tb;
    cD(GRD.ID) =c_D;
    cZ(GRD.ID) =c_Z;
    cZl(GRD.ID) =c_Zl;

    %%
    f1 = figure('Units','inches','Position',[1 3 5 5]);

    subplot('Position',[0.05 0.68 0.45 0.3])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,cTp)
    cmocean('thermal')
    clim([0 35])
    colorbar
    title([mod ' Tp'])
    patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    subplot('Position',[0.5 0.68 0.45 0.3])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,cTb)
    cmocean('thermal')
    clim([0 35])
    colorbar
    title([mod ' Tb'])
    h2=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    subplot('Position',[0.05 0.37 0.45 0.3])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(cZ))
    cmocean('tempo')
    clim([0 2])
    colorbar
    title([mod ' Zmeso (gWW/m2)'])
    h3=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    subplot('Position',[0.5 0.37 0.45 0.3])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(cD))
    cmocean('tempo')
    clim([-2.5 1.5])
    colorbar
    title([mod ' Det (gWW/m2/d)'])
    h4=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    subplot('Position',[0.05 0.06 0.45 0.3])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,log10(cZl))
    cmocean('tempo')
    clim([-2 1])
    colorbar
    title([mod ' Zmeso loss (gWW/m2/d)'])
    h5=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    print('-dpng',[ppath 'Map_',mod,'_1850_daily_interp_forcings.png'])

end

