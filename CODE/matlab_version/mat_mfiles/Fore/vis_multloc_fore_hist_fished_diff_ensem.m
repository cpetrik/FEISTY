% Map best and worst for each functional type

clear all
close all

nfile = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
load([nfile 'LHS_param5_mid5_bestAIC_params_multFup_neg_multPneg.mat'],'params');
load([nfile 'simnames_ensem5_mid5_bestAIC_multFup_multPneg.mat']);

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' ...
    'param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'...
    'full_runs/'];

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

%% Original parameters
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath = ['/Volumes/GFDL/NC/Matlab_new_size/',...
    'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/'];
load([fpath 'Time_Means_Historic_Forecast_',harv,'_' cfile '.mat']);

%% Ensemble parameter sets
epath = ['/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/',...
    'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([epath 'Historic_All_fish03_ensem5_mid5_bestAIC_multFup_multPneg.mat']);
load([epath 'Forecast_All_fish03_ensem5_mid5_bestAIC_multFup_multPneg.mat']);

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%% FEISTY Orig Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';

% Hindcast
load([fpath 'Means_Historic_' harv '_' cfile '.mat'],...
    'sf_mean50','sp_mean50','sd_mean50',...
    'mf_mean50','mp_mean50','md_mean50',...
    'lp_mean50','ld_mean50','b_mean50');

HF = sf_mean50 +mf_mean50 ;
HP = sp_mean50 +mp_mean50 +lp_mean50 ;
HD = sd_mean50 +md_mean50 +ld_mean50 ;

[hi,hj]=size(geolon_t);
hF=NaN*ones(hi,hj);
hP=NaN*ones(hi,hj);
hD=NaN*ones(hi,hj);
hB =NaN*ones(hi,hj);
hF(ID) = HF;
hP(ID) = HP;
hD(ID) = HD;
hB(ID) = b_mean50;

clear sf_mean50 sp_mean50 sd_mean50 mf_mean50 mp_mean50 md_mean50 lp_mean50 ld_mean50 b_mean50

% Forecast
load([fpath 'Means_fore_',harv,'_' cfile '.mat'],...
    'sf_mean50','sp_mean50','sd_mean50',...
    'mf_mean50','mp_mean50','md_mean50',...
    'lp_mean50','ld_mean50','b_mean50');

FF = sf_mean50 +mf_mean50 ;
FP = sp_mean50 +mp_mean50 +lp_mean50 ;
FD = sd_mean50 +md_mean50 +ld_mean50 ;

[ni,nj]=size(geolon_t);
cF=NaN*ones(ni,nj);
cP=NaN*ones(ni,nj);
cD=NaN*ones(ni,nj);
cB=NaN*ones(ni,nj);
cF(ID) = FF;
cP(ID) = FP;
cD(ID) = FD;
cB(ID) = b_mean50;

clear sf_mean50 sp_mean50 sd_mean50 mf_mean50 mp_mean50 md_mean50 lp_mean50 ld_mean50 b_mean50

%%
cAll = cF+cP+cD;
hAll = hF+hP+hD;

pdiffF = (cF-hF) ./ hF;
pdiffP = (cP-hP) ./ hP;
pdiffD = (cD-hD) ./ hD;
pdiffB = (cB-hB) ./ hB;
pdiffAll = (cAll-hAll) ./ hAll;

% diffF(hF(:)<1e-6) = nan;
% diffP(hP(:)<1e-6) = nan;
% diffD(hD(:)<1e-6) = nan;
% diffB(Hb(:)<1e-6) = nan;
% diffAll(hAll(:)<1e-6) = nan;

%% Maps
% All 3 F on subplots
figure(1)
% hist F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(hF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar('Position',[0.05 0.51 0.4 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
title('F 1951-2000')

% fore F
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(cF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('F 2051-2100')

% pdiff F
subplot('Position',[0.5 0.25 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffF))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar('Position',[0.55 0.3 0.4 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
title('F Percent Difference')
stamp('Orig')
print('-dpng',[pp 'Hist_Fore_',harv,'_global_F_pdiff_subplot.png'])

%%
% All 3 P on subplots
figure(2)
% hist P
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(hP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar('Position',[0.05 0.51 0.4 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
title('P 1951-2000')

% fore P
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(cP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('P 2051-2100')

% pdiff P
subplot('Position',[0.5 0.25 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffP))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar('Position',[0.55 0.3 0.4 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
title('P Percent Difference')
stamp('Orig')
print('-dpng',[pp 'Hist_Fore_',harv,'_global_P_pdiff_subplot.png'])

%%
% All 3 D on subplots
figure(3)
% hist D
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(hD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar('Position',[0.05 0.51 0.4 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
title('D 1951-2000')

% fore D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(cD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('D 2051-2100')

% pdiff D
subplot('Position',[0.5 0.25 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffD))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar('Position',[0.55 0.3 0.4 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
title('D Percent Difference')
stamp('Orig')
print('-dpng',[pp 'Hist_Fore_',harv,'_global_D_pdiff_subplot.png'])

%%
% All 3 all on subplots
figure(4)
% hist All
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(hAll))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar('Position',[0.05 0.51 0.4 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
title('All 1951-2000')

% fore All
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(cAll))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('All 2051-2100')

% pdiff All
subplot('Position',[0.5 0.25 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(pdiffAll))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar('Position',[0.55 0.3 0.4 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
title('All Percent Difference')
stamp('Orig')
print('-dpng',[pp 'Hist_Fore_',harv,'_global_All_pdiff_subplot.png'])


%% -----------------------------------------------------------------------
%  Ensemble members with most and least of diff fishes
%  -----------------------------------------------------------------------
%sims = [3;4;7;11;12;13];
sims = [4;8;10];
for i=1:length(sims)
    k = sims(i);
    
    %% Historic
    HF = hSsF(:,k) + hSmF(:,k);
    HP = hSsP(:,k) + hSmP(:,k) + hSlP(:,k);
    HD = hSsD(:,k) + hSmD(:,k) + hSlD(:,k);
    
    [hi,hj] = size(geolon_t);
    hF = NaN*ones(hi,hj);
    hP = NaN*ones(hi,hj);
    hD = NaN*ones(hi,hj);
    hB = NaN*ones(hi,hj);
    hF(ID) = HF;
    hP(ID) = HP;
    hD(ID) = HD;
    hB(ID) = hSB(:,k);
    
    % Forecast
    FF = fSsF(:,k) + fSmF(:,k);
    FP = fSsP(:,k) + fSmP(:,k) + fSlP(:,k);
    FD = fSsD(:,k) + fSmD(:,k) + fSlD(:,k);
    
    [ni,nj] = size(geolon_t);
    cF = NaN*ones(ni,nj);
    cP = NaN*ones(ni,nj);
    cD = NaN*ones(ni,nj);
    cB = NaN*ones(ni,nj);
    cF(ID) = FF;
    cP(ID) = FP;
    cD(ID) = FD;
    cB(ID) = fSB(:,k);
    
    %%
    cAll = cF+cP+cD;
    hAll = hF+hP+hD;
    
    pdiffF = (cF-hF) ./ hF;
    pdiffP = (cP-hP) ./ hP;
    pdiffD = (cD-hD) ./ hD;
    pdiffB = (cB-hB) ./ hB;
    pdiffAll = (cAll-hAll) ./ hAll;
    
    % diffF(hF(:)<1e-6) = nan;
    % diffP(hP(:)<1e-6) = nan;
    % diffD(hD(:)<1e-6) = nan;
    % diffB(Hb(:)<1e-6) = nan;
    % diffAll(hAll(:)<1e-6) = nan;
    
    %% Maps
    % All 3 F on subplots
    figure(1)
    clf
    % hist F
    subplot('Position',[0 0.51 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,log10(hF))
    colormap('jet')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-2 2]);
    colorbar('Position',[0.05 0.51 0.4 0.025],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title('F 1951-2000')
    
    % fore F
    subplot('Position',[0 0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,log10(cF))
    colormap('jet')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-2 2]);
    set(gcf,'renderer','painters')
    title('F 2051-2100')
    
    % pdiff F
    subplot('Position',[0.5 0.25 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,(pdiffF))
    cmocean('balance')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-0.5 0.5]);
    colorbar('Position',[0.55 0.3 0.4 0.025],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title('F Percent Difference')
    stamp(num2str(k))
    print('-dpng',[ppath 'Hist_Fore_',harv,'_global_F_pdiff_subplot_',snames{k},'.png'])
    
    %%
    % All 3 P on subplots
    figure(2)
    clf
    % hist P
    subplot('Position',[0 0.51 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,log10(hP))
    colormap('jet')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-2 2]);
    colorbar('Position',[0.05 0.51 0.4 0.025],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title('P 1951-2000')
    
    % fore P
    subplot('Position',[0 0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,log10(cP))
    colormap('jet')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-2 2]);
    set(gcf,'renderer','painters')
    title('P 2051-2100')
    
    % pdiff P
    subplot('Position',[0.5 0.25 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,(pdiffP))
    cmocean('balance')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-0.5 0.5]);
    colorbar('Position',[0.55 0.3 0.4 0.025],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title('P Percent Difference')
    stamp(num2str(k))
    print('-dpng',[ppath 'Hist_Fore_',harv,'_global_P_pdiff_subplot_',snames{k},'.png'])
    
    %%
    % All 3 D on subplots
    figure(3)
    clf
    % hist D
    subplot('Position',[0 0.51 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,log10(hD))
    colormap('jet')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-2 2]);
    colorbar('Position',[0.05 0.51 0.4 0.025],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title('D 1951-2000')
    
    % fore D
    subplot('Position',[0 0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,log10(cD))
    colormap('jet')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-2 2]);
    set(gcf,'renderer','painters')
    title('D 2051-2100')
    
    % pdiff D
    subplot('Position',[0.5 0.25 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,(pdiffD))
    cmocean('balance')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-0.5 0.5]);
    colorbar('Position',[0.55 0.3 0.4 0.025],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title('D Percent Difference')
    stamp(num2str(k))
    print('-dpng',[ppath 'Hist_Fore_',harv,'_global_D_pdiff_subplot_',snames{k},'.png'])
    
    %%
    % All 3 all on subplots
    figure(4)
    clf
    % hist All
    subplot('Position',[0 0.51 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,log10(hAll))
    colormap('jet')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-2 2]);
    colorbar('Position',[0.05 0.51 0.4 0.025],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title('All 1951-2000')
    
    % fore All
    subplot('Position',[0 0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,log10(cAll))
    colormap('jet')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-2 2]);
    set(gcf,'renderer','painters')
    title('All 2051-2100')
    
    % pdiff All
    subplot('Position',[0.5 0.25 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,(pdiffAll))
    cmocean('balance')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-0.5 0.5]);
    colorbar('Position',[0.55 0.3 0.4 0.025],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title('All Percent Difference')
    stamp(num2str(k))
    print('-dpng',[ppath 'Hist_Fore_',harv,'_global_All_pdiff_subplot_',snames{k},'.png'])
    
end
