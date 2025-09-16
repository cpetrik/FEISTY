% Visualize output of Spinup 
% 50 years using no-movement spinup as ICs
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
exper = 'Spinup1988_move_mort_v21_dt6h_All_fish03_';
load([fpath 'Means_' exper cfile '.mat'],'time','lyr',...
    'sf_smean','sp_smean','sd_smean',...
    'mf_smean','mp_smean','md_smean',...
    'b_smean','lp_smean','ld_smean');

%%
allF = sf_smean+mf_smean;
allP = sp_smean+mp_smean+lp_smean;
allD = sd_smean+md_smean+ld_smean;
allS = sf_smean+sp_smean+sd_smean;
allM = mf_smean+mp_smean+md_smean;
allL = lp_smean+ld_smean;
allFish = allF+allP+allD;

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

yr = 1:(length(time)/12);

%% F
ZF = NaN*ones(ni,nj);
ZF(GRD.ID) = allF(:,1);

figure(1)
h = imagescn(geolat_t',geolon_t',log10(ZF)');
cb = colorbar;
ylabel(cb,'log10 abund (g m^-^2)')
cmocean dense
title('Forage fish')
%title(datestr(datenum(0,0,1),'dd'))
clim([-2 2])

hold on
% he = earthimage;
% uistack(he,'bottom')

gif([ppath exper 'F.gif'],'frame',gcf,'delaytime',1/12,'nodither')

for k=2:50
    ZF = NaN*ones(ni,nj);
    ZF(GRD.ID) = allF(:,k);

    yd = yr(k);
    h.CData = real(log10(ZF)');
    title(['Year ' num2str(k)])
    gif
end

%% P
ZP = NaN*ones(ni,nj);
ZP(GRD.ID) = allP(:,1);

figure(2)
h = imagescn(geolat_t',geolon_t',log10(ZP)');
cb = colorbar;
ylabel(cb,'log10 abund (g m^-^2)')
cmocean dense
title('Forage fish')
%title(datestr(datenum(0,0,1),'dd'))
clim([-2 2])

hold on
% he = earthimage;
% uistack(he,'bottom')

gif([ppath exper 'P.gif'],'frame',gcf,'delaytime',1/12,'nodither')

for k=2:50
    ZP = NaN*ones(ni,nj);
    ZP(GRD.ID) = allP(:,k);

    yd = yr(k);
    h.CData = real(log10(ZP)');
    title(['Year ' num2str(k)])
    gif
end

%% D
ZD = NaN*ones(ni,nj);
ZD(GRD.ID) = allD(:,1);

figure(6)
h = imagescn(geolat_t',geolon_t',log10(ZD)');
cb = colorbar;
ylabel(cb,'log10 abund (g m^-^2)')
cmocean dense
title('Forage fish')
%title(datestr(datenum(0,0,1),'dd'))
clim([-2 2])

hold on
% he = earthimage;
% uistack(he,'bottom')

gif([ppath exper 'D.gif'],'frame',gcf,'delaytime',1/12,'nodither')

for k=2:50
    ZD = NaN*ones(ni,nj);
    ZD(GRD.ID) = allD(:,k);

    yd = yr(k);
    h.CData = real(log10(ZD)');
    title(['Year ' num2str(k)])
    gif
end

%% S
ZS = NaN*ones(ni,nj);
ZS(GRD.ID) = allS(:,1);

figure(7)
h = imagescn(geolat_t',geolon_t',log10(ZS)');
cb = colorbar;
ylabel(cb,'log10 abund (g m^-^2)')
cmocean dense
title('Forage fish')
%title(datestr(datenum(0,0,1),'dd'))
clim([-2 2])

hold on
% he = earthimage;
% uistack(he,'bottom')

gif([ppath exper 'S.gif'],'frame',gcf,'delaytime',1/12,'nodither')

for k=2:50
    ZS = NaN*ones(ni,nj);
    ZS(GRD.ID) = allS(:,k);

    yd = yr(k);
    h.CData = real(log10(ZS)');
    title(['Year ' num2str(k)])
    gif
end

%% M
ZM = NaN*ones(ni,nj);
ZM(GRD.ID) = allM(:,1);

figure(8)
h = imagescn(geolat_t',geolon_t',log10(ZM)');
cb = colorbar;
ylabel(cb,'log10 abund (g m^-^2)')
cmocean dense
title('Forage fish')
%title(datestr(datenum(0,0,1),'dd'))
clim([-2 2])

hold on
% he = earthimage;
% uistack(he,'bottom')

gif([ppath exper 'M.gif'],'frame',gcf,'delaytime',1/12,'nodither')

for k=2:50
    ZM = NaN*ones(ni,nj);
    ZM(GRD.ID) = allM(:,k);

    yd = yr(k);
    h.CData = real(log10(ZM)');
    title(['Year ' num2str(k)])
    gif
end

%% L
ZL = NaN*ones(ni,nj);
ZL(GRD.ID) = allL(:,1);

figure(9)
h = imagescn(geolat_t',geolon_t',log10(ZL)');
cb = colorbar;
ylabel(cb,'log10 abund (g m^-^2)')
cmocean dense
title('Forage fish')
%title(datestr(datenum(0,0,1),'dd'))
clim([-2 2])

hold on
% he = earthimage;
% uistack(he,'bottom')

gif([ppath exper 'L.gif'],'frame',gcf,'delaytime',1/12,'nodither')

for k=2:50
    ZL = NaN*ones(ni,nj);
    ZL(GRD.ID) = allL(:,k);

    yd = yr(k);
    h.CData = real(log10(ZL)');
    title(['Year ' num2str(k)])
    gif
end

%% ALL
ZAll = NaN*ones(ni,nj);
ZAll(GRD.ID) = allFish(:,1);

figure(10)
h = imagescn(geolat_t',geolon_t',log10(ZAll)');
cb = colorbar;
ylabel(cb,'log10 abund (g m^-^2)')
cmocean dense
title('Forage fish')
%title(datestr(datenum(0,0,1),'dd'))
clim([-2 2])

hold on
% he = earthimage;
% uistack(he,'bottom')

gif([ppath exper 'ALL.gif'],'frame',gcf,'delaytime',1/12,'nodither')

for k=2:50
    ZAll = NaN*ones(ni,nj);
    ZAll(GRD.ID) = allFish(:,k);

    yd = yr(k);
    h.CData = real(log10(ZAll)');
    title(['Year ' num2str(k)])
    gif
end

