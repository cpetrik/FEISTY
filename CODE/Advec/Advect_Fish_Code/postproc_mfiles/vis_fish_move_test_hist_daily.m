% Visualize advection test cases

clear 
close all

%%
vpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CORE-forced/';

%1-D
load([vpath 'Data_grid_ocean_cobalt_ESM2Mcore.mat'],'GRD');
grid = GRD;
clear GRD

%2-D
load([vpath 'Data_hindcast_grid_cp2D.mat'],'GRD')

load([vpath 'ocean_cobalt_grid.mat'],'geolon_t','geolat_t');

%%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

spath = ['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/CORE/'];

cname='CORE_1988_pelfish_dt1d_velDY';

load([spath 'AdvectPred_' cname '.mat'],'MFbiov','LPbiov');

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/';
ppath = [pp cfile '/CORE/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%% Grid instead of vectors 
[ni,nj] = size(GRD.area);
[nid,nd] = size(MFbiov);

MFbio2 = NaN*ones(ni,nj,nd);
LPbio2 = NaN*ones(ni,nj,nd);

for d = 1:nd
    Fbio = NaN*ones(ni,nj);
    Fbio(grid.ID) = (MFbiov(:,d));
    MFbio2(:,:,d) = Fbio;

    Pbio = NaN*ones(ni,nj);
    Pbio(grid.ID) = (LPbiov(:,d));
    LPbio2(:,:,d) = Pbio;
end


%% plot info
% Land
surf_tmask = GRD.mask;
lmask = surf_tmask;
lmask(lmask==0) = 999;
lmask(lmask==1) = NaN;

%axes labels
xt = -250:50:50;
xl = xt;
xl(xl<-180) = xl(xl<-180) + 350;

t = 1:72.75:nd;
%t = 1:3:15;
t = round(t);

% colors
cmB=cbrewer('seq','Blues',50);

%% Global flat
for n=1:length(t)
    B1 = MFbio2(:,:,t(n));
    
    figure(1)
    surf(geolon_t,geolat_t,B1);
    view(2);
    shading interp;
    colorbar;
    clim([0 20]);
    %colormap('jet')
    colormap(cmB)
    title(['Day ' num2str(t(n)) ' Year 1'])
    print('-dpng',[ppath 'advec_test_MF_' cname '_' num2str(t(n)) '.png'])
end


%% Arctic projection
for n=1:length(t)
    B1 = MFbio2(:,:,t(n));
    
    figure(2)
    m_proj('stereographic','lat',90,'long',30,'radius',30);
    m_pcolor(geolon_t,geolat_t,B1);
    shading interp
    colorbar
    %colormap('jet')
    colormap(cmB)
    clim([0 20])
    m_grid('xtick',6,'tickdir','out','ytick',[70 80],'linest','-');
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    title(['Day ' num2str(t(n)) ' Year 1'])
    print('-dpng',[ppath 'advec_test_MF_' cname '_arcticproj_' num2str(t(n)) '.png'])
end

%% Global flat
for n=1:length(t)
    B2 = LPbio2(:,:,t(n));
    
    figure(1)
    surf(geolon_t,geolat_t,B2);
    view(2);
    shading interp;
    colorbar;
    clim([0 20]);
    %colormap('jet')
    colormap(cmB)
    title(['Day ' num2str(t(n)) ' Year 1'])
    print('-dpng',[ppath 'advec_test_LP_' cname '_' num2str(t(n)) '.png'])
end


%% Arctic projection
for n=1:length(t)
    B2 = LPbio2(:,:,t(n));
    
    figure(2)
    m_proj('stereographic','lat',90,'long',30,'radius',30);
    m_pcolor(geolon_t,geolat_t,B2);
    shading interp
    colorbar
    %colormap('jet')
    colormap(cmB)
    clim([0 20])
    m_grid('xtick',6,'tickdir','out','ytick',[70 80],'linest','-');
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    title(['Day ' num2str(t(n)) ' Year 1'])
    print('-dpng',[ppath 'advec_test_LP' cname '_arcticproj_' num2str(t(n)) '.png'])
end
