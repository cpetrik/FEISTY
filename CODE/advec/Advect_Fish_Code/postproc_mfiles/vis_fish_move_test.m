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

%cname='Atl_even_dt1d_velMO_b100_swim10';
%cname='Atl_even_dt1d_velMO_b100_advectonly';
%cname='Atl_evenFish_randPrey_dt1d_velMO_b100';
%cname='Arctic_evenFish_randPrey_dt1d_velMO_b100_swim01';
%cname='Arctic_even_dt1d_velMO_b100_swim10';
%cname = 'Arctic_even_dt1d_velMO_b100_swim01';
%cname='Equat_even_dt1d_velMO_b100_swim10';
%cname='Global_evenfish_evenprey_dt1d_velDAY_b100_swim01';
cname='Global_ICs_MF_LZ_dt1d_velDAY_swim10';

load([spath 'AdvectPred_' cname '.mat'],'biov','preyv');

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/';
ppath = [pp cfile '/CORE/advect_test_figs/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%% Conservation of mass
[nid,nd] = size(biov);

% Grid with area instead of vectors of area
[ni,nj] = size(GRD.area);
bio2 = NaN*ones(ni,nj,nd);
for d = 1:nd
    bio = NaN*ones(ni,nj);
    bio(grid.ID) = (biov(:,d));
    bio2(:,:,d) = bio;
end

mass = bio2 .* repmat(GRD.area,1,1,nd);
totb = squeeze(nansum(nansum(mass,1)));
cons = 100*(totb(end)-totb(1))/totb(1)

massv = biov .* repmat(GRD.area(grid.ID),1,nd);
totbv = sum(massv,1,'omitnan');
cons2 = 100*(totbv(end)-totbv(1))/totbv(1)

%
yrs=[1:length(totb)]/365;
figure(11)
plot(yrs,totb,'LineWidth',2)
%xlim([0 12])
%xlabel('days')
xlabel('Year')
%title(cname)
ylabel('Total number of particles')
print('-dpng',[ppath 'advec_test_' cname '_totb.png'])

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
cmB=cbrewer('seq','Blues',50,'PCHIP');


%% Global flat
for n=1:length(t)
    B1 = bio2(:,:,t(n));
    
    figure
    surf(geolon_t,geolat_t,B1);
    view(2);
    shading interp;
    colorbar;
    clim([0 10]);
    %clim([0 20]);
    %clim([0 150]);
    %colormap('jet')
    colormap(cmB)
    title(['Day ' num2str(t(n)) ' Year 1'])
    print('-dpng',[ppath 'advec_test_' cname '_' num2str(t(n)) '.png'])
end


%% Arctic projection
% for n=1:length(t)
%     B1 = bio2(:,:,t(n));
% 
%     figure
%     m_proj('stereographic','lat',90,'long',30,'radius',30);
%     m_pcolor(geolon_t,geolat_t,B1);
%     shading interp
%     colorbar
%     %colormap('jet')
%     colormap(cmB)
%     clim([0 20]);
%     %clim([0 150]);
%     m_grid('xtick',6,'tickdir','out','ytick',[70 80],'linest','-');
%     m_coast('patch',[.7 .7 .7],'edgecolor','k');
%     title(['Day ' num2str(t(n)) ' Year 1'])
%     print('-dpng',[ppath 'advec_test_' cname '_arcticproj_' num2str(t(n)) '.png'])
% end
% 
% %% Antarctic projection
% for n=1:length(t)
%     B1 = bio2(:,:,t(n));
% 
%     figure
%     m_proj('stereographic','lat',-90,'long',30,'radius',50);
%     m_pcolor(geolon_t,geolat_t,B1);
%     shading interp
%     colorbar
%     %colormap('jet')
%     colormap(cmB)
%     clim([0 20]);
%     %clim([0 150]);
%     m_grid('xtick',12,'tickdir','out','ytick',[-50 -60 -70],'linest','-');
%     m_coast('patch',[.7 .7 .7],'edgecolor','k');
%     title(['Day ' num2str(t(n)) ' Year 1'])
%     print('-dpng',[ppath 'advec_test_' cname '_Spoleproj_' num2str(t(n)) '.png'])
% end