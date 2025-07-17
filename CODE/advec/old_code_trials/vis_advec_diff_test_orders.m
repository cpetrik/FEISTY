% Visualize advection test cases w/diffusion
% Test 1st order upwind, 2nd order centered, 4th order centered 
% advection schemes from MOM5

clear all
close all

dpath = '/Volumes/GFDL/CSV/advect_tests/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/Advec-Diff/figs_tests/';

cname = 'AtlArc_even_dt100hr_esm2m2000_vel_allgrid';
load([dpath 'Matlab_adv_diff_1st_2nd_4th_K600_' cname '.mat'],'biov1','biov2','biov4');

grid = csvread('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/grid_csv.csv');
load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/gridspec_forecast.mat');
load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/Data_hindcast_grid_cp2D.mat')

%% Conservation of mass
% 1st order upwind
[nid,nd] = size(biov1);

% Grid with area instead of vectors of area
[ni,nj] = size(GRD.area);
bio1 = NaN*ones(ni,nj,nd);
for d = 1:nd
    bio = NaN*ones(ni,nj);
    bio(grid(:,1)) = (biov1(:,d));
    bio1(:,:,d) = bio;
end

mass1 = bio1 .* repmat(GRD.area,1,1,nd);
totb1 = squeeze(nansum(nansum(mass1,1)));
cons1 = 100*(totb1(end)-totb1(1))/totb1(1)

%
yrs=[1:length(totb1)]/365;
figure(51)
plot(yrs,totb1,'LineWidth',2)
%xlim([0 12])
%xlabel('days')
xlabel('Year')
%title(cname)
ylabel('Total number of particles 1st')
%print('-dpng',[fpath 'advec_diff_test_K600_' cname '_totb_1st.png'])
print('-dpng',[fpath 'diff_test_K600_' cname '_totb.png'])

%% 2nd order
[nid,nd] = size(biov2);

% Grid with area instead of vectors of area
[ni,nj] = size(GRD.area);
bio2 = NaN*ones(ni,nj,nd);
for d = 1:nd
    bio = NaN*ones(ni,nj);
    bio(grid(:,1)) = (biov2(:,d));
    bio2(:,:,d) = bio;
end

mass2 = bio2 .* repmat(GRD.area,1,1,nd);
totb2 = squeeze(nansum(nansum(mass2,1)));
cons2 = 100*(totb2(end)-totb2(1))/totb2(1)

%
yrs=[1:length(totb2)]/365;
figure(52)
plot(yrs,totb2,'LineWidth',2)
%xlim([0 12])
%xlabel('days')
xlabel('Year')
%title(cname)
ylabel('Total number of particles 2nd')
print('-dpng',[fpath 'advec_diff_test_' cname '_totb_2nd.png'])

%% 4th order
[nid,nd] = size(biov4);

% Grid with area instead of vectors of area
[ni,nj] = size(GRD.area);
bio4 = NaN*ones(ni,nj,nd);
for d = 1:nd
    bio = NaN*ones(ni,nj);
    bio(grid(:,1)) = (biov4(:,d));
    bio4(:,:,d) = bio;
end

mass4 = bio4 .* repmat(GRD.area,1,1,nd);
totb4 = squeeze(nansum(nansum(mass4,1)));
cons4 = 100*(totb4(end)-totb4(1))/totb4(1)

%
yrs=[1:length(totb4)]/365;
figure(54)
plot(yrs,totb4,'LineWidth',2)
%xlim([0 12])
%xlabel('days')
xlabel('Year')
%title(cname)
ylabel('Total number of particles 4th')
print('-dpng',[fpath 'advec_diff_test_' cname '_totb_4th.png'])

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

% Land
surf_tmask = tmask(:,:,1);
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
    B1 = bio1(:,:,t(n));
%     B2 = bio2(:,:,t(n));
%     B4 = bio4(:,:,t(n));
    
    figure(1)
    clf
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,B1)
    colormap(cmB)
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colorbar;
    caxis([0 150]);
    %colormap('jet')
    colormap(cmB)
    title(['Day ' num2str(t(n)) ' Year 1'])
    print('-dpng',[fpath 'advec_diff_test_1st_K600_' cname '_' num2str(t(n)) '.png'])
    
%     figure(2)
%     clf
%     axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
%     surfm(geolat_t,geolon_t,B2)
%     colormap(cmB)
%     load coast;                     %decent looking coastlines
%     h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%     colorbar;
%     caxis([0 150]);
%     %colormap('jet')
%     colormap(cmB)
%     title(['Day ' num2str(t(n)) ' Year 1'])
%     print('-dpng',[fpath 'advec_diff_test_2nd_' cname '_' num2str(t(n)) '.png'])
%     
%     figure(4)
%     clf
%     axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
%     surfm(geolat_t,geolon_t,B4)
%     colormap(cmB)
%     load coast;                     %decent looking coastlines
%     h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%     colorbar;
%     caxis([0 150]);
%     %colormap('jet')
%     colormap(cmB)
%     title(['Day ' num2str(t(n)) ' Year 1'])
%     print('-dpng',[fpath 'advec_diff_test_4th_' cname '_' num2str(t(n)) '.png'])
end


%% Arctic projection
for n=1:length(t)
    B1 = bio1(:,:,t(n));
%     B2 = bio2(:,:,t(n));
%     B4 = bio4(:,:,t(n));
    
    figure(1)
    clf
    m_proj('stereographic','lat',90,'long',30,'radius',30);
    m_pcolor(geolon_t,geolat_t,B1);
    shading interp
    colorbar
    %colormap('jet')
    colormap(cmB)
    caxis([0 150])
    m_grid('xtick',6,'tickdir','out','ytick',[70 80],'linest','-');
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    title(['Day ' num2str(t(n)) ' Year 1'])
    print('-dpng',[fpath 'advec_diff_test_1st_K600_' cname '_arcticproj_' num2str(t(n)) '.png'])
    
%     figure(2)
%     clf
%     m_proj('stereographic','lat',90,'long',30,'radius',30);
%     m_pcolor(geolon_t,geolat_t,B2);
%     shading interp
%     colorbar
%     %colormap('jet')
%     colormap(cmB)
%     caxis([0 150])
%     m_grid('xtick',6,'tickdir','out','ytick',[70 80],'linest','-');
%     m_coast('patch',[.7 .7 .7],'edgecolor','k');
%     title(['Day ' num2str(t(n)) ' Year 1'])
%     print('-dpng',[fpath 'advec_diff_test_2nd_' cname '_arcticproj_' num2str(t(n)) '.png'])
% 
%     figure(4)
%     clf
%     m_proj('stereographic','lat',90,'long',30,'radius',30);
%     m_pcolor(geolon_t,geolat_t,B4);
%     shading interp
%     colorbar
%     %colormap('jet')
%     colormap(cmB)
%     caxis([0 150])
%     m_grid('xtick',6,'tickdir','out','ytick',[70 80],'linest','-');
%     m_coast('patch',[.7 .7 .7],'edgecolor','k');
%     title(['Day ' num2str(t(n)) ' Year 1'])
%     print('-dpng',[fpath 'advec_diff_test_4th_' cname '_arcticproj_' num2str(t(n)) '.png'])

end

%% Antarctic projection
for n=1:length(t)
    B1 = bio2(:,:,t(n));
    
    figure
    m_proj('stereographic','lat',-90,'long',30,'radius',50);
    m_pcolor(geolon_t,geolat_t,B1);
    shading interp
    colorbar
    %colormap('jet')
    colormap(cmB)
    caxis([0 150])
    m_grid('xtick',12,'tickdir','out','ytick',[-50 -60 -70],'linest','-');
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    title(['Day ' num2str(t(n)) ' Year 1'])
    print('-dpng',[fpath 'advec_diff_test_' cname '_Spoleproj_' num2str(t(n)) '.png'])
end