% Visualize advection test cases

clear all
close all

%dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/advect_tests/';
dpath = '/Volumes/GFDL/CSV/advect_tests/';
%dpath = '/Volumes/GFDL/NC/AdvectTests/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/advect_tests/';

bio = csvread([dpath 'bio_2Ddiff_test_Atl_dt1hr_v22.csv']);
cname = 'diff_test_Atl_dt1hr_v1';

grid = csvread('grid_csv.csv');
load('gridspec_forecast.mat');

[nd,nid] = size(bio);

% Conservation of mass
area = grid(:,5);
area = area';
mass = bio .* repmat(area,nd,1);
totb = sum(mass,2);
figure(10)
subplot(2,2,1)
plot(totb)
cons = 100*(totb(end)-totb(1))/totb(1)

%
yrs=[1:length(totb)]/365;
figure(11)
plot(yrs,totb,'LineWidth',2)
%xlim([0 12])
%xlabel('days')
xlabel('Year')
%title(cname)
ylabel('Total number of particles')
print('-dpng',[fpath 'advec_test_' cname '_totb.png'])

%% plot info
% Land
surf_tmask = tmask(:,:,1);
lmask = surf_tmask;
lmask(lmask==0) = 999;
lmask(lmask==1) = NaN;

%axes labels
xt = -250:50:50;
xl = xt;
xl(xl<-180) = xl(xl<-180) + 350;


%% Global flat
%t = 1:72.75:nd;
t = 1:3:15;
t = round(t);

%Global map
for n=1:length(t)
    B1 = NaN*ones(size(geolat_t));
    B1(grid(:,1))=bio(t(n),:);
    
    figure
    surf(geolon_t,geolat_t,B1);
    view(2);
    shading flat;
    colorbar;
    caxis([0 1e5]);
    colormap('jet')
    title(['Day ' num2str(t(n)) ' Year 1'])
    %print('-dpng',[fpath 'advec_test_' cname '_' num2str(t(n)) '.png'])
    print('-dpng',[fpath cname '_' num2str(t(n)) '.png'])
end


%% Arctic projection
for n=1:length(t)
    B1(grid(:,1))=bio(t(n),:);
    
    figure
    m_proj('stereographic','lat',90,'long',30,'radius',30);
    m_pcolor(geolon_t,geolat_t,B1);
    shading flat
    colorbar
    colormap('jet')
    caxis([0 1e5])
    m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    title(['Day ' num2str(t(n)) ' Year 1'])
    %print('-dpng',[fpath 'advec_test_' cname '_arcticproj_' num2str(t(n)) '.png'])
    print('-dpng',[fpath cname '_arcticproj_' num2str(t(n)) '.png'])
end

%% Antarctic projection
for n=1:length(t)
    B1(grid(:,1))=bio(t(n),:);
    
    figure
    m_proj('stereographic','lat',-90,'long',30,'radius',50);
    m_pcolor(geolon_t,geolat_t,B1);
    shading interp
    colorbar
    colormap('jet')
    caxis([0 2e6])
    m_grid('xtick',12,'tickdir','out','ytick',[-50 -60 -70],'linest','-');
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    title(['Day ' num2str(t(n)) ' Year 1'])
    print('-dpng',[fpath 'advec_test_' cname '_Spoleproj_' num2str(t(n)) '.png'])
end