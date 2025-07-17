% Visualize advection test cases

clear all
close all

dpath = '/Volumes/GFDL/CSV/advect_tests/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/advect_tests/';

%biov = csvread([dpath 'POEM_adv_diff_Global_even_dt1hr_vel_daily_b100_core2000.csv']);
cname = 'Global_even_dt1hr_esm2m2000_vel_b100_area';
load([dpath 'Matlab_adv_diff_' cname '.mat']);

grid = csvread('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/grid_csv.csv');
load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/gridspec_forecast.mat');
load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/Data_hindcast_grid_cp2D.mat')


%% Conservation of mass
[nid,nd] = size(biov);

% Grid with area instead of vectors of area
[ni,nj] = size(GRD.area);

bio2 = NaN*ones(ni,nj,nd);
for d = 1:nd
    bio = NaN*ones(ni,nj);
    bio(grid(:,1)) = (biov(:,d));
    bio2(:,:,d) = bio;
end

mass = bio2 .* repmat(GRD.area,1,1,nd);
totb = squeeze(nansum(nansum(mass,1)));
figure(10)
subplot(2,2,1)
plot(totb)
cons = 100*(totb(end)-totb(1))/totb(1) %gives a different result from all the others

%%
yrs=[1:length(totb)]/365;
figure(11)
plot(yrs,totb,'LineWidth',2)
%xlim([0 12])
%xlabel('days')
xlabel('Year')
%title(cname)
ylabel('Total number of particles')
print('-dpng',[fpath 'advec_diff_test_' cname '_totb.png'])

%% plot info
% Land
lmask = GRD.mask;
lmask(lmask==0) = 999;
lmask(lmask==1) = NaN;

%axes labels
xt = -250:50:50;
xl = xt;
xl(xl<-180) = xl(xl<-180) + 350;

t = 1:72.75:nd;
%t = 1:3:15;
t = round(t);

%% Global flat
for n=1:length(t)
%     B1 = NaN*ones(size(geolat_t));
%     B1(grid(:,1))=bio(:,t(n));
    B1 = bio2(:,:,t(n));
    
    figure
    surf(geolon_t,geolat_t,B1);
    view(2);
    shading interp;
    colorbar;
    caxis([0 1e2]);
    colormap('jet')
    title(['Day ' num2str(t(n)) ' Year 1'])
    print('-dpng',[fpath 'advec_diff_test_' cname '_' num2str(t(n)) '.png'])
end


%% Arctic projection
for n=1:length(t)
%     B1(grid(:,1))=bio(:,t(n));
    B1 = bio2(:,:,t(n));
    
    figure
    m_proj('stereographic','lat',90,'long',30,'radius',30);
    m_pcolor(geolon_t,geolat_t,B1);
    shading interp
    colorbar
    colormap('jet')
    caxis([0 1e2])
    m_grid('xtick',6,'tickdir','out','ytick',[70 80],'linest','-');
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    title(['Day ' num2str(t(n)) ' Year 1'])
    print('-dpng',[fpath 'advec_diff_test_' cname '_arcticproj_' num2str(t(n)) '.png'])
end

%% Antarctic projection
for n=1:length(t)
%     B1(grid(:,1))=bio(:,t(n));
    B1 = bio2(:,:,t(n));
    
    figure
    m_proj('stereographic','lat',-90,'long',30,'radius',50);
    m_pcolor(geolon_t,geolat_t,B1);
    shading interp
    colorbar
    colormap('jet')
    caxis([0 1e2])
    m_grid('xtick',12,'tickdir','out','ytick',[-50 -60 -70],'linest','-');
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    title(['Day ' num2str(t(n)) ' Year 1'])
    print('-dpng',[fpath 'advec_diff_test_' cname '_Spoleproj_' num2str(t(n)) '.png'])
end