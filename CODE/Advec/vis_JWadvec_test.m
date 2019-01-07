% Visualize advection test cases

clear all
close all

dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/advect_tests/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/advect_tests/';

bio = csvread([dpath 'bio_JWadvect_test_Atl_noswim_1d.csv']);
cname = 'Atl_noswim_1d';

grid = csvread('grid_csv.csv');
load('gridspec_forecast.mat');

% Conservation of mass

totb = sum(bio,2);
%
figure(10)
subplot(2,2,1)
plot(totb)
100*(totb(1)-totb(end))/totb(1)
%           dt=1d       dt=6hr      dt=3hr      dt=1hr  
%Eq =                                           -4.6e+93
%Pac =                                          -1.3e+24
%Atl =      -5.7e+95%   -1.2e+106   -4.4e+93    -9.6e+83
%no reflect                         -2.7e+09
%2D noreflec                        -9.1e+07
%full2D nor                         -2.5e+03
%Arctic =                                       -7.5e+29
%Antarc =                                       -37.9%
%WInd =                      
%EInd =             


%
B1=NaN*ones(size(geolat_t));
B2=B1;
B3=B1;
B4=B1;
B5=B1;

%for n=1:size(bio,1);
    B1(grid(:,1))=bio(1,:);
    B2(grid(:,1))=bio(73,:); %73
    B3(grid(:,1))=bio(146,:); %146
    B4(grid(:,1))=bio(219,:); %219
    B5(grid(:,1))=bio(366,:); %365
%end

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


%% Start
figure(11)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=0')
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,B1); hold on;
shading flat
% Link them together
linkaxes([ax1,ax2])
% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
% Give each one its own colormap
colormap(ax1,'gray')
colormap(ax2,'jet')
caxis([-1 2])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_1.png'])

%
figure(12)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=73') %73
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,B2); hold on;
shading flat
% Link them together
linkaxes([ax1,ax2])
% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
% Give each one its own colormap
colormap(ax1,'gray')
colormap(ax2,'jet')
caxis([-1 2])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_2.png'])

figure(13)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=146') 
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,B3); hold on;
shading flat
% Link them together
linkaxes([ax1,ax2])
% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
% Give each one its own colormap
colormap(ax1,'gray')
colormap(ax2,'jet')
caxis([-1 2])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_3.png'])

figure(14)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=219') %219
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,B4); hold on;
shading flat
% Link them together
linkaxes([ax1,ax2])
% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
% Give each one its own colormap
colormap(ax1,'gray')
colormap(ax2,'jet')
caxis([-1 2])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_4.png'])

%
figure(15)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=365')
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,B5); hold on;
shading flat
% Link them together
linkaxes([ax1,ax2])
% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
% Give each one its own colormap
colormap(ax1,'gray')
colormap(ax2,'jet')
caxis([-1 2])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_5.png'])

%% Arctic projection
figure(16)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,B1);
shading flat
colorbar
colormap('jet')
caxis([0 10])
title('Day=1')
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_1_arcticproj.png'])

figure(17)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,B2);
shading flat
colorbar
colormap('jet')
caxis([0 10])
title('Day=73')
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_2_arcticproj.png'])

figure(18)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,B3);
shading flat
colorbar
colormap('jet')
caxis([0 10])
title('Day=146')
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_3_arcticproj.png'])

figure(19)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,B4);
shading flat
colorbar
colormap('jet')
caxis([0 10])
title('Day=219')
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_4_arcticproj.png'])

figure(20)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,B5);
shading flat
colorbar
colormap('jet')
caxis([0 10])
title('Day=365')
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_5_arcticproj.png'])

%% Antarctic projection
figure(26)
m_proj('stereographic','lat',-90,'long',30,'radius',50);
m_pcolor(geolon_t,geolat_t,B1);
shading flat
colorbar
colormap('jet')
caxis([0 10])
title('Day=1')
m_grid('xtick',12,'tickdir','out','ytick',[-50 -60 -70],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_1_Spoleproj.png'])

figure(27)
m_proj('stereographic','lat',-90,'long',30,'radius',50);
m_pcolor(geolon_t,geolat_t,B2);
shading flat
colorbar
colormap('jet')
caxis([0 10])
title('Day=73')
m_grid('xtick',12,'tickdir','out','ytick',[-50 -60 -70],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_2_Spoleproj.png'])

figure(28)
m_proj('stereographic','lat',-90,'long',30,'radius',50);
m_pcolor(geolon_t,geolat_t,B3);
shading flat
colorbar
colormap('jet')
caxis([0 10])
title('Day=146')
m_grid('xtick',12,'tickdir','out','ytick',[-50 -60 -70],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_3_Spoleproj.png'])

figure(29)
m_proj('stereographic','lat',-90,'long',30,'radius',50);
m_pcolor(geolon_t,geolat_t,B4);
shading flat
colorbar
colormap('jet')
caxis([0 10])
title('Day=219')
m_grid('xtick',12,'tickdir','out','ytick',[-50 -60 -70],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_4_Spoleproj.png'])

figure(30)
m_proj('stereographic','lat',-90,'long',30,'radius',50);
m_pcolor(geolon_t,geolat_t,B5);
shading flat
colorbar
colormap('jet')
caxis([0 10])
title('Day=365')
m_grid('xtick',12,'tickdir','out','ytick',[-50 -60 -70],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_5_Spoleproj.png'])

