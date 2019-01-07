% Visualize advection test cases

clear all
close all

%dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/advect_tests/';
dpath = '/Volumes/GFDL/CSV/advect_tests/';
%dpath = '/Volumes/GFDL/NC/AdvectTests/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/advect_tests/';

bio = csvread([dpath 'bio_2Dadvect_swim_deep_test_Atl_vel0_dt1hr_j2_nodiv_divdepth3_passQ_depdiv0_v2_sep.csv']);
cname = 'swim_deep_test_Atl_vel0_dt1hr_j2_nodiv_divdepth3_passQ_depdiv0_v2_sep';

grid = csvread('grid_csv.csv');
load('gridspec_forecast.mat');

% Conservation of mass
area = grid(:,5);
area = area';
mass = bio .* repmat(area,365,1);
totb = sum(mass,2);
figure(10)
subplot(2,2,1)
plot(totb)
cons = 100*(totb(end)-totb(1))/totb(1)

%%
yrs=[1:length(totb)]/365;
figure(11)
plot(yrs,totb,'LineWidth',2)
%xlim([0 12])
%xlabel('days')
xlabel('Year')
%title(cname)
ylabel('Total number of particles')
print('-dpng',[fpath 'advec_test_' cname '_totb.png'])

%%
B1=NaN*ones(size(geolat_t));
B2=B1;
B3=B1;
B4=B1;
B5=B1;
B6=B1;
B7=B1;
B8=B1;
B9=B1;

B1(grid(:,1))=mass(1,:);
B2(grid(:,1))=mass(73,:); %73
B3(grid(:,1))=mass(146,:); %146
B4(grid(:,1))=mass(219,:); %219
B5(grid(:,1))=mass(365,:); %365
% B6(grid(:,1))=mass(456,:); %456
% B7(grid(:,1))=mass(547,:); %547
% B8(grid(:,1))=mass(638,:); %638
% B9(grid(:,1))=mass(730,:); %730

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

%%
if n == 1
    figure(1)
    surf(geolon_t,geolat_t,TF); view(2); shading interp; caxis([0 1]);
    %pause
end

%% Start
figure(1)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=1 Year 1')
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
caxis([0 5e16])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_1.png'])

figure(2)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=73 Year 1') %73
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
caxis([0 5e16])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_2.png'])

figure(3)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=146 Year 1')
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
caxis([0 5e16])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_3.png'])

figure(4)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=219 Year 1') %219
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
caxis([0 5e16])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_4.png'])

%
figure(5)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=365 Year 1')
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
caxis([0 5e16])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_5.png'])

%% 
% figure(6)
% % Create two axes
% ax1 = axes;
% pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
% shading flat
% title('Day=73 Year 1') %456
% view(2)
% ax2 = axes;
% pcolor(ax2,geolon_t,geolat_t,B6); hold on;
% shading flat
% % Link them together
% linkaxes([ax1,ax2])
% % Hide the top axes
% ax2.Visible = 'off';
% ax2.XTick = [];
% ax2.YTick = [];
% % Give each one its own colormap
% colormap(ax1,'gray')
% colormap(ax2,'jet')
% caxis([0 5e16])
% % Then add colorbars and get everything lined up
% set([ax1,ax2],'Position',[.1 .11 .75 .815]);
% cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
% set(ax1,'XTick',xt,'XTickLabel',xl)
% print('-dpng',[fpath 'advec_test_' cname '_6.png'])
% 
% figure(7)
% % Create two axes
% ax1 = axes;
% pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
% shading flat
% title('Day=146 Year 1') %547
% view(2)
% ax2 = axes;
% pcolor(ax2,geolon_t,geolat_t,B7); hold on;
% shading flat
% % Link them together
% linkaxes([ax1,ax2])
% % Hide the top axes
% ax2.Visible = 'off';
% ax2.XTick = [];
% ax2.YTick = [];
% % Give each one its own colormap
% colormap(ax1,'gray')
% colormap(ax2,'jet')
% caxis([0 5e16])
% % Then add colorbars and get everything lined up
% set([ax1,ax2],'Position',[.1 .11 .75 .815]);
% cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
% set(ax1,'XTick',xt,'XTickLabel',xl)
% print('-dpng',[fpath 'advec_test_' cname '_7.png'])
% 
% figure(8)
% % Create two axes
% ax1 = axes;
% pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
% shading flat
% title('Day=219 Year 2') %638
% view(2)
% ax2 = axes;
% pcolor(ax2,geolon_t,geolat_t,B8); hold on;
% shading flat
% % Link them together
% linkaxes([ax1,ax2])
% % Hide the top axes
% ax2.Visible = 'off';
% ax2.XTick = [];
% ax2.YTick = [];
% % Give each one its own colormap
% colormap(ax1,'gray')
% colormap(ax2,'jet')
% caxis([0 5e16])
% % Then add colorbars and get everything lined up
% set([ax1,ax2],'Position',[.1 .11 .75 .815]);
% cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
% set(ax1,'XTick',xt,'XTickLabel',xl)
% print('-dpng',[fpath 'advec_test_' cname '_8.png'])
% 
% figure(9)
% % Create two axes
% ax1 = axes;
% pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
% shading flat
% title('Day=365 Year 1') %730
% view(2)
% ax2 = axes;
% pcolor(ax2,geolon_t,geolat_t,B9); hold on;
% shading flat
% % Link them together
% linkaxes([ax1,ax2])
% % Hide the top axes
% ax2.Visible = 'off';
% ax2.XTick = [];
% ax2.YTick = [];
% % Give each one its own colormap
% colormap(ax1,'gray')
% colormap(ax2,'jet')
% caxis([0 5e16])
% % Then add colorbars and get everything lined up
% set([ax1,ax2],'Position',[.1 .11 .75 .815]);
% cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
% set(ax1,'XTick',xt,'XTickLabel',xl)
% print('-dpng',[fpath 'advec_test_' cname '_9.png'])

%% Arctic projection
figure(16)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,B1);
shading flat
colorbar
colormap('jet')
caxis([0 5e16])
title('Day=1')
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_arcticproj_1.png'])

figure(17)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,B2);
shading flat
colorbar
colormap('jet')
caxis([0 5e16])
title('Day=73')
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_arcticproj_2.png'])

figure(18)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,B3);
shading flat
colorbar
colormap('jet')
caxis([0 5e16])
title('Day=146')
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_arcticproj_3.png'])

figure(19)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,B4);
shading flat
colorbar
colormap('jet')
caxis([0 5e16])
title('Day=219')
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_arcticproj_4.png'])

figure(20)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,B5);
shading flat
colorbar
colormap('jet')
caxis([0 5e16])
title('Day=365')
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_arcticproj_5.png'])

%% Antarctic projection
figure(26)
m_proj('stereographic','lat',-90,'long',30,'radius',50);
m_pcolor(geolon_t,geolat_t,B1);
shading flat
colorbar
colormap('jet')
caxis([0 5e16])
title('Day=1')
m_grid('xtick',12,'tickdir','out','ytick',[-50 -60 -70],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_Spoleproj_1.png'])

figure(27)
m_proj('stereographic','lat',-90,'long',30,'radius',50);
m_pcolor(geolon_t,geolat_t,B2);
shading flat
colorbar
colormap('jet')
caxis([0 5e16])
title('Day=73')
m_grid('xtick',12,'tickdir','out','ytick',[-50 -60 -70],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_Spoleproj_2.png'])

figure(28)
m_proj('stereographic','lat',-90,'long',30,'radius',50);
m_pcolor(geolon_t,geolat_t,B3);
shading flat
colorbar
colormap('jet')
caxis([0 5e16])
title('Day=146')
m_grid('xtick',12,'tickdir','out','ytick',[-50 -60 -70],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_Spoleproj_3.png'])

figure(29)
m_proj('stereographic','lat',-90,'long',30,'radius',50);
m_pcolor(geolon_t,geolat_t,B4);
shading flat
colorbar
colormap('jet')
caxis([0 5e16])
title('Day=219')
m_grid('xtick',12,'tickdir','out','ytick',[-50 -60 -70],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_Spoleproj_4.png'])

figure(30)
m_proj('stereographic','lat',-90,'long',30,'radius',50);
m_pcolor(geolon_t,geolat_t,B5);
shading flat
colorbar
colormap('jet')
caxis([0 5e16])
title('Day=365')
m_grid('xtick',12,'tickdir','out','ytick',[-50 -60 -70],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_Spoleproj_5.png'])

