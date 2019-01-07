% Visualize advection test cases

clear all
close all

%dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/advect_tests/';
dpath = '/Volumes/GFDL/NC/AdvectTests/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/advect_tests/';

bio = csvread([dpath 'bio_2Dadvect_test_global_vel200_dt1hr_j2_30yrs.csv']);
cname = 'Global_vel200_dt1hr_j2_30yrs';

grid = csvread('grid_csv.csv');
load('gridspec_forecast.mat');

%% Conservation of mass

totb = sum(bio,2);
figure(10)
%subplot(2,2,1)
plot(totb)
100*(totb(1)-totb(end))/totb(1)

%%
yrs=[1:length(totb)]/365;
figure(11)
plot(yrs,totb,'LineWidth',2)
%xlim([0 12])
%xlabel('days')
xlabel('Year')
title('Global seed of 1e-3')
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

%for n=1:size(bio,1);
B1(grid(:,1))=bio(1,:);
B2(grid(:,1))=bio(73,:); %73
B3(grid(:,1))=bio(146,:); %146
B4(grid(:,1))=bio(219,:); %219
B5(grid(:,1))=bio(365,:); %365
B6(grid(:,1))=bio(10731,:); %456
B7(grid(:,1))=bio(10804,:); %547
B8(grid(:,1))=bio(10877,:); %638
B9(grid(:,1))=bio(10950,:); %730
%end

% plot info
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
caxis([0 2e3])
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
caxis([0 2e3])
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
caxis([0 2e3])
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
caxis([0 2e3])
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
caxis([0 2e3])
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
caxis([0 2e3])
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
caxis([0 2e3])
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
caxis([0 2e3])
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
caxis([0 2e3])
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
caxis([0 2e3])
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
caxis([0 2e3])
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
caxis([0 2e3])
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
caxis([0 2e3])
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
caxis([0 2e3])
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
caxis([0 2e3])
title('Day=365')
m_grid('xtick',12,'tickdir','out','ytick',[-50 -60 -70],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
print('-dpng',[fpath 'advec_test_' cname '_5_Spoleproj.png'])

%% Years 10 and 20

C1=NaN*ones(size(geolat_t));
C2=C1;
C3=C1;
C4=C1;
C5=C1;
C6=C1;
C7=C1;
C8=C1;
C9=C1;
C11=C1;
C12=C1;
C13=C1;
C14=C1;

C1(grid(:,1))=bio(3286,:); %yr20,d1
C2(grid(:,1))=bio(3377,:); %yr20,d92
C3(grid(:,1))=bio(3468,:); %yr20,d183
C4(grid(:,1))=bio(3559,:); %yr20,d274
C5(grid(:,1))=bio(3650,:); %yr20,d365
C6(grid(:,1))=bio(7027,:); %yr20,d92
C7(grid(:,1))=bio(7118,:); %yr20,d183
C8(grid(:,1))=bio(7209,:); %yr20,d274
C9(grid(:,1))=bio(7300,:); %yr20,d365

C11(grid(:,1))=bio(1552,:); %yr5,d92
C12(grid(:,1))=bio(1643,:); %yr5,d183
C13(grid(:,1))=bio(1734,:); %yr5,d274
C14(grid(:,1))=bio(1825,:); %yr5,d365

%% Year 5
figure(2)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=91 Year 5') 
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,C11); hold on;
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
caxis([0 2e3])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_6.png'])

figure(3)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=183 Year 5')
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,C12); hold on;
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
caxis([0 2e3])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_7.png'])

figure(4)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=274 Year 5') %219
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,C13); hold on;
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
caxis([0 2e3])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_8.png'])

%
figure(5)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=365 Year 5')
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,C14); hold on;
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
caxis([0 2e3])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_9.png'])

%% Year 10
figure(2)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=91 Year 10') 
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,C2); hold on;
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
caxis([0 2e3])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_10.png'])

figure(3)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=183 Year 10')
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,C3); hold on;
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
caxis([0 2e3])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_11.png'])

figure(4)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=274 Year 10') %219
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,C4); hold on;
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
caxis([0 2e3])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_12.png'])

%
figure(5)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=365 Year 10')
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,C5); hold on;
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
caxis([0 2e3])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_13.png'])

%% Year 20
figure(6)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=92 Year 20') %456
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,C6); hold on;
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
caxis([0 2e3])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_14.png'])

figure(7)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=183 Year 20') %547
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,C7); hold on;
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
caxis([0 2e3])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_15.png'])

figure(8)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=274 Year 20') %638
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,C8); hold on;
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
caxis([0 2e3])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_16.png'])

figure(9)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=365 Year 20') %730
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,C9); hold on;
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
caxis([0 2e3])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_17.png'])

%% Year 30
figure(6)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=92 Year 30') %456
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,B6); hold on;
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
caxis([0 2e3])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_18.png'])

figure(7)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=183 Year 30') %547
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,B7); hold on;
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
caxis([0 2e3])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_19.png'])

figure(8)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=274 Year 30') %638
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,B8); hold on;
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
caxis([0 2e3])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_20.png'])

figure(9)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=365 Year 30') %730
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,B9); hold on;
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
caxis([0 2e3])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_' cname '_21.png'])

