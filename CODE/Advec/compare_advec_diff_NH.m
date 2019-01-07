% Compare N Atl advection versions

clear all
close all

fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/advect_tests/';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'AREA_OCN','dat','dxtn','dyte','ht','geolon_t','geolat_t');

cname1='AtlNH_dt1hr_velH50_Jan88_b100';
cname2='AtlNH_dt1hr_velH50_Jan88_b100_v2';

%% end distributions
load(['/Volumes/GFDL/CSV/advect_tests/POEM_adv_diff_test_' cname1 '.mat'],'TF2');
TF = TF2;
clear TF2
load(['/Volumes/GFDL/CSV/advect_tests/POEM_adv_diff_test_' cname2 '.mat'],'TF2');

area = AREA_OCN;
area = (area)*510072000*1e6;
area = max(area,1);

%Natasha's region
[ni,nj] = size(geolon_t);
iids = [205:296];
jids = [149:178];
TF1 = zeros(ni,nj);
TF1(iids,jids) = TF;

land = find(ht<=0);
TF1(land) = nan;
TF2(land) = nan;

%% Plot
figure(1)
surf(geolon_t,geolat_t,TF1); view(2); shading interp;
colorbar
caxis([0 100]); 
axis([min(geolon_t(:)) max(geolon_t(:)) min(geolat_t(:)) max(geolat_t(:))])
title('Small grid 1988 day 365');
print('-dpng',[fpath 'POEM_adv_diff_test_tracer_' cname1 '_day365.png'])

figure(2)
surf(geolon_t,geolat_t,TF2); view(2); shading interp;
colorbar
caxis([0 100]); 
axis([min(geolon_t(:)) max(geolon_t(:)) min(geolat_t(:)) max(geolat_t(:))])
title('Whole grid 1988 day 365');
print('-dpng',[fpath 'POEM_adv_diff_test_tracer_' cname2 '_day365.png'])

%%
iids2 = [206:295];
jids2 = [150:177];

TFdiff=TF1-TF2;
TFdiff2=TFdiff(iids2,jids2);
TFdiff(land) = nan;

TFdiff3=100*(TF1-TF2)./TF1;
TFdiff3=TFdiff3(iids2,jids2);

figure(3)
pcolor(TFdiff2')
colorbar
caxis([-400 400]);
colormap('jet')
shading flat
title('Raw difference (Small grid - Whole grid)')
set(gca,'XTickLabel',geolon_t(iids2,1),'YTickLabel',geolat_t(jids2,1))
print('-dpng',[fpath 'POEM_adv_diff_test_tracer_' cname2 '_Diffv1-v2.png'])

figure(4)
pcolor(TFdiff3')
colorbar
caxis([-100 100]);
colormap('jet')
shading flat
title('Percent difference (Small grid - Whole grid)/Small grid')
set(gca,'XTickLabel',geolon_t(iids2,1),'YTickLabel',geolat_t(jids2,1))
print('-dpng',[fpath 'POEM_adv_diff_test_tracer_' cname2 '_Pdiffv1-v2.png'])

figure(5)
surf(geolon_t,geolat_t,TFdiff); view(2); shading interp;
colorbar
caxis([-300 300]); 
title('Small - Whole');
