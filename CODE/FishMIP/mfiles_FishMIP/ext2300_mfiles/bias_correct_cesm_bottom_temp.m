% SSP 534-over theta-bot off from tob in SSP 585
% bias-correct with delta method using SSP585 tob 1st year
% same method as Daniele

clear 
close all

pp = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/WGs/2300/testing_forcing/';

%% 
spath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp585/';

load([spath 'cesm2_ssp585_temp_btm_monthly_2015_2299.mat']);
load('/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/gridspec_cesm2_cmip6_2300.mat')

temp_btm(temp_btm > 1.0e19) = nan;

ssp585_Tb = temp_btm;
ssp585_yr = yr;

clear temp_btm yr

%%
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp534over/';

load([fpath 'cesm2_ssp534-over_temp_btm_monthly_2040_2299.mat']);

temp_btm(temp_btm > 1.0e19) = nan;

ssp534_Tb = temp_btm;
ssp534_yr = yr;

clear temp_btm yr

%% find 2040
[iyr,yid] = intersect(ssp585_yr,ssp534_yr);

ssp585_yr(yid(1:12))
ssp534_yr(1:12)

%% bias
t2040 = mean(ssp585_Tb(:,:,yid(1:12)),3,'omitnan');
tob = mean(ssp534_Tb(:,:,1:12),3,'omitnan');
tdiff = tob - t2040;

%% map
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

% 3 figure subplot
figure(1)
subplot('Position',[0 0.53 0.5 0.5])
%585
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,t2040)
cmocean('thermal')
clim([0 10]);
set(gcf,'renderer','painters')
title('CESM 585 tob')

%534
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tob)
cmocean('thermal')
clim([0 10]);
set(gcf,'renderer','painters')
title('CESM 534 theta-bot')

%Diff
subplot('Position',[0.25 0.0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tdiff)
cmocean('balance')
clim([-10 10]);
colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('534 theta-bot - 585 tob')
stamp('')
print('-dpng',[pp 'CESM 2040 tob diff.png'])

%% Save
[ni,nj,nt] = size(ssp534_Tb);
tdiff_mat = repmat(tdiff,1,1,nt);
temp_btm = ssp534_Tb - tdiff_mat;

save([fpath 'cesm2_ssp534-over_temp_btm_corrected_monthly_2040_2299.mat'],'temp_btm',...
    'ssp534_yr','lat','lon','units');

%%
m585 = squeeze(mean(ssp585_Tb,1,'omitnan'));
m534 = squeeze(mean(ssp534_Tb,1,'omitnan'));
mtb = squeeze(mean(temp_btm,1,'omitnan'));

m585 = squeeze(mean(m585,1,'omitnan'));
m534 = squeeze(mean(m534,1,'omitnan'));
mtb = squeeze(mean(mtb,1,'omitnan'));

figure(2)
plot(ssp585_yr(1:12:end),m585(1:12:end),'r','LineWidth',2); hold on
plot(ssp534_yr(1:12:end),m534(1:12:end),'color',[0 0.75 0.5],'LineWidth',2);
plot(ssp534_yr(1:12:end),mtb(1:12:end),'color',[0 0 0.5],'LineWidth',2);
title(['CESM mean btm temp'])
legend('585','534','534 corrected')
legend('location','northwest')

print('-dpng',[pp 'CESM_bottom_temp_bias_corrected.png'])




