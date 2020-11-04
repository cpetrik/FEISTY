% CMIP6 IPSL output
% Redid 100m integration

clear all
close all

%% Paths

ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/preindust/';
hpath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
spath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ssp126/';
rpath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ssp585/';
ppath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP6/';

load('/Volumes/MIP/Fish-MIP/CMIP6/IPSL/Data_grid_ipsl.mat','GRD');
load('/Volumes/MIP/Fish-MIP/CMIP6/IPSL/gridspec_ipsl_cmip6.mat');

CGRD = GRD;
clear GRD

%% Old Units
%zoo: mol C m-3

% meso zoo: from molC m-3 to g(WW) m-2
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
% mult by 10 m depth interval for m-3 to m-2
% D_Zm(j,:) = yi * 12.01 * 9.0 * 10;

% Hist
load([hpath 'ipsl_hist_zmeso100_monthly_1950_2014.mat'],'zmeso_100','runs');

zmeso_100(zmeso_100 > 1.0e19) = nan;
zoo = double(zmeso_100) * 12.01 * 9.0 * 10;

%
t=1:length(runs);
mo=t/12;
mo=mo+1950;
hmo = mo;

%map
zoo_hist1 = nanmean(zoo,3);
%ts
zoo_ts1 = squeeze(nanmean(nanmean(zoo,1),2));

clear zmeso_100 zoo 

%% New Units
%zoo: mol C m-2

% meso zoo: from molC m-3 to g(WW) m-2
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
% D_Zm(j,:) = yi * 12.01 * 9.0;

% Hist
load([hpath 'ipsl_hist_zmeso_100_monthly_1950_2014.mat'],'zmeso_100');

zmeso_100(zmeso_100 > 1.0e19) = nan;
zoo = double(zmeso_100) * 12.01 * 9.0;
%map
zoo_hist2 = nanmean(zoo,3);
%ts
zoo_ts2 = squeeze(nanmean(nanmean(zoo,1),2));

clear zmeso_100 zoo 

%%
clatlim=[-90 90];
clonlim=[-180 180];

%% Hist
figure(1)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(zoo_hist1))
cmocean('tempo')
caxis([0 2])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('IPSL Hist mesoz old')
% load coast;
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(zoo_hist2))
cmocean('tempo')
caxis([0 2])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('mesoz new')

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(zoo_hist1 - zoo_hist2))
cmocean('balance')
caxis([-20 20])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('mesoz old - new')

%subplot('Position',[0.5 0 0.5 0.5])

print('-dpng',[ppath 'Map_IPSL_Hist_diff_zmeso100.png'])

%% Timeseries
yr = (t/12)+1950;
figure(2)
plot(yr,zoo_ts1,'r'); hold on;
plot(yr,zoo_ts2,'b'); hold on;
legend('old','new')
title('IPSL')
ylabel('mean zoo (gWW m^-^2)')
print('-dpng',[ppath 'ts_IPSL_Hist_diff_zmeso100.png'])

%% Quantiles

q_tp(1,:) = quantile(tp_pi1(:),[0.05 0.25 0.5 0.75 0.95]);
q_tp(2,:) = quantile(tp_pi2(:),[0.05 0.25 0.5 0.75 0.95]);
q_tp(3,:) = quantile(tp_hist(:),[0.05 0.25 0.5 0.75 0.95]);
q_tp(4,:) = quantile(tp_126(:),[0.05 0.25 0.5 0.75 0.95]);
q_tp(5,:) = quantile(tp_585(:),[0.05 0.25 0.5 0.75 0.95]);

q_tb(1,:) = quantile(tb_pi1(:),[0.05 0.25 0.5 0.75 0.95]);
q_tb(2,:) = quantile(tb_pi2(:),[0.05 0.25 0.5 0.75 0.95]);
q_tb(3,:) = quantile(tb_hist(:),[0.05 0.25 0.5 0.75 0.95]);
q_tb(4,:) = quantile(tb_126(:),[0.05 0.25 0.5 0.75 0.95]);
q_tb(5,:) = quantile(tb_585(:),[0.05 0.25 0.5 0.75 0.95]);

q_zoo(1,:) = quantile(zoo_pi1(:),[0.05 0.25 0.5 0.75 0.95]);
q_zoo(2,:) = quantile(zoo_pi2(:),[0.05 0.25 0.5 0.75 0.95]);
q_zoo(3,:) = quantile(zoo_hist(:),[0.05 0.25 0.5 0.75 0.95]);
q_zoo(4,:) = quantile(zoo_126(:),[0.05 0.25 0.5 0.75 0.95]);
q_zoo(5,:) = quantile(zoo_585(:),[0.05 0.25 0.5 0.75 0.95]);

q_det(1,:) = quantile(det_pi1(:),[0.05 0.25 0.5 0.75 0.95]);
q_det(2,:) = quantile(det_pi2(:),[0.05 0.25 0.5 0.75 0.95]);
q_det(3,:) = quantile(det_hist(:),[0.05 0.25 0.5 0.75 0.95]);
q_det(4,:) = quantile(det_126(:),[0.05 0.25 0.5 0.75 0.95]);
q_det(5,:) = quantile(det_585(:),[0.05 0.25 0.5 0.75 0.95]);

% save('/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ipsl_cm4_input_means.mat',...
%     'q_tp','q_tb','q_zoo','q_det','-append')
% 
% TP = array2table(q_tp,'RowNames',{'pi1','pi2','hist','ssp126','ssp585'},...
%     'VariableNames',{'5th','25th','50th','75th','95th'});
% writetable(TP,'/Volumes/MIP/Fish-MIP/CMIP6/IPSL/Tp_quantiles.csv','Delimiter',',','WriteRowNames',true)
% 
% TB = array2table(q_tb,'RowNames',{'pi1','pi2','hist','ssp126','ssp585'},...
%     'VariableNames',{'5th','25th','50th','75th','95th'});
% writetable(TB,'/Volumes/MIP/Fish-MIP/CMIP6/IPSL/Tb_quantiles.csv','Delimiter',',','WriteRowNames',true)
% 
% TZ = array2table(q_zoo,'RowNames',{'pi1','pi2','hist','ssp126','ssp585'},...
%     'VariableNames',{'5th','25th','50th','75th','95th'});
% writetable(TZ,'/Volumes/MIP/Fish-MIP/CMIP6/IPSL/Zoo_quantiles.csv','Delimiter',',','WriteRowNames',true)
% 
% TD = array2table(q_det,'RowNames',{'pi1','pi2','hist','ssp126','ssp585'},...
%     'VariableNames',{'5th','25th','50th','75th','95th'});
% writetable(TD,'/Volumes/MIP/Fish-MIP/CMIP6/IPSL/Det_quantiles.csv','Delimiter',',','WriteRowNames',true)
% 

