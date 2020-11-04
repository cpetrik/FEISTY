% CMIP6 IPSL output

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

%% Units
%poc flux: mmol C m-2 s-1
%zoo: mol C m-2
%tp: degC
%tb: degC

% meso zoo: from molC m-2 to g(WW) m-2
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
%D_Zm(j,:) = yi * 12.01 * 9.0;

% detrital flux to benthos: from molC m-2 s-1 to g(WW) m-2 d-1
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
% 60*60*24 sec in a day
%D_det(j,:) = yi * 12.01 * 9.0 * 60 * 60 * 24;

%% Preindust
load([ipath 'ipsl_pi_temp_100_monthly_1950_2100.mat'],'temp_100');
load([ipath 'ipsl_pi_temp_btm_monthly_1950_2100.mat'],'temp_btm');
load([ipath 'ipsl_pi_zmeso_100_monthly_1950_2100.mat'],'zmeso_100');
load([ipath 'ipsl_pi_det_btm_monthly_1950_2100.mat']); %,'det_btm'

temp_100(temp_100 > 1.0e19) = nan;
temp_btm(temp_btm > 1.0e19) = nan;
zmeso_100(zmeso_100 > 1.0e19) = nan;
det_btm(det_btm > 1.0e19) = nan;

% flip vint and vavg
temp_100 = fliplr(temp_100);
zmeso_100 = fliplr(zmeso_100);

tp = double(temp_100);
tb = double(temp_btm);
zoo = double(zmeso_100) * 12.01 * 9.0;
det = double(det_btm) * 12.01 * 9.0 * 60 * 60 * 24;

%
t=time+1;
mo=t/12;
mo=mo+1600;

pmo = mo(runs);

yr1=find(pmo>1990 & pmo<=2000); 
yr2=find(pmo>2090 & pmo<=2100);

tp_pi1 = nanmean(tp(:,:,yr1),3);
tb_pi1 = nanmean(tb(:,:,yr1),3);
zoo_pi1 = nanmean(zoo(:,:,yr1),3);
det_pi1 = nanmean(det(:,:,yr1),3);

tp_pi2 = nanmean(tp(:,:,yr2),3);
tb_pi2 = nanmean(tb(:,:,yr2),3);
zoo_pi2 = nanmean(zoo(:,:,yr2),3);
det_pi2 = nanmean(det(:,:,yr2),3);

clear temp_100 temp_btm zmeso_100 det_btm runs
clear tp tb zoo det

%% Hist
load([hpath 'ipsl_hist_temp_100_monthly_1950_2014.mat'],'temp_100');
load([hpath 'ipsl_hist_temp_btm_monthly_1950_2014.mat'],'temp_btm');
load([hpath 'ipsl_hist_zmeso_100_monthly_1950_2014.mat'],'zmeso_100');
load([hpath 'ipsl_hist_det_btm_monthly_1950_2014.mat']);

temp_100(temp_100 > 1.0e19) = nan;
temp_btm(temp_btm > 1.0e19) = nan;
zmeso_100(zmeso_100 > 1.0e19) = nan;
det_btm(det_btm > 1.0e19) = nan;

% flip vint and vavg
temp_100 = fliplr(temp_100);

tp = double(temp_100);
tb = double(temp_btm);
zoo = double(zmeso_100) * 12.01 * 9.0;
det = double(det_btm) * 12.01 * 9.0 * 60 * 60 * 24;

%
t=1:length(runs);
mo=t/12;
mo=mo+1950;

hmo = mo;

yr1=find(hmo>1990 & hmo<=2000); 
yr2=find(hmo>2090 & hmo<=2100);

tp_hist = nanmean(tp(:,:,yr1),3);
tb_hist = nanmean(tb(:,:,yr1),3);
zoo_hist = nanmean(zoo(:,:,yr1),3);
det_hist = nanmean(det(:,:,yr1),3);

clear temp_100 temp_btm zmeso_100 det_btm
clear tp tb zoo det

%% SSP 126
load([spath 'ipsl_ssp126_temp_100_monthly_2015_2100.mat'],'temp_100');
load([spath 'ipsl_ssp126_temp_btm_monthly_2015_2100.mat'],'temp_btm');
load([spath 'ipsl_ssp126_zmeso_100_monthly_2015_2100.mat'],'zmeso_100');
load([spath 'ipsl_ssp126_det_btm_monthly_2015_2100.mat']);

temp_100(temp_100 > 1.0e19) = nan;
temp_btm(temp_btm > 1.0e19) = nan;
zmeso_100(zmeso_100 > 1.0e19) = nan;
det_btm(det_btm > 1.0e19) = nan;

% flip vint and vavg
zmeso_100 = fliplr(zmeso_100);

tp = double(temp_100);
tb = double(temp_btm);
zoo = double(zmeso_100) * 12.01 * 9.0;
det = double(det_btm) * 12.01 * 9.0 * 60 * 60 * 24;

%
t=1:length(time);
mo=t/12;
mo=mo+2015;

smo = mo;

yr1=find(smo>1990 & smo<=2000); 
yr2=find(smo>2090 & smo<=2100);

tp_126 = nanmean(tp(:,:,yr2),3);
tb_126 = nanmean(tb(:,:,yr2),3);
zoo_126 = nanmean(zoo(:,:,yr2),3);
det_126 = nanmean(det(:,:,yr2),3);

clear temp_100 temp_btm zmeso_100 det_btm
clear tp tb zoo det

%% SSP 585
load([rpath 'ipsl_ssp585_temp_100_monthly_2015_2100.mat'],'temp_100');
load([rpath 'ipsl_ssp585_temp_btm_monthly_2015_2100.mat'],'temp_btm');
load([rpath 'ipsl_ssp585_zmeso_100_monthly_2015_2100.mat'],'zmeso_100');
load([rpath 'ipsl_ssp585_det_btm_monthly_2015_2100.mat']);

temp_100(temp_100 > 1.0e19) = nan;
temp_btm(temp_btm > 1.0e19) = nan;
zmeso_100(zmeso_100 > 1.0e19) = nan;
det_btm(det_btm > 1.0e19) = nan;

tp = double(temp_100);
tb = double(temp_btm);
zoo = double(zmeso_100) * 12.01 * 9.0;
det = double(det_btm) * 12.01 * 9.0 * 60 * 60 * 24;

%
t=1:length(time);
mo=t/12;
mo=mo+2015;

rmo = mo;

yr1=find(rmo>1990 & rmo<=2000); 
yr2=find(rmo>2090 & rmo<=2100);

tp_585 = nanmean(tp(:,:,yr2),3);
tb_585 = nanmean(tb(:,:,yr2),3);
zoo_585 = nanmean(zoo(:,:,yr2),3);
det_585 = nanmean(det(:,:,yr2),3);

clear temp_100 temp_btm zmeso_100 det_btm
clear tp tb zoo det

%%
save('/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ipsl_cm4_input_means.mat',...
    'tp_pi1','tb_pi1','zoo_pi1','det_pi1',...
    'tp_pi2','tb_pi2','zoo_pi2','det_pi2',...
    'tp_hist','tb_hist','zoo_hist','det_hist',...
    'tp_126','tb_126','zoo_126','det_126',...
    'tp_585','tb_585','zoo_585','det_585')

%%
clatlim=[-90 90];
clonlim=[-180 180];

%% PI Hist
figure(1)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tp_pi1)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('IPSL PI 1990-2000 Tp')
% load coast;
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tb_pi1)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('Tb')

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(zoo_pi1))
cmocean('tempo')
caxis([0 2])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('log_1_0 Zoo (g m^-^2)')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(det_pi1))
cmocean('tempo')
caxis([-2 1])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('log_1_0 Det (g m^-^2)')
print('-dpng',[ppath 'Map_IPSL_Pre_1990_2000.png'])

%% PI Fore
figure(2)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tp_pi2)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('IPSL PI 2090-2100 Tp')
% load coast;
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tb_pi2)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('Tb')

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(zoo_pi2))
cmocean('tempo')
caxis([0 2])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('log_1_0 Zoo (g m^-^2)')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(det_pi2))
cmocean('tempo')
caxis([-2 1])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('log_1_0 Det (g m^-^2)')
print('-dpng',[ppath 'Map_IPSL_Pre_2090_2100.png'])

%% Hist
figure(3)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tp_hist)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('IPSL Hist 1990-2000 Tp')
% load coast;
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tb_hist)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('Tb')

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(zoo_hist))
cmocean('tempo')
caxis([0 2])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('log_1_0 Zoo (g m^-^2)')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(det_hist))
cmocean('tempo')
caxis([-2 1])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('log_1_0 Det (g m^-^2)')
print('-dpng',[ppath 'Map_IPSL_Hist_1990_2000.png'])

%% SSP 126
figure(4)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tp_126)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('IPSL SSP126 2090-2100 Tp')
% load coast;
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tb_126)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('Tb')

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(zoo_126))
cmocean('tempo')
caxis([0 2])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('log_1_0 Zoo (g m^-^2)')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(det_126))
cmocean('tempo')
caxis([-2 1])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('log_1_0 Det (g m^-^2)')
print('-dpng',[ppath 'Map_IPSL_SSP126_2090_2100.png'])

%% SSP 585
figure(5)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tp_585)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('IPSL SSP585 2090-2100 Tp')
% load coast;
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tb_585)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('Tb')

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(zoo_585))
cmocean('tempo')
caxis([0 2])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('log_1_0 Zoo (g m^-^2)')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(det_585))
cmocean('tempo')
caxis([-2 1])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('log_1_0 Det (g m^-^2)')
print('-dpng',[ppath 'Map_IPSL_SSP585_2090_2100.png'])

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

save('/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ipsl_cm4_input_means.mat',...
    'q_tp','q_tb','q_zoo','q_det','-append')

TP = array2table(q_tp,'RowNames',{'pi1','pi2','hist','ssp126','ssp585'},...
    'VariableNames',{'5th','25th','50th','75th','95th'});
writetable(TP,'/Volumes/MIP/Fish-MIP/CMIP6/IPSL/Tp_quantiles.csv','Delimiter',',','WriteRowNames',true)

TB = array2table(q_tb,'RowNames',{'pi1','pi2','hist','ssp126','ssp585'},...
    'VariableNames',{'5th','25th','50th','75th','95th'});
writetable(TB,'/Volumes/MIP/Fish-MIP/CMIP6/IPSL/Tb_quantiles.csv','Delimiter',',','WriteRowNames',true)

TZ = array2table(q_zoo,'RowNames',{'pi1','pi2','hist','ssp126','ssp585'},...
    'VariableNames',{'5th','25th','50th','75th','95th'});
writetable(TZ,'/Volumes/MIP/Fish-MIP/CMIP6/IPSL/Zoo_quantiles.csv','Delimiter',',','WriteRowNames',true)

TD = array2table(q_det,'RowNames',{'pi1','pi2','hist','ssp126','ssp585'},...
    'VariableNames',{'5th','25th','50th','75th','95th'});
writetable(TD,'/Volumes/MIP/Fish-MIP/CMIP6/IPSL/Det_quantiles.csv','Delimiter',',','WriteRowNames',true)


