% Use interpolated time series from GFDL to calc annual means
% Historic CMIP6 scenario 1880-2014
% Use MZ and LZ instead of one mesozoo
% For comparing to Kathryn's

clear
close all

%%
% fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/GFDL/hist/';
% gpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/GFDL/';
fpath='/project/Feisty/Fish-MIP/CMIP6/GFDL/hist/';
gpath='/project/Feisty/Fish-MIP/CMIP6/GFDL/';

load([gpath 'Data_grid_gfdl.mat'],'GRD');
load([gpath 'gridspec_gfdl_cmip6.mat']);

%%
yrs = 1950:2014;
nyrs = length(yrs);

tp_amean = nan*ones(GRD.N,nyrs);
tb_amean = tp_amean;
mz_amean = tp_amean;
lz_amean = tp_amean;
mzhp_amean = tp_amean;
lzhp_amean = tp_amean;
det_amean = tp_amean;

%%
for y = 1:nyrs
    YR = yrs(y)

    load([fpath 'Data_gfdl_cmip6_hist_2meso_daily_',num2str(YR),'.mat'],'ESM');

    %%
    tp_amean(:,y) = squeeze(mean(ESM.Tp,2,'omitnan'));
    tb_amean(:,y) = squeeze(mean(ESM.Tb,2,'omitnan'));
    mz_amean(:,y) = squeeze(mean(ESM.Zm,2,'omitnan'));
    lz_amean(:,y) = squeeze(mean(ESM.Zl,2,'omitnan'));
    mzhp_amean(:,y) = squeeze(mean(ESM.dZm,2,'omitnan'));
    lzhp_amean(:,y) = squeeze(mean(ESM.dZl,2,'omitnan'));
    det_amean(:,y) = squeeze(mean(ESM.det,2,'omitnan'));

end

%%
%Time
tp_tmean = squeeze(mean(tp_amean,1,'omitnan'));
tb_tmean = squeeze(mean(tb_amean,1,'omitnan'));
mz_tmean = squeeze(mean(mz_amean,1,'omitnan'));
lz_tmean = squeeze(mean(lz_amean,1,'omitnan'));
mzhp_tmean = squeeze(mean(mzhp_amean,1,'omitnan'));
lzhp_tmean = squeeze(mean(lzhp_amean,1,'omitnan'));
det_tmean = squeeze(mean(det_amean,1,'omitnan'));

%Space
tp_smean = squeeze(mean(tp_amean,2,'omitnan'));
tb_smean = squeeze(mean(tb_amean,2,'omitnan'));
mz_smean = squeeze(mean(mz_amean,2,'omitnan'));
lz_smean = squeeze(mean(lz_amean,2,'omitnan'));
mzhp_smean = squeeze(mean(mzhp_amean,2,'omitnan'));
lzhp_smean = squeeze(mean(lzhp_amean,2,'omitnan'));
det_smean = squeeze(mean(det_amean,2,'omitnan'));

%%
save([fpath 'Time_Means_gfdl_cmip6_hist_2meso_daily.mat'],'yrs',...
    'tb_tmean','tp_tmean','mz_tmean',...
    'mzhp_tmean','lz_tmean','lzhp_tmean',...
    'det_tmean')

save([fpath 'Space_Means_gfdl_cmip6_hist_2meso_daily.mat'],'yrs',...
    'tb_smean','tp_smean','mz_smean',...
    'mzhp_smean','lz_smean','lzhp_smean',...
    'det_smean')

save([fpath 'Annual_Means_gfdl_cmip6_hist_2meso_daily.mat'],'yrs',...
    'tb_amean','tp_amean','mz_amean',...
    'mzhp_amean','lz_amean','lzhp_amean',...
    'det_amean')

