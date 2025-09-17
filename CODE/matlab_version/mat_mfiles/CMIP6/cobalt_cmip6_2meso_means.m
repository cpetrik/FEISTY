% CMIP6 COBALT outputs for FEISTY
% 2-mesozoo version

clear
close all

%%
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/GFDL/hist/';

% Units
%poc flux: mol C m-2 s-1
%zoo: mol N m-2 and mol N m-2 s-1
%tp: degC
%tb: degC

load([fpath 'gfdl_hist_temp_100_monthly_1950_2014.mat'],'temp_100');
load([fpath 'gfdl_hist_temp_btm_monthly_1950_2014.mat'],'temp_btm');
load([fpath 'gfdl_hist_det_btm_monthly_1950_2014.mat']);

load([fpath 'gfdl_hist_int100_2zmeso_reorient_monthly_1950_2014.mat']);

temp_100 = double(temp_100);
temp_btm = double(temp_btm);
det_btm = double(det_btm);

%%
gpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/GFDL/';
load([gpath 'Data_grid_gfdl.mat'],'GRD');
load([gpath 'gridspec_gfdl_cmip6.mat']);

%%
yr65 = find(year>=1950);
time = time(yr65);
yr = year(yr65);

mos = length(time);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;

yrs = floor(yr(1)):yr(end);

%% convert units to what FEISTY uses
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/FishMIP/FishMIP6/';

mdz_100 = mdz_100 * (106/16) * 12.01 * 9.0; 
lgz_100 = lgz_100 * (106/16) * 12.01 * 9.0; 
hp_mdz_100 = hp_mdz_100 * (106/16) * 12.01 * 9.0 * 60 * 60 * 24; 
hp_lgz_100 = hp_lgz_100 * (106/16) * 12.01 * 9.0 * 60 * 60 * 24; 
det_btm = det_btm * 12.01 * 9.0 * 60 * 60 * 24; 

%% Take means for visualization

%Time
tp_tmean = squeeze(mean(temp_100,1,'omitnan'));
tb_tmean = squeeze(mean(temp_btm,1,'omitnan'));
mz_tmean = squeeze(mean(mdz_100,1,'omitnan'));
lz_tmean = squeeze(mean(lgz_100,1,'omitnan'));
mzhp_tmean = squeeze(mean(hp_mdz_100,1,'omitnan'));
lzhp_tmean = squeeze(mean(hp_lgz_100,1,'omitnan'));
det_tmean = squeeze(mean(det_btm,1,'omitnan'));

tp_tmean = squeeze(mean(tp_tmean,1,'omitnan'));
tb_tmean = squeeze(mean(tb_tmean,1,'omitnan'));
mz_tmean = squeeze(mean(mz_tmean,1,'omitnan'));
lz_tmean = squeeze(mean(lz_tmean,1,'omitnan'));
mzhp_tmean = squeeze(mean(mzhp_tmean,1,'omitnan'));
lzhp_tmean = squeeze(mean(lzhp_tmean,1,'omitnan'));
det_tmean = squeeze(mean(det_tmean,1,'omitnan'));

%Space
tp_smean = squeeze(mean(temp_100,3,'omitnan'));
tb_smean = squeeze(mean(temp_btm,3,'omitnan'));
mz_smean = squeeze(mean(mdz_100,3,'omitnan'));
lz_smean = squeeze(mean(lgz_100,3,'omitnan'));
mzhp_smean = squeeze(mean(hp_mdz_100,3,'omitnan'));
lzhp_smean = squeeze(mean(hp_lgz_100,3,'omitnan'));
det_smean = squeeze(mean(det_btm,3,'omitnan'));

%% Annual means
nt = length(time);
nyr = nt/12;
% st=1:12:length(time);
% en=12:12:length(time);
st=1:12:nt;
en=12:12:nt;

[ni,nj] = size(LAT);

tp_amean = nan*ones(ni,nj,length(nyr));
tb_amean = tp_amean;
mz_amean = tp_amean;
lz_amean = tp_amean;
mzhp_amean = tp_amean;
lzhp_amean = tp_amean;
det_amean = tp_amean;

for n=1:length(st)
    % mean 
    tp_amean(:,:,n)=mean(temp_100(:,:,st(n):en(n)),3,'omitnan');
    tb_amean(:,:,n)=mean(temp_btm(:,:,st(n):en(n)),3,'omitnan');
    mz_amean(:,:,n)=mean(mdz_100(:,:,st(n):en(n)),3,'omitnan');
    lz_amean(:,:,n)=mean(lgz_100(:,:,st(n):en(n)),3,'omitnan');
    mzhp_amean(:,:,n)=mean(hp_mdz_100(:,:,st(n):en(n)),3,'omitnan');
    lzhp_amean(:,:,n)=mean(hp_lgz_100(:,:,st(n):en(n)),3,'omitnan');
    det_amean(:,:,n)=mean(det_btm(:,:,st(n):en(n)),3,'omitnan');
    
end

%%
save([fpath 'Time_Means_gfdl_cmip6_hist_2meso.mat'],'time',...
    'tb_tmean','tp_tmean','mz_tmean',...
    'mzhp_tmean','lz_tmean','lzhp_tmean',...
    'det_tmean')

save([fpath 'Space_Means_gfdl_cmip6_hist_2meso.mat'],'time',...
    'tb_smean','tp_smean','mz_smean',...
    'mzhp_smean','lz_smean','lzhp_smean',...
    'det_smean')

save([fpath 'Annual_Means_gfdl_cmip6_hist_2meso.mat'],'time',...
    'tb_amean','tp_amean','mz_amean',...
    'mzhp_amean','lz_amean','lzhp_amean',...
    'det_amean')
