% Get mean biomasses of temp, det, zoop

clear 
close all

%%
fpath='/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/IPSLdown/';

load([fpath 'feisty_ipsl_gridspec.mat'],'TMO','TMO_bnds','tmo_units','BATHY')
load([fpath 'Data_grid_nemuro_ipsl.mat'],'GRD');

%%
time = (TMO/365)+1900;
yrs= 1980:2100;
Tdays=1:365;

mos = length(time);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;

%% Units
%poc_bot: mmolN/m3/d %Might be mmolN/m2 now?
%lzoo: mmolN/m2
%pzoo: mmolN/m2
%tp: degC
%tb: degC

load([fpath 'feisty_ipsl_lzoo_1980-2100.mat'],'LZOO_INT_200M',...
    'lz_long_name','lz_units');

load([fpath 'feisty_ipsl_pon_1980-2100.mat'],'PON_BOT',...
    'pon_long_name','pon_units');

load([fpath 'feisty_ipsl_pzoo_1980-2100.mat'],'PZOO_INT_200M',...
    'pz_long_name','pz_units');

load([fpath 'feisty_ipsl_tb_1980-2100.mat'],'TEMP_BOT','tb_long_name','tb_units');

load([fpath 'feisty_ipsl_tp_1980-2100.mat'],'TEMP_AVG_200M',...
    'tp_long_name','tp_units','LAT','LON');

[ni,nj,nt] = size(TEMP_BOT);

%%
LZOO_INT_200M(LZOO_INT_200M<0) = 0.0;
PZOO_INT_200M(PZOO_INT_200M<0) = 0.0;
PON_BOT(PON_BOT<0) = 0.0;

%%
figure
pcolor(PON_BOT(:,:,1)); shading flat;
colorbar
%caxis([0 2])

%%
WID = GRD.ID;  %1560 of these NaNs - maybe FW bodies on land?
NID = GRD.N;

Tp = reshape(TEMP_AVG_200M,ni*nj,nt);
Tb = reshape(TEMP_BOT,ni*nj,nt);
Zm = reshape(LZOO_INT_200M,ni*nj,nt);
Zl = reshape(PZOO_INT_200M,ni*nj,nt);
det= reshape(PON_BOT,ni*nj,nt);

Tp = Tp(WID,:);
Tb = Tb(WID,:);
Zm = Zm(WID,:);
Zl = Zl(WID,:);
det= det(WID,:);

%% Units
Zm = Zm * 1e-3 * (106.0/16.0) * 12.01 * 9.0;
Zl = Zl * 1e-3 * (106.0/16.0) * 12.01 * 9.0;
det = det * 1e-3  * (106.0/16.0) * 12.01 * 9.0 * 20;

%% annual means
D_Tp  = nan*ones(NID,nyrs);
D_Tb  = nan*ones(NID,nyrs);
D_Zm  = nan*ones(NID,nyrs);
D_Zl  = nan*ones(NID,nyrs);
D_det = nan*ones(NID,nyrs);

for y = 1:nyrs
    D_Tp(:,y)  = mean(Tp(:,mstart(y):mend(y)),2,"omitnan");
    D_Tb(:,y)  = mean(Tb(:,mstart(y):mend(y)),2,"omitnan");
    D_Zm(:,y)  = mean(Zm(:,mstart(y):mend(y)),2,"omitnan");
    D_Zl(:,y)  = mean(Zl(:,mstart(y):mend(y)),2,"omitnan");
    D_det(:,y) = mean(det(:,mstart(y):mend(y)),2,"omitnan");
    
end

%%

ptemp_tmean_ipsl = double(mean(D_Tp,1,"omitnan"));
btemp_tmean_ipsl = double(mean(D_Tb,1,"omitnan"));
det_tmean_ipsl = double(mean(D_det,1,"omitnan"));
mz_tmean_ipsl = double(mean(D_Zm,1,"omitnan"));
lz_tmean_ipsl = double(mean(D_Zl,1,"omitnan"));

thist = find(yrs>=1980 & yrs<2001);
tssp = find(yrs>=2080 & yrs<2101);

ptemp_hist_mean_ipsl = double(mean(D_Tp,2,"omitnan"));
btemp_hist_mean_ipsl = double(mean(D_Tb,2,"omitnan"));
det_hist_mean_ipsl = double(mean(D_det,2,"omitnan"));
mz_hist_mean_ipsl = double(mean(D_Zm,2,"omitnan"));
lz_hist_mean_ipsl = double(mean(D_Zl,2,"omitnan"));

ptemp_ssp_mean_ipsl = double(mean(D_Tp,2,"omitnan"));
btemp_ssp_mean_ipsl = double(mean(D_Tb,2,"omitnan"));
det_ssp_mean_ipsl = double(mean(D_det,2,"omitnan"));
mz_ssp_mean_ipsl = double(mean(D_Zm,2,"omitnan"));
lz_ssp_mean_ipsl = double(mean(D_Zl,2,"omitnan"));

%%
save([fpath 'nemuro_ipsl_means.mat'],'WID',...
    'D_det','D_Tp','D_Tb','D_Zm','D_Zl',...,
    'ptemp_hist_mean_ipsl','btemp_hist_mean_ipsl','det_hist_mean_ipsl',...
    'mz_hist_mean_ipsl','lz_hist_mean_ipsl',...
    'ptemp_ssp_mean_ipsl','btemp_ssp_mean_ipsl','det_ssp_mean_ipsl',...
    'mz_ssp_mean_ipsl','lz_ssp_mean_ipsl',...
    'ptemp_tmean_ipsl','btemp_tmean_ipsl','det_tmean_ipsl',...
    'mz_tmean_ipsl','lz_tmean_ipsl');






