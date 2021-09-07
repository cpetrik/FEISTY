% map interan var of inputs from CORE 

clear all
close all

fpath = '/Volumes/MIP/GCM_DATA/CORE-forced/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% grid info
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t','AREA_OCN');
load([cpath 'lme_mask_esm2m.mat']);

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);
tlme = lme_mask_esm2m';
LAT = double(geolat_t);
LON = double(geolon_t);

%% Forcing
%load([fpath 'cobalt_core_variance_1950_2007.mat'])
load([fpath 'ocean_cobalt_temp100_monthly_1950_2007.mat'],'tp_100');
load([fpath 'ocean_cobalt_temp_btm_monthly_1950_2007.mat'],'tb');
load([fpath 'ocean_cobalt_mz100_monthly_1950_2007.mat'],'mz_100');
load([fpath 'ocean_cobalt_lz100_monthly_1950_2007.mat'],'lz_100');
load([fpath 'ocean_cobalt_hploss_mz100_monthly_1950_2007.mat'],'hploss_mz_100');
load([fpath 'ocean_cobalt_hploss_lz100_monthly_1950_2007.mat'],'hploss_lz_100');
load([fpath 'ocean_cobalt_fndet_btm_monthly_1950_2007.mat']); %,'det_btm'

%% change units
tp_100 = double(tp_100);
tbtm = double(tb);
mz_100 = double(mz_100) * (106.0/16.0) * 12.01 * 9.0;
lz_100 = double(lz_100) * (106.0/16.0) * 12.01 * 9.0;
det_btm = double(det_btm) * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
hploss_mz_100 = double(hploss_mz_100) * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
hploss_lz_100 = double(hploss_lz_100) * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;

z_100 = mz_100 + lz_100;
zloss_100 = hploss_mz_100 + hploss_lz_100;

%% annual means
[ni,nj,nt] = size(tp_100);
st=1:12:nt;
en=12:12:nt;
nyr = nt/12;

tp = nan*ones(ni,nj,nyr);
tb = nan*ones(ni,nj,nyr);
det = nan*ones(ni,nj,nyr);
zoo = nan*ones(ni,nj,nyr);
zlos = nan*ones(ni,nj,nyr);
for n=1:nyr
    tp(:,:,n)=nanmean(tp_100(:,:,st(n):en(n)),3);
    tb(:,:,n)=nanmean(tbtm(:,:,st(n):en(n)),3);
    det(:,:,n)=nanmean(det_btm(:,:,st(n):en(n)),3);
    zoo(:,:,n)=nanmean(z_100(:,:,st(n):en(n)),3);
    zlos(:,:,n)=nanmean(zloss_100(:,:,st(n):en(n)),3);
end

%% anomalies
atp = tp - nanmean(tp,3);
atb = tb - nanmean(tb,3);
adet = det - nanmean(det,3);
azoo = zoo - nanmean(zoo,3);
azlos = zlos - nanmean(zlos,3);

save([fpath 'CORE_interann_mean_forcings_anom.mat'],...
    'tp','tb','det','zoo','zlos',...
    'atp','atb','adet','azoo','azlos');

%% var by grid cell
vtp = var(atp,0,3,'omitnan');
vtb = var(atb,0,3,'omitnan');
vdet = var(adet,0,3,'omitnan');
vzoo = var(azoo,0,3,'omitnan');
vzlos = var(azlos,0,3,'omitnan');

%% var by lme
atp2 = reshape(atp,ni*nj,nyr);
atb2 = reshape(atb,ni*nj,nyr);
adet2 = reshape(adet,ni*nj,nyr);
amz = reshape(azoo,ni*nj,nyr);
amzl = reshape(azlos,ni*nj,nyr);

lme_tp_var1 = NaN*ones(66,1);
lme_tb_var1 = NaN*ones(66,1);
lme_det_var1 = NaN*ones(66,1);
lme_mz_var1 = NaN*ones(66,1);
lme_mzl_var1 = NaN*ones(66,1);
lme_tp_var2 = NaN*ones(66,1);
lme_tb_var2 = NaN*ones(66,1);
lme_det_var2 = NaN*ones(66,1);
lme_mz_var2 = NaN*ones(66,1);
lme_mzl_var2 = NaN*ones(66,1);

for L=1:66
    lid = find(tlme==L);
    
    ltp = atp2(lid,:);
    ltb = atb2(lid,:);
    ldet = adet2(lid,:);
    lmz = amz(lid,:);
    lmzl = amzl(lid,:);
    
    lme_tp_var1(L,1) = var(ltp(:),'omitnan');
    lme_tb_var1(L,1) = var(ltb(:),'omitnan');
    lme_det_var1(L,1) = var(ldet(:),'omitnan');
    lme_mz_var1(L,1) = var(lmz(:),'omitnan');
    lme_mzl_var1(L,1) = var(lmzl(:),'omitnan');
    
    lme_tp_var2(L,1) = nanmean(var(ltp,0,2));
    lme_tb_var2(L,1) = nanmean(var(ltb,0,2));
    lme_det_var2(L,1) = nanmean(var(ldet,0,2));
    lme_mz_var2(L,1) = nanmean(var(lmz,0,2));
    lme_mzl_var2(L,1) = nanmean(var(lmzl,0,2));
    
end
% var1 == var2

%% map info
clatlim=[-90 90];
clonlim=[-280 80];
load coastlines

cmYOR=cbrewer('seq','YlOrRd',66,'pchip');
cmYOB=cbrewer('seq','YlOrBr',66,'pchip');
cmOR=cbrewer('seq','OrRd',66,'pchip');

%% map
figure(1)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,vtp)
colormap(cmYOR)
caxis([0 1])
colorbar%('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'var Tp','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,vtb)
colormap(cmYOR)
caxis([0 0.3])
colorbar%('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'var Tb','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,vzoo)
colormap(cmYOR)
caxis([0 1.5])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'var Zoo','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,vdet)
colormap(cmYOR)
caxis([0 0.02])
colorbar%('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'var Det','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,vzlos)
colormap(cmYOR)
caxis([0 0.02])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'var Zoo loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[pp 'Map_CORE_interann_var_forcings.png'])

%% save
save([fpath 'CORE_interann_var_forcings.mat'],...
    'vtp','vtb','vdet','vzoo','vzlos',...
    'lme_tp_var1','lme_tb_var1','lme_det_var1','lme_mz_var1','lme_mzl_var1');



