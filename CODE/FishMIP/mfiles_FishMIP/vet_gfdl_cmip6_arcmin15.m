% Save time chunks to calc change in future 
% GFDL mesoz prepared by ISIMIP team

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/';
hpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
ppath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/Phase3/quarter_degree_vetting/';

%% zmeso200, sst, schl
load([hpath 'gfdl_hist_zmeso_200_monthly_1950_2014.mat'])
load([hpath 'gfdl_hist_sst_monthly_1950_2014.mat'])
load([hpath 'gfdl_hist_surf_chl_monthly_1950_2014.mat'])

glat = lat;
glon = lon;
zmeso200_units = units_vint;
schl_units = units;

%% NPP
load('/Volumes/MIP/Fish-MIP/CMIP6/GFDL/gridspec_gfdl_cmip6.mat',...
    'LAT','LON');
load([hpath 'gfdl_hist_npp_monthly_1951_2014.mat'],'npp','units');

npp_units = units;
gLAT = LAT;
gLON = LON;

%% MLD
load([hpath 'gfdl_hist_mld_monthly_1965_2014.mat'],'mld','yr',...
    'long_name','standard_name','units');

mld_units = units;

%% intpoc, tbot, phyc, phydiat, phydiaz, zooc, zmeso, zmicro
load([hpath 'gfdl_hist_int_poc_monthly_1951_2014.mat']);
intpoc_units = units;

load([hpath 'gfdl_hist_tbot_monthly_1951_2014.mat']);
tbot_units = units;

load([hpath 'gfdl_hist_phyc_vint_monthly_1951_2014.mat']);
phyc_units = units;

load([hpath 'gfdl_hist_phydiat_vint_monthly_1951_2014.mat']);
diat_units = units;

load([hpath 'gfdl_hist_phydiaz_vint_monthly_1951_2014.mat']);
diaz_units = units;

load([hpath 'gfdl_hist_zooc_vint_monthly_1951_2014.mat']);
zooc_units = units;

load([hpath 'gfdl_hist_zmeso_vint_monthly_1951_2014.mat']);
zmeso_units = units;

load([hpath 'gfdl_hist_zmicro_vint_monthly_1951_2014.mat']);
zmicro_units = units;

%% saved
load([ppath 'vet_gfdl_hist_cmip6.mat']);

%% Nans
zmeso_200(zmeso_200>=1e19) = NaN;
sst(sst>=1e19) = NaN;
schl(schl>=1e19) = NaN;
npp(npp(:)>=1e20) = NaN;
int_poc(int_poc>=1e19) = NaN;
mld(mld>=1e19) = NaN;
phyc_vint(phyc_vint>=1e19) = NaN;
phydiat_vint(phydiat_vint>=1e19) = NaN;
phydiaz_vint(phydiaz_vint>=1e19) = NaN;
t_bot(t_bot>=1e19) = NaN;
zooc_vint(zooc_vint>=1e19) = NaN;
zmeso_vint(zmeso_vint>=1e19) = NaN;
zmicro_vint(zmicro_vint>=1e19) = NaN;

%% Doubles
zmeso_200=double(zmeso_200);
sst=double(sst);
schl=double(schl);
npp=double(npp);
int_poc=double(int_poc);
mld=double(mld);
phyc_vint=double(phyc_vint);
phydiat_vint=double(phydiat_vint);
phydiaz_vint=double(phydiaz_vint);
t_bot=double(t_bot);
zooc_vint=double(zooc_vint);
zmeso_vint=double(zmeso_vint);
zmicro_vint=double(zmicro_vint);

%% Convert chl, phyc, or zooc? - NO
%chl is kg/m3 in 1/4 deg
%npp is molC/m2/s in 1/4 deg
%phyc-vint is molC/m2 in 1/4 deg
%zooc-vint is molC/m2 in 1/4 deg

% to mgC/m2
%all models in molC: 12.01 g C in 1 mol C
%1e3 mg in 1 g
%gzmo_all = gzmo_all * 12.01 * 1e3;

%% space means
mzmeso_200=nanmean(zmeso_200,3);
msst=nanmean(sst,3);
mschl=nanmean(schl,3);
mnpp=nanmean(npp,3);
mint_poc=nanmean(int_poc,3);
mmld=nanmean(mld,3);
mphyc_vint=nanmean(phyc_vint,3);
mphydiat_vint=nanmean(phydiat_vint,3);
mphydiaz_vint=nanmean(phydiaz_vint,3);
mt_bot=nanmean(t_bot,3);
mzooc_vint=nanmean(zooc_vint,3);
mzmeso_vint=nanmean(zmeso_vint,3);
mzmicro_vint=nanmean(zmicro_vint,3);

%% Ranges
%quantiles
%qCAN(1,:) = quantile(cnpp(:),[0.05 0.25 0.5 0.75 0.95]);

qzmeso_200 = quantile(mzmeso_200(:),[0.0 1]);
qsst = quantile(msst(:),[0.0 1]);
qschl = quantile(mschl(:),[0.0 1]);
qnpp = quantile(mnpp(:),[0.0 1]);
qint_poc = quantile(mint_poc(:),[0.0 1]);
qmld = quantile(mmld(:),[0.0 1]);
qphyc_vint = quantile(mphyc_vint(:),[0.0 1]);
qphydiat_vint = quantile(mphydiat_vint(:),[0.0 1]);
qphydiaz_vint = quantile(mphydiaz_vint(:),[0.0 1]);
qt_bot = quantile(mt_bot(:),[0.0 1]);
qzooc_vint = quantile(mzooc_vint(:),[0.0 1]);
qzmeso_vint = quantile(mzmeso_vint(:),[0.0 1]);
qzmicro_vint = quantile(mzmicro_vint(:),[0.0 1]);

% Tq = array2table(qCAN,'VariableNames',{'5th','25th','50th','75th','95th'},...
%     'RowNames',{'CAN','CNRM','GFDL','IPSL','UKESM'});
% writetable(Tq,'npp_quantiles_rawunits.csv','WriteRowNames',true)

%%
clatlim=[-90 90];
clonlim=[-280 80];

mlonlim=[-180 180];

load coastlines;

%% MAPS
f5 = figure('Units','inches','Position',[1 3 7.5 10]);

%1 - Intpoc
subplot('Position',[0.01 0.8 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,mint_poc)
%cmocean('matter')
colormap('jet')
caxis([6e-5 0.007])
colorbar
text(-2.5,2.75,'int-poc','HorizontalAlignment','left')
text(-2.5,2.25,num2str(qint_poc),'HorizontalAlignment','left')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%2 - NPP
subplot('Position',[0.01 0.625 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,mnpp)
% cmocean('algae')
colormap('jet')
caxis([3.2e-9 1.2e-6])
colorbar
text(-2.5,2.75,'intpp','HorizontalAlignment','left')
text(-2.5,2.25,num2str(qnpp),'HorizontalAlignment','left')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%3 - phyc-vint
subplot('Position',[0.01 0.45 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,mphyc_vint)
% cmocean('algae')
colormap('jet')
caxis([0.0008 0.27])
colorbar
text(-2.5,2.75,'phyc-vint','HorizontalAlignment','left')
text(-2.5,2.25,num2str(qphyc_vint),'HorizontalAlignment','left')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%4 - phydiat-vint
subplot('Position',[0.01 0.275 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,mphydiat_vint)
% cmocean('algae')
colormap('jet')
caxis([0.0006 0.13])
colorbar
text(-2.5,2.75,'phydiat-vint','HorizontalAlignment','left')
text(-2.5,2.25,num2str(qphydiat_vint),'HorizontalAlignment','left')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%5 - phydiat-vint
% subplot('Position',[0.01 0.1 0.3 0.175])
% axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(gLAT,gLON,mphydiat_vint)
% cmocean('algae')
% colormap('jet')
% caxis([0.001 0.2])
% colorbar
% text(-1.75,1.75,['phydiat-vint ' num2str(qphydiat_vint)],'HorizontalAlignment','center')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% colorbar('Position',[0.0195 0.075 0.275 0.02],'orientation','horizontal');

% New col (middle)
%1 - phydiaz-vint
subplot('Position',[0.32 0.8 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,mphydiaz_vint)
% cmocean('algae')
% caxis([3e-14 0.0055])
colormap('jet')
caxis([2.5e-5 0.023])
colorbar
text(-2.5,2.75,'phydiaz-vint','HorizontalAlignment','left')
text(-2.5,2.25,num2str(qphydiaz_vint),'HorizontalAlignment','left')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%2 - Chl 
subplot('Position',[0.32 0.625 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,mschl)
% cmocean('algae')
colormap('jet')
caxis([4.5e-13 1.2e-6])
colorbar
text(-2.5,2.75,'schl','HorizontalAlignment','left')
text(-2.5,2.25,num2str(qschl),'HorizontalAlignment','left')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%3 - SST
subplot('Position',[0.32 0.45 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,msst)
% cmocean('thermal')
colormap('jet')
caxis([-1.8 30])
colorbar
text(-2.5,2.75,'tos','HorizontalAlignment','left')
text(-2.5,2.25,num2str(qsst),'HorizontalAlignment','left')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%4 - t bot
subplot('Position',[0.32 0.275 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,mt_bot)
% cmocean('thermal')
colormap('jet')
caxis([-2 31])
colorbar
text(-2.5,2.75,'tob','HorizontalAlignment','left')
text(-2.5,2.25,num2str(qt_bot),'HorizontalAlignment','left')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% %5 - zooc
% subplot('Position',[0.32 0.1 0.3 0.175])
% axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(gLAT,gLON,mzooc_vint)
% cmocean('tempo')
% colormap('jet')
% caxis([-100 100])
% text(-1.95,1.75,['zooc-vint ' num2str(qzooc_vint)],'HorizontalAlignment','center')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% 6 - zoo
subplot('Position',[0.63 0.8 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,mzooc_vint)
% cmocean('tempo')
% caxis([0.0055 0.065])
colormap('jet')
caxis([0.009 0.28])
colorbar
text(-2.5,2.75,'zooc-vint','HorizontalAlignment','left')
text(-2.5,2.25,num2str(qzooc_vint),'HorizontalAlignment','left')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%7 - zmeso
subplot('Position',[0.63 0.625 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,mzmeso_vint)
% cmocean('tempo')
% caxis([0.0055 0.065])
colormap('jet')
caxis([0.002 0.15])
colorbar
text(-2.5,2.75,'zmeso-vint','HorizontalAlignment','left')
text(-2.5,2.25,num2str(qzmeso_vint),'HorizontalAlignment','left')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%8 - zmicro
subplot('Position',[0.63 0.45 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,mzmicro_vint)
% cmocean('tempo')
% caxis([0.0055 0.065])
colormap('jet')
caxis([0.0065 0.125])
colorbar
text(-2.5,2.75,'zmicro-vint','HorizontalAlignment','left')
text(-2.5,2.25,num2str(qzmicro_vint),'HorizontalAlignment','left')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%9 - MLD
subplot('Position',[0.63 0.275 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,mmld)
% cmocean('deep')
colormap('jet')
caxis([0.5 300])
colorbar
text(-2.5,2.75,'MLD','HorizontalAlignment','left')
text(-2.5,2.25,num2str(qmld),'HorizontalAlignment','left')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% 
% %10 - Zoo uk
% subplot('Position',[0.63 0.1 0.3 0.175])
% axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(gLAT,gLON,pdiff_uz2)
% cmocean('tempo')
% colormap('jet')
% caxis([-100 100])
% text(-1.95,1.75,['int-poc ' num2str(qint_poc)],'HorizontalAlignment','center')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% colorbar('Position',[0.475 0.075 0.3 0.02],'orientation','horizontal');

print('-dpng',[ppath 'Map_GFDL_CMIP6_1deg_from_ISIMIP_v2_jet.png'])




