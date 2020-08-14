% Compare daily interp zooplankton and detritus from
% CMIP6 COBALT output and ESM2M output
% to check unit conversion

clear all
close all

%% Paths

cpath = '/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/preindust/';
npath = '/Volumes/FEISTY/GCM_DATA/PreIndust/';
gpath = '/Volumes/FEISTY/POEM_JLD/esm2m_hist/';

ppath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP6/';

%% ESM2M ----------------------------------------------------------------
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/Data_grid_hindcast_NOTflipped.mat','GRD');
load([gpath 'Data_ESM2Mhist_1949.mat'])

EGRD = GRD;
clear GRD

%% Annual mean of daily values from interp
% WID NOT EQUAL TO GRD.ID!!!

e_Tp = mean(COBALT.Tp,2);
e_Tb = mean(COBALT.Tb,2);
e_D = mean(COBALT.det,2);
e_Z = mean(COBALT.Zm,2)+mean(COBALT.Zl,2);

[ni,nj]=size(geolon_t);
eTp=NaN*ones(ni,nj);
eTb=NaN*ones(ni,nj);
eD=NaN*ones(ni,nj);
eZ=NaN*ones(ni,nj);

eTp(EGRD.ID)=e_Tp;
eTb(EGRD.ID)=e_Tb;
eD(EGRD.ID) =e_D;
eZ(EGRD.ID) =e_Z;

%%
elatlim=[-90 90];
elonlim=[-280 80];

figure(1)
subplot(2,2,1)
axesm ('Robinson','MapLatLimit',elatlim,'MapLonLimit',elonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,eTp)
cmocean('thermal')
caxis([0 35])
title('ESM2M Tp')
% load coast;                     
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot(2,2,2)
axesm ('Robinson','MapLatLimit',elatlim,'MapLonLimit',elonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,eTb)
cmocean('thermal')
caxis([0 35])
title('ESM2M Tb')

subplot(2,2,3)
axesm ('Robinson','MapLatLimit',elatlim,'MapLonLimit',elonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,eZ)
cmocean('tempo')
caxis([0 20])
title('ESM2M Z')

subplot(2,2,4)
axesm ('Robinson','MapLatLimit',elatlim,'MapLonLimit',elonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,eD)
cmocean('tempo')
caxis([0 3])
title('ESM2M D')
stamp('')
print('-dpng',[ppath 'Map_ESM2M_Hist_1949_from_daily_interp_forcings.png'])

%% Annual mean of monthly values from netcdf
% ncdisp([npath 'ocean_cobalt_btm.186101-200512.fndet_btm.nc.nc'])
% 
% ncid = netcdf.open([npath 'ocean_cobalt_btm.186101-200512.fndet_btm.nc.nc'],'NC_NOWRITE');
% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% for i = 1:nvars
%     varname = netcdf.inqVar(ncid, i-1);
%     eval([ varname ' = netcdf.getVar(ncid,i-1);']);
%     eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
% end
% netcdf.close(ncid);

%'ocean.186101-200512.temp100.nc'
%'ocean_cobalt_btm.186101-200512.btm_temp.nc'
%'ocean_cobalt_biomass_100.186101-200512.nmdz_100.nc'
%'ocean_cobalt_biomass_100.186101-200512.nlgz_100.nc'

%% CM4 -----------------------------------------------------------------

%% Annual mean of daily values from interp
load('/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/Data_grid_gfdl.mat','GRD');
load('/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/gridspec_gfdl.mat');
load([cpath 'Data_gfdl_spinup_daily_1949.mat'])

CGRD = GRD;
clear GRD

c_Tp = mean(ESM.Tp,2);
c_Tb = mean(ESM.Tb,2);
c_D = mean(ESM.det,2);
c_Z = mean(ESM.Zm,2);

[mi,mj]=size(geolon);
cTp=NaN*ones(mi,mj);
cTb=NaN*ones(mi,mj);
cD=NaN*ones(mi,mj);
cZ=NaN*ones(mi,mj);

cTp(CGRD.ID)=c_Tp;
cTb(CGRD.ID)=c_Tb;
cD(CGRD.ID) =c_D;
cZ(CGRD.ID) =c_Z;

%%
clatlim=[-90 90];
clonlim=[-180 180];
geolat = double(geolat);

figure(2)
subplot(2,2,1)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,cTp)
cmocean('thermal')
caxis([0 35])
title('CM4 Tp')
% load coast;                     
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot(2,2,2)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,cTb)
cmocean('thermal')
caxis([0 35])
title('CM4 Tb')

subplot(2,2,3)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,cZ)
cmocean('tempo')
caxis([0 20])
title('CM4 Z')

subplot(2,2,4)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,cD)
cmocean('tempo')
caxis([0 3])
title('CM4 D')
print('-dpng',[ppath 'Map_CM4_Pre_1949_from_daily_interp_forcings.png'])

%% Annual mean of monthly values from matfiles

load([cpath 'gfdl_pi_temp100_monthly_1850_1949.mat']);
load([cpath 'gfdl_pi_zmeso100_monthly_1850_1949.mat']);
load([cpath 'gfdl_pi_temp_btm_monthly_1850_1949.mat']);
load([cpath 'gfdl_pi_det_btm_monthly_1850_1949.mat']);
load('/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/gridspec_gfdl.mat');

lyr = (length(y)-11):(length(y));

%%
cTp = double(mean(temp_100(:,:,lyr),3));
cTb = double(mean(temp_btm(:,:,lyr),3));
cZ = double(mean(zmeso_100(:,:,lyr),3));
cD = double(mean(det_btm(:,:,lyr),3));

%%
clatlim=[-90 90];
clonlim=[-180 180];
geolat = double(geolat);

figure(3)
subplot(2,2,1)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,cTp)
cmocean('thermal')
caxis([0 35])
title('CM4 Tp')
% load coast;                     
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot(2,2,2)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,cTb)
cmocean('thermal')
caxis([0 35])
title('CM4 Tb')

subplot(2,2,3)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,cZ)
cmocean('tempo')
caxis([0 0.03])
title('CM4 Z')

subplot(2,2,4)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,cD)
cmocean('tempo')
caxis([0 5e-7])
title('CM4 D')
print('-dpng',[ppath 'Map_CM4_Pre_1949_from_monthly_intmat_forcings.png'])

%% Annual mean of monthly values from netcdf

