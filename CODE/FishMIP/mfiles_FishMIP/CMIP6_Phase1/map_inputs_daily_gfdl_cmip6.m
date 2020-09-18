% CMIP6 COBALT output 

clear all
close all

%% Paths

cpath = '/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/preindust/';
ppath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP6/';

%% CM4 -----------------------------------------------------------------

% Annual mean of daily values from interp
load('/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/Data_grid_gfdl.mat','GRD');
load('/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/gridspec_gfdl_cmip6.mat');
load([cpath 'Data_gfdl_spinup_daily_1949.mat'])

CGRD = GRD;
clear GRD

%%
c_Tp = mean(ESM.Tp,2);
c_Tb = mean(ESM.Tb,2);
c_D = mean(ESM.det,2);
c_Z = mean(ESM.Zm,2);

[mi,mj]=size(LON);
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

figure(2)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,cTp)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('GFDL Tp')
% load coast;                     
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,cTb)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('GFDL Tb')

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,cZ)
cmocean('tempo')
caxis([0 20])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('GFDL Zoo')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,cD)
cmocean('tempo')
caxis([0 3])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('GFDL Det')
print('-dpng',[ppath 'Map_GFDL_Pre_1949_from_daily_interp_forcings.png'])

% WORKS NOW!!!

