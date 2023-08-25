% CESM FOSI output 

clear all
close all

%% Paths

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load([fpath 'gridspec_POP_gx1v6.mat']);
load([fpath 'Data_grid_POP_gx1v6.mat'],'GRD');
load([fpath 'Data_cesm_fosi_daily_1.mat']);

%%
c_Tp = mean(ESM.Tp,2);
c_Tb = mean(ESM.Tb,2);
c_D = mean(ESM.det,2);
c_Z = mean(ESM.Zm,2);
c_Zl = mean(ESM.dZm,2);

[mi,mj]=size(TLONG);
cTp=NaN*ones(mi,mj);
cTb=NaN*ones(mi,mj);
cD=NaN*ones(mi,mj);
cZ=NaN*ones(mi,mj);
cZl=NaN*ones(mi,mj);

cTp(GRD.ID)=c_Tp;
cTb(GRD.ID)=c_Tb;
cD(GRD.ID) =c_D;
cZ(GRD.ID) =c_Z;
cZl(GRD.ID) =c_Zl;

%%
clatlim=[-90 90];
clonlim=[-180 180];

figure(2)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cTp)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('Tp')
% load coast;                     
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cTb)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('Tb')

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZ))
cmocean('tempo')
caxis([0 2])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('log_1_0 mesoZoo g m^-^2')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cD))
cmocean('tempo')
caxis([-0.5 1.5])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('log_1_0 Det g m^-^2 d^-^1')
print('-dpng',[pp 'Map_CESM_FOSI_yr1_from_daily_interp_forcings.png'])



