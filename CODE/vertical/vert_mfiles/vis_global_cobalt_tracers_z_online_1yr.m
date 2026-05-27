% Look at COBALT tracers from online sim
% 12 mo of output

clear
close all

%%
fpath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/cobalt_feisty/';

gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

cfile ='NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODe/Figs/Matlab_New_sizes/',...
    cfile,'/Cobalt_Feisty/'];

exper = 'OM4_05_COBALTv3_FEISTYon_021326';
Yr = '1994';
mod = [exper '_' Yr];

%%
load([gpath 'grid_OM4_05_COBALTv3.mat'],'wet',...
    'z_l_units','z_l_long_name','z_l','geolon','geolat')

dz = diff(z_l);
%dz_mat = repmat(dz,1,12);

%% colors
cm10=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0.5 0.5 0.5; ...    %med grey
    0 0 0];...      %black
    
set(groot,'defaultAxesColorOrder',cm10);

%%
[ni,nj]=size(geolon);
geolon = double(geolon);
geolat = double(geolat);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;  

%% 
load([fpath Yr '0101.ocean_cobalt_tracers_month_z.mat'])

%% molN/kg to gC/m3
NtoC= 1035 * (106/16) * 12.01;

%%
% z_l_ts = repmat(z_l,1,length(tmos));
% [z_l2,tts] = meshgrid(tmos,-1*z_l);

%% Time series
tmos = 1:12;

figure(1)
plot(tmos,log10(NtoC*tDe(tmos)+eps),'color',cm10(8,:)); hold on;
plot(tmos,log10(NtoC*tDi(tmos)+eps),'color',cm10(1,:)); hold on;
plot(tmos,log10(NtoC*tSp(tmos)+eps),'color',cm10(3,:)); hold on;
plot(tmos,log10(NtoC*tLp(tmos)+eps),'color',cm10(5,:)); hold on; 
plot(tmos,log10(NtoC*tSz(tmos)+eps),'color',cm10(2,:)); hold on;
legend({'Det','Diaz','SP','LP','SZ'})
legend('location','eastoutside')
title('log_1_0 Integrated Biomass (gC m^-^3)')
xlabel('Months')
stamp(Yr)
print('-dpng',[ppath mod '_ts_mean_ocean_tracers_z.png'])

%% Vert distrib
figure(2)
subplot(1,2,1)
plot(log10(NtoC*vDe(:,1)+eps),-1*z_l,'color',cm10(8,:)); hold on;
plot(log10(NtoC*vDi(:,1)+eps),-1*z_l,'color',cm10(1,:)); hold on;
plot(log10(NtoC*vSp(:,1)+eps),-1*z_l,'color',cm10(3,:)); hold on;
plot(log10(NtoC*vLp(:,1)+eps),-1*z_l,'color',cm10(5,:)); hold on; 
plot(log10(NtoC*vSz(:,1)+eps),-1*z_l,'color',cm10(2,:)); hold on;
legend({'Det','Diaz','SP','LP','SZ'})
legend('location','east')
title('log_1_0 Mean Biomass (gC m^-^3)')
ylabel('Depth (m)')

subplot(1,2,2)
plot((NtoC*vDe(1:8,1)),-1*z_l(1:8),'color',cm10(8,:)); hold on;
plot((NtoC*vDi(1:8,1)),-1*z_l(1:8),'color',cm10(1,:)); hold on;
plot((NtoC*vSp(1:8,1)),-1*z_l(1:8),'color',cm10(3,:)); hold on;
plot((NtoC*vLp(1:8,1)),-1*z_l(1:8),'color',cm10(5,:)); hold on; 
plot((NtoC*vSz(1:8,1)),-1*z_l(1:8),'color',cm10(2,:)); hold on;
% legend({'Det','Diaz','SP','LP','SZ'})
% legend('location','east')
title('Mean Biomass (gC m^-^3)')
ylabel('Depth (m)')
stamp(Yr)
print('-dpng',[ppath mod '_vert_mean_ocean_tracers_z.png'])

%% Maps
% Det
figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NtoC*sDe(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3 0]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean detritus (gC m^-^2)')
stamp(Yr)
print('-dpng',[ppath mod '_global_Det.png'])

%% Diaz
figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NtoC*sDi(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3 0]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean diazotrophs (gC m^-^2)')
stamp(Yr)
print('-dpng',[ppath mod '_global_Diaz.png'])

%% Sp
figure(5)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NtoC*sSp(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 0]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean small phyto (gC m^-^2)')
stamp(Yr)
print('-dpng',[ppath mod '_global_Sp.png'])

%% Lp
figure(6)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NtoC*sLp(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 0.5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean large phyto (gC m^-^2)')
stamp(Yr)
print('-dpng',[ppath mod '_global_Lp.png'])

%% Sz
figure(7)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NtoC*sSz(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-0.5 0.5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean small zoo (gC m^-^2)')
stamp(Yr)
print('-dpng',[ppath mod '_global_Sz.png'])

%% Vert distrib over time
% figure(1)
% subplot(3,3,1)
% pcolor(z_l2(1:10,:),tts(1:10,:),log10(mNO2(1:10,:)+eps));
% shading flat;
% colorbar
% %clim([-10 -6])
% title('NO3')
% 
% subplot(3,3,2)
% pcolor(z_l2(1:10,:),tts(1:10,:),log10(mNH2(1:10,:)+eps));
% shading flat;
% colorbar
% %clim([-12 -9])
% title('NH4')
% 
% subplot(3,3,3)
% pcolor(z_l2(1:10,:),tts(1:10,:),log10(mO2(1:10,:)+eps));
% shading flat;
% colorbar
% %clim([-14 -9])
% title('O2')
% 
% subplot(3,3,4)
% pcolor(z_l2(1:10,:),tts(1:10,:),log10(mB2(1:10,:)+eps));
% shading flat;
% colorbar
% clim([-10 -6])
% title('Bact')
% ylabel('Depth (m)')
% 
% subplot(3,3,5)
% pcolor(z_l2(1:10,:),tts(1:10,:),log10(mDP2(1:10,:)+eps));
% shading flat;
% colorbar
% clim([-12 -9])
% title('Diaz')
% 
% subplot(3,3,6)
% pcolor(z_l2(1:10,:),tts(1:10,:),log10(mSp2(1:10,:)+eps));
% shading flat;
% colorbar
% clim([-10 -6])
% title('Sp')
% 
% subplot(3,3,7)
% pcolor(z_l2(1:10,:),tts(1:10,:),log10(mMP2(1:10,:)+eps));
% shading flat;
% colorbar
% clim([-10 -6])
% title('MP')
% 
% subplot(3,3,8)
% pcolor(z_l2(1:10,:),tts(1:10,:),log10(mLp2(1:10,:)+eps));
% shading flat;
% colorbar
% clim([-10 -6])
% title('Lp')
% xlabel('Time (mo)')
% 
% subplot(3,3,9)
% pcolor(z_l2(1:10,:),tts(1:10,:),log10(mSz2(1:10,:)+eps));
% shading flat;
% colorbar
% clim([-10 -6])
% title('Sz')
% xlabel('Time (mo)')
% 
%print('-dpng',[ppath mod '_depth_ts_phyto_sz.png'])



