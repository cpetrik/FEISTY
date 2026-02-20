% Look at COBALT only tracers from online sim

clear
close all

%%
fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/',...
    'FEISTY/CODE/Figs/GCM_figs/OM4_05_COBALTv3/'];

exper = 'OM4_05_COBALTv3_FEISTYoff';

%%
load([fpath 'grid_OM4_05_COBALTv3.mat'],'wet',...
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
load([fpath 'ocean_cobalt_ocean_tracers_z.199001-199412_means.mat'])

%% molN/kg to gC/m3
NtoC= 1035 * (106/16) * 12.01;

%%
% z_l_ts = repmat(z_l,1,length(tmos));
% [z_l2,tts] = meshgrid(tmos,-1*z_l);

%% Time series
tmos = 1:60;

figure(1)
% plot(tmos,log10(NtoC*tDE(tmos)+eps),'color',cm10(8,:)); hold on;
% plot(tmos,log10(NtoC*tDI(tmos)+eps),'color',cm10(1,:)); hold on;
plot(tmos,log10(NtoC*tSP(tmos)+eps),'color',cm10(3,:)); hold on;
plot(tmos,log10(NtoC*tLP(tmos)+eps),'color',cm10(5,:)); hold on; 
plot(tmos,log10(NtoC*tSZ(tmos)+eps),'color',cm10(2,:)); hold on;
%legend({'Det','Diaz','SP','LP','SZ'})
legend({'SP','LP','SZ'})
legend('location','eastoutside')
title('log_1_0 Integrated Biomass (gC m^-^3)')
xlabel('Months')
stamp('')
print('-dpng',[ppath exper '_ts_mean_ocean_tracers_z.png'])

%% Vert distrib
figure(2)
subplot(1,2,1)
% plot(log10(NtoC*vDE(:,1)+eps),-1*z_l,'color',cm10(8,:)); hold on;
% plot(log10(NtoC*vDI(:,1)+eps),-1*z_l,'color',cm10(1,:)); hold on;
plot(log10(NtoC*vSP(:,1)+eps),-1*z_l,'color',cm10(3,:)); hold on;
plot(log10(NtoC*vLP(:,1)+eps),-1*z_l,'color',cm10(5,:)); hold on; 
plot(log10(NtoC*vSZ(:,1)+eps),-1*z_l,'color',cm10(2,:)); hold on;
legend({'SP','LP','SZ'})
%legend({'Det','Diaz','SP','LP','SZ'})
legend('location','east')
title('log_1_0 Mean Biomass (gC m^-^3)')
ylabel('Depth (m)')

subplot(1,2,2)
% plot((NtoC*vDE(1:8,1)),-1*z_l(1:8),'color',cm10(8,:)); hold on;
% plot((NtoC*vDI(1:8,1)),-1*z_l(1:8),'color',cm10(1,:)); hold on;
plot((NtoC*vSP(1:8,1)),-1*z_l(1:8),'color',cm10(3,:)); hold on;
plot((NtoC*vLP(1:8,1)),-1*z_l(1:8),'color',cm10(5,:)); hold on; 
plot((NtoC*vSZ(1:8,1)),-1*z_l(1:8),'color',cm10(2,:)); hold on;
% legend({'Det','Diaz','SP','LP','SZ'})
% legend('location','east')
title('Mean Biomass (gC m^-^3)')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath exper '_vert_mean_ocean_tracers_z.png'])

%% Maps
% Det
% figure(3)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat,geolon,log10(squeeze(NtoC*sDE(:,:,1))))
% cmocean('matter')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-3 0]);
% hcb = colorbar('h');
% set(gcf,'renderer','painters')
% title('log10 mean detritus (gC m^-^2)')
% stamp(cfile)
% print('-dpng',[ppath exper '_global_Det.png'])

%% Diaz
% figure(4)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat,geolon,log10(squeeze(NtoC*sDI(:,:,1))))
% cmocean('matter')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-3 0]);
% hcb = colorbar('h');
% set(gcf,'renderer','painters')
% title('log10 mean diazotrophs (gC m^-^2)')
% stamp(cfile)
% print('-dpng',[ppath exper '_global_Diaz.png'])

%% SP
figure(5)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NtoC*sSP(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 0]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean small phyto (gC m^-^2)')
stamp(cfile)
print('-dpng',[ppath exper '_global_SP.png'])

%% LP
figure(6)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NtoC*sLP(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 0.5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean large phyto (gC m^-^2)')
stamp(cfile)
print('-dpng',[ppath exper '_global_LP.png'])

%% SZ
figure(7)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NtoC*sSZ(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-0.5 0.5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean small zoo (gC m^-^2)')
stamp(cfile)
print('-dpng',[ppath exper '_global_SZ.png'])

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
% pcolor(z_l2(1:10,:),tts(1:10,:),log10(mSP2(1:10,:)+eps));
% shading flat;
% colorbar
% clim([-10 -6])
% title('SP')
% 
% subplot(3,3,7)
% pcolor(z_l2(1:10,:),tts(1:10,:),log10(mMP2(1:10,:)+eps));
% shading flat;
% colorbar
% clim([-10 -6])
% title('MP')
% 
% subplot(3,3,8)
% pcolor(z_l2(1:10,:),tts(1:10,:),log10(mLP2(1:10,:)+eps));
% shading flat;
% colorbar
% clim([-10 -6])
% title('LP')
% xlabel('Time (mo)')
% 
% subplot(3,3,9)
% pcolor(z_l2(1:10,:),tts(1:10,:),log10(mSZ2(1:10,:)+eps));
% shading flat;
% colorbar
% clim([-10 -6])
% title('SZ')
% xlabel('Time (mo)')
% 
%print('-dpng',[ppath exper '_OSP_depth_ts_phyto_nuts.png'])



