% Look at COBALT forcing from offline sim

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
load([fpath 'ocean_cobalt_feisty_forcing_z.199001-199912_means.mat'])

%% molN/kg to gC/m3
S2D= 60*60*24;
%MZ and LZ molN/kg
N2Ckg= 1035 * (106/16) * 12.01;
%HPloss 'mol N kg-1 s-1'
N2CkgD = N2Ckg * S2D;
%fn_tot_btm = 'mol m-2 s-1'
N2CmD = (106/16) * 12.01 * S2D;

%%
% z_l_ts = repmat(z_l,1,length(tmos));
% [z_l2,tts] = meshgrid(tmos,-1*z_l);

%% Time series
tmos = 1:120;

figure(1)
plot(tmos,log10(N2CmD*tDE(tmos)+eps),'color',cm10(1,:)); hold on;
plot(tmos,log10(N2Ckg*tMZ(tmos)+eps),'color',cm10(4,:)); hold on;
plot(tmos,log10(N2Ckg*tLZ(tmos)+eps),'color',cm10(5,:)); hold on;
plot(tmos,log10(N2CkgD*tMH(tmos)+eps),'color',cm10(6,:)); hold on; 
plot(tmos,log10(N2CkgD*tLH(tmos)+eps),'color',cm10(7,:)); hold on;
legend({'DetBtm','MZbio','LZbio','MZhploss','LZhploss'})
legend('location','eastoutside')
title('log_1_0 Integrated Biomass/Flux (gC m^-^2 or gC m^-^2 d^-^1)')
xlabel('Months')
stamp('')
print('-dpng',[ppath exper '_ts_mean_feisty_forcing_z.png'])

%% Vert distrib
figure(2)
subplot(1,2,1)
plot(log10(N2Ckg*vMZ(:,1)+eps),-1*z_l,'color',cm10(4,:)); hold on;
plot(log10(N2Ckg*vLZ(:,1)+eps),-1*z_l,'color',cm10(5,:)); hold on;
plot(log10(N2CkgD*vMH(:,1)+eps),-1*z_l,'color',cm10(6,:)); hold on; 
plot(log10(N2CkgD*vLH(:,1)+eps),-1*z_l,'color',cm10(7,:)); hold on;
legend({'MZbio','LZbio','MZhploss','LZhploss'})
legend('location','east')
title('log_1_0 Mean Biomass/Flux (gC m^-^3 or gC m^-^3 d^-^1)')
ylabel('Depth (m)')

subplot(1,2,2)
plot((N2Ckg*vMZ(1:8,1)),-1*z_l(1:8),'color',cm10(4,:)); hold on;
plot((N2Ckg*vLZ(1:8,1)),-1*z_l(1:8),'color',cm10(5,:)); hold on;
plot((N2CkgD*vMH(1:8,1)),-1*z_l(1:8),'color',cm10(6,:)); hold on; 
plot((N2CkgD*vLH(1:8,1)),-1*z_l(1:8),'color',cm10(7,:)); hold on;
% legend({'Det','Diaz','SP','LP','SZ'})
% legend('location','east')
title('Mean Biomass/Flux (gC m^-^3 or gC m^-^3 d^-^1)')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath exper '_vert_mean_feisty_forcing_z.png'])

%% Maps
% Det
figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(N2CmD*sDE(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3 -1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean btm detritus flux (gC m^-^2 d^-^1)')
stamp('')
print('-dpng',[ppath exper '_global_DetBtm.png'])

%% MZ
figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(N2Ckg*sMZ(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 0.5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean medium zoo (gC m^-^2)')
stamp('')
print('-dpng',[ppath exper '_global_MZ.png'])

%% LZ
figure(5)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(N2Ckg*sLZ(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean large zoo (gC m^-^2)')
stamp('')
print('-dpng',[ppath exper '_global_LZ.png'])

%% MZ HPloss
figure(6)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(N2CkgD*sMH(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-4 -2]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean MZ HPloss (gC m^-^2 d^-^1)')
stamp('')
print('-dpng',[ppath exper '_global_MZhploss.png'])

%% LZ HPloss
figure(7)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(N2CkgD*sLH(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-4 -1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean LZ HPloss (gC m^-^2 d^-^1)')
stamp('')
print('-dpng',[ppath exper '_global_LZhploss.png'])

%% Vert distrib over time
figure(1)
subplot(3,3,1)
pcolor(z_l2(1:10,:),tts(1:10,:),log10(mNO2(1:10,:)+eps));
shading flat;
colorbar
%clim([-10 -6])
title('NO3')

subplot(3,3,2)
pcolor(z_l2(1:10,:),tts(1:10,:),log10(mNH2(1:10,:)+eps));
shading flat;
colorbar
%clim([-12 -9])
title('NH4')

subplot(3,3,3)
pcolor(z_l2(1:10,:),tts(1:10,:),log10(mO2(1:10,:)+eps));
shading flat;
colorbar
%clim([-14 -9])
title('O2')

subplot(3,3,4)
pcolor(z_l2(1:10,:),tts(1:10,:),log10(mB2(1:10,:)+eps));
shading flat;
colorbar
clim([-10 -6])
title('Bact')
ylabel('Depth (m)')

subplot(3,3,5)
pcolor(z_l2(1:10,:),tts(1:10,:),log10(mDP2(1:10,:)+eps));
shading flat;
colorbar
clim([-12 -9])
title('Diaz')

subplot(3,3,6)
pcolor(z_l2(1:10,:),tts(1:10,:),log10(mSP2(1:10,:)+eps));
shading flat;
colorbar
clim([-10 -6])
title('SP')

subplot(3,3,7)
pcolor(z_l2(1:10,:),tts(1:10,:),log10(mMP2(1:10,:)+eps));
shading flat;
colorbar
clim([-10 -6])
title('MP')

subplot(3,3,8)
pcolor(z_l2(1:10,:),tts(1:10,:),log10(mLP2(1:10,:)+eps));
shading flat;
colorbar
clim([-10 -6])
title('LP')
xlabel('Time (mo)')

subplot(3,3,9)
pcolor(z_l2(1:10,:),tts(1:10,:),log10(mSZ2(1:10,:)+eps));
shading flat;
colorbar
clim([-10 -6])
title('SZ')
xlabel('Time (mo)')


%print('-dpng',[ppath exper '_OSP_depth_ts_phyto_nuts.png'])



