% Look at COBALT-FEISTY forcing from online sim

clear
close all

%%
fpath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/cobalt_feisty/';

gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

cfile ='NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/',...
    cfile,'/Cobalt_Feisty/'];

exper = 'OM4_05_COBALTv3_FEISTYon_021326';

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
load([fpath 'ocean_cobalt_feisty_forcing_1990_means.mat'])

%% molN/kg to gC/m3
S2D= 60*60*24;
%MZ and LZ molN/kg
N2Ckg= 1035 * (106/16) * 12.01;
%HPloss 'mol N kg-1 s-1'
N2CkgD = N2Ckg * S2D;
%fn_tot_btm = 'mol m-2 s-1'
N2CmD = (106/16) * 12.01 * S2D;

%%
tmos = 1:12;
z_l_ts = repmat(z_l,1,length(tmos));
[z_l2,tts] = meshgrid(tmos,-1*z_l);

%% Time series

figure(1)
plot(tmos,log10(N2CmD*tDet(tmos)+eps),'color',cm10(1,:)); hold on;
plot(tmos,log10(N2Ckg*tMz(tmos)+eps),'color',cm10(4,:)); hold on;
plot(tmos,log10(N2Ckg*tLz(tmos)+eps),'color',cm10(5,:)); hold on;
plot(tmos,log10(N2CkgD*tMhp(tmos)+eps),'color',cm10(6,:)); hold on; 
plot(tmos,log10(N2CkgD*tLhp(tmos)+eps),'color',cm10(7,:)); hold on;
legend({'DetBtm','MZbio','LZbio','MZhploss','LZhploss'})
legend('location','eastoutside')
title('log_1_0 Integrated Biomass/Flux (gC m^-^2 or gC m^-^2 d^-^1)')
xlabel('Months')
stamp('')
print('-dpng',[ppath exper '_ts_mean_feisty_forcing_z.png'])

%% Vert distrib
figure(2)
subplot(1,2,1)
plot(log10(N2Ckg*vMz1(:,1)+eps),-1*z_l,'color',cm10(4,:)); hold on;
plot(log10(N2Ckg*vLz1(:,1)+eps),-1*z_l,'color',cm10(5,:)); hold on;
plot(log10(N2CkgD*vMhp1(:,1)+eps),-1*z_l,'color',cm10(6,:)); hold on; 
plot(log10(N2CkgD*vLhp1(:,1)+eps),-1*z_l,'color',cm10(7,:)); hold on;
legend({'MZbio','LZbio','MZhploss','LZhploss'})
legend('location','east')
title('log_1_0 Mean Biomass/Flux (gC m^-^3 or gC m^-^3 d^-^1)')
ylabel('Depth (m)')

subplot(1,2,2)
plot((N2Ckg*vMz1(1:8,1)),-1*z_l(1:8),'color',cm10(4,:)); hold on;
plot((N2Ckg*vLz1(1:8,1)),-1*z_l(1:8),'color',cm10(5,:)); hold on;
plot((N2CkgD*vMhp1(1:8,1)),-1*z_l(1:8),'color',cm10(6,:)); hold on; 
plot((N2CkgD*vLhp1(1:8,1)),-1*z_l(1:8),'color',cm10(7,:)); hold on;
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
surfm(geolat,geolon,log10(squeeze(N2CmD*sDet)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3 -1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean btm detritus flux (gC m^-^2 d^-^1)')
stamp(cfile)
print('-dpng',[ppath exper '_global_DetBtm.png'])

%% MZ
figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(N2Ckg*sMz)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 0.5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean medium zoo (gC m^-^2)')
stamp(cfile)
print('-dpng',[ppath exper '_global_MZ.png'])

%% LZ
figure(5)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(N2Ckg*sLz)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean large zoo (gC m^-^2)')
stamp(cfile)
print('-dpng',[ppath exper '_global_LZ.png'])

%% MZ HPloss
figure(6)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(N2CkgD*sMhp)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-4 -2]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean MZ HPloss (gC m^-^2 d^-^1)')
stamp(cfile)
print('-dpng',[ppath exper '_global_MZhploss.png'])

%% LZ HPloss
figure(7)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(N2CkgD*sLhp)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-4 -1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean LZ HPloss (gC m^-^2 d^-^1)')
stamp(cfile)
print('-dpng',[ppath exper '_global_LZhploss.png'])

%% Vert distrib over time
figure(8)
subplot(2,2,1)
pcolor(z_l2(1:10,:),tts(1:10,:),log10(N2Ckg*vMz2(1:10,:)+eps));
shading flat;
colorbar
clim([-10 -6])
title('MZ bio (gC m^-^3)')

subplot(2,2,2)
pcolor(z_l2(1:10,:),tts(1:10,:),log10(N2Ckg*vLz2(1:10,:)+eps));
shading flat;
colorbar
clim([-10 -6])
title('LZ bio')

subplot(2,2,3)
pcolor(z_l2(1:10,:),tts(1:10,:),log10(N2CkgD*vMhp2(1:10,:)+eps));
shading flat;
colorbar
clim([-10 -6])
title('MZ hploss (gC m^-^3 d^-^1)')

subplot(2,2,4)
pcolor(z_l2(1:10,:),tts(1:10,:),log10(N2CkgD*vLhp2(1:10,:)+eps));
shading flat;
colorbar
clim([-10 -6])
title('LZ hploss')
ylabel('Depth (m)')
xlabel('Time (mo)')

%print('-dpng',[ppath exper '_OSP_depth_ts_phyto_nuts.png'])



