% Look at COBALT-FEISTY tracers from online sim

clear
close all

%%
fpath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/';

gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

cfile ='NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/',...
    cfile,'/Cobalt_Feisty/'];

%%
% load([gpath 'grid_OM4_05_COBALTv3.mat'],'wet','z_l_units','z_l_long_name','z_l')
% 
% dz = diff(z_i);
% dz_mat = repmat(dz,1,12);

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
load([fpath 'ocean_cobalt_tracers_month_z.199001-199412_means.nc'])

%%
z_l_ts = repmat(z_l,1,length(tmos));
[z_l2,tts] = meshgrid(tmos,-1*z_l);

%% Plots
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


print('-dpng',[ppath exper '_OSP_depth_ts_phyto_nuts.png'])

%%
figure(3)
subplot(3,1,1)
plot(tmos,log10(tB2+eps),'color',cm10(8,:)); hold on;
plot(tmos,log10(tDP2+eps),'color',cm10(1,:)); hold on;
plot(tmos,log10(tSP2+eps),'color',cm10(3,:)); hold on;
plot(tmos,log10(tMP2+eps),'color',cm10(4,:)); hold on;
plot(tmos,log10(tLP2+eps),'color',cm10(5,:)); hold on; 
plot(tmos,log10(tSZ2+eps),'color',cm10(2,:)); hold on;
legend({'B','Di','SP','MP','LP','SZ'})
legend('location','eastoutside')
title('log_1_0 Integrated Biomass')

subplot(3,1,2)
plot(tmos,log10(tNO2+eps),'color',cm10(9,:)); hold on;
plot(tmos,log10(tNH2+eps),'color',cm10(10,:)); hold on;
legend({'NO3','NH4'})
legend('location','eastoutside')
xlabel('Time (d)')
title('log_1_0 N weighted mean')
ylim([-8 -4.5])

subplot(3,1,3)
plot(tmos,tO2,'color',cm10(2,:)); hold on;
xlabel('Time (d)')
title('O2 weighted mean')
print('-dpng',[ppath exper '_OSP_ts_mean_phyto_nuts_subplot.png'])

