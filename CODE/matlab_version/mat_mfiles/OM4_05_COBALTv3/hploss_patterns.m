% Look at COBALT-FEISTY forcing from online sim
% See if relationships between HPloss and Z biomass and temp

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

%% temp-dep fn

tTdep = exp(0.63*(tTp-10));
vTdep = exp(0.63*(vTp-10));
sTdep = exp(0.63*(sTp-10));
mTdep = exp(0.63*(mTP-10));


%% Start simple - time
figure(1)
subplot(1,2,1)
scatter(tMz,tMhp,'k','filled')
ylabel('MZhploss')
xlabel('MZbio')

subplot(1,2,2)
scatter(tTp,tMhp,'k','filled')
ylabel('MZhploss')
xlabel('Tpel')


figure(2)
subplot(1,2,1)
scatter(tLz,tLhp,'k','filled')
ylabel('LZhploss')
xlabel('LZbio')

subplot(1,2,2)
scatter(tTp,tLhp,'k','filled')
ylabel('LZhploss')
xlabel('Tpel')

figure(3)
subplot(1,2,1)
scatter(tMz,tMhp./tTdep,'k','filled')
ylabel('MZhploss/Tdep')
xlabel('MZbio')

subplot(1,2,2)
scatter(tLz,tLhp./tTdep,'k','filled')
ylabel('LZhploss/Tdep')
xlabel('LZbio')

%% Vertical

%% Space
figure(7)
subplot(1,2,1)
scatter(sMz,sMhp,'k')
ylabel('MZhploss')
xlabel('MZbio')

subplot(1,2,2)
scatter(sTp,sMhp,'k')
ylabel('MZhploss')
xlabel('Tpel')


figure(8)
subplot(1,2,1)
scatter(sLz,sLhp,'k')
ylabel('LZhploss')
xlabel('LZbio')

subplot(1,2,2)
scatter(sTp,sLhp,'k')
ylabel('LZhploss')
xlabel('Tpel')

figure(9)
subplot(1,2,1)
scatter(sMz,sMhp./sTdep,'k')
ylabel('MZhploss/Tdep')
xlabel('MZbio')

subplot(1,2,2)
scatter(sLz,sLhp./sTdep,'k')
ylabel('LZhploss/Tdep')
xlabel('LZbio')

%% Space and time