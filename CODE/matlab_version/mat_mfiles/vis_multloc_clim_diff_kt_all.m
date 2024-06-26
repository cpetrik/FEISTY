% Visualize output of POEM
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/Bio_rates/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

close all

% plot info
[ni,nj]=size(lon);
geolat_t=lat;
geolon_t=lon;
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));


% colors
load('MyColormaps.mat')
cm9=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...  %b
    0 0 0];...      %black
    
cm21=[1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    0 1 1;...     %c
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.75 0.75 0.75;... %lt grey
    0.5 0.5 0.5;...    %med grey
    49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    1 1 0;...      %yellow
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255]; %peach

set(groot,'defaultAxesColorOrder',cm9);

%% Plots in space

cfileA = 'Dc_enc70-b200_cm20_m-b175-k06_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
fpathA=['/Volumes/GFDL/CSV/Matlab_new_size/' cfileA '/'];
load([fpathA 'Clim_means_All_fish03.mat'],'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Afish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Aall=NaN*ones(ni,nj);
Aall(ID)=Afish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

cfileB = 'Dc_enc70-b200_cm20_m-b175-k08_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
fpathB=['/Volumes/GFDL/CSV/Matlab_new_size/' cfileB '/'];
load([fpathB 'Clim_means_All_fish03.mat'],'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Bfish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Ball=NaN*ones(ni,nj);
Ball(ID)=Bfish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

cfileC = 'Dc_enc70-b200_cm20_m-b175-k10_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
fpathC=['/Volumes/GFDL/CSV/Matlab_new_size/' cfileC '/'];
load([fpathC 'Clim_means_All_fish03.mat'],'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Cfish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Call=NaN*ones(ni,nj);
Call(ID)=Cfish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

cfileD = 'Dc_enc70-b200_cm20_m-b175-k12_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
fpathD=['/Volumes/GFDL/CSV/Matlab_new_size/' cfileD '/'];
load([fpathD 'Clim_means_All_fish03.mat'],'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Dfish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Dall=NaN*ones(ni,nj);
Dall(ID)=Dfish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

%% Loop over kts
% POEM file info
harv = 'All_fish03';
kays = 0.0405:0.01:0.1205;
skays={'.04','.05','.06','.07','.08','.09','.10','.11','.12'};

fish = NaN*ones(length(ID),length(kays));

for i=1:length(kays)
    kt=kays(i);
    tkfn = num2str(100+int64(100*kt));
    cfile = ['Dc_enc70-b200_cm20_m-b175-k',tkfn(2:end),...
        '_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100'];
    fpath=['/Volumes/GFDL/CSV/Matlab_new_size/' cfile '/'];
    load([fpath 'Clim_means_All_fish03.mat'],...
    'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');

    fish(:,i) = sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+...
        lp_mean+ld_mean;

    clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean
end
%%
A1=NaN*ones(ni,nj);
A2=NaN*ones(ni,nj);
A3=NaN*ones(ni,nj);
A4=NaN*ones(ni,nj);
A5=NaN*ones(ni,nj);
A6=NaN*ones(ni,nj);
A7=NaN*ones(ni,nj);
A8=NaN*ones(ni,nj);
A9=NaN*ones(ni,nj);
A1(ID)=fish(:,1);
A2(ID)=fish(:,2);
A3(ID)=fish(:,3);
A4(ID)=fish(:,4);
A5(ID)=fish(:,5);
A6(ID)=fish(:,6);
A7(ID)=fish(:,7);
A8(ID)=fish(:,8);
A9(ID)=fish(:,9);

%% 
figure(1)
%A
subplot('Position',[0 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(A3))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'A')

%B
subplot('Position',[0.5 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(A5))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'B')

%C
subplot('Position',[0 0.1 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(A7))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'C')

%D
subplot('Position',[0.5 0.1 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(A9))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.25 0.075 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
text(-2.75,1.75,'D')
%stamp([harv '_' cfile])
print('-dpng',[pp 'Climatol_' harv '_All_comp_kt_6-12.png'])

%%
figure(3)
%A
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(A1))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'A')
text(-0.75,1.75,'kt=0.04')
%B
subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(A3))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'B')
text(-0.75,1.75,'kt=0.06')
%C
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(A5))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'C')
text(-0.75,1.75,'kt=0.08')
%D
subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(A7))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'D')
text(-0.75,1.75,'kt=0.10')
%E
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(A9))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'E')
text(-0.75,1.75,'kt=0.12')
colorbar('Position',[0.25 0.04 0.34 0.025],'orientation','horizontal')
stamp('')
%F
% subplot('Position',[0.41 0.06 0.4 0.3])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(A6))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% set(gcf,'renderer','painters')
print('-dpng',[pp 'Climatol_' harv '_All_comp_kt_4-12.png'])

