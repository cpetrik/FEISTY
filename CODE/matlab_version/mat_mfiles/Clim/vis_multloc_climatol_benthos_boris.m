% Visualize output of FEISTY
% ESM2.6 Climatology of 5 yrs
% 150 years
% Last year saved as mat file

clear all
close all

Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';

cpath = ['/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%Orig: cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
cfile='NoDc_enc70-b200_m4-b175-k086_c20-b250_nmort1_2B_BE75_noCC_RE00100';
%harv = 'All_fish03';
harv = 'pristine_vGEB';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%load([fpath 'Means_Climatol_' harv '_' cfile '.mat']);
load([fpath 'Climatol_' harv '.mat']);

%% Means
%Time
sd_tmean=mean(Spinup_Sml_d.bio,1);
md_tmean=mean(Spinup_Med_d.bio,1);
ld_tmean=mean(Spinup_Lrg_d.bio,1);
sb_tmean=mean(Spinup_SBent.bio,1);
mb_tmean=mean(Spinup_MBent.bio,1);

% Last year means
[id,nt] = size(Spinup_SBent.bio);
time=1:nt;
lyr=time((end-12+1):end);

sd_mean=mean(Spinup_Sml_d.bio(:,lyr),2);
md_mean=mean(Spinup_Med_d.bio(:,lyr),2);
ld_mean=mean(Spinup_Lrg_d.bio(:,lyr),2);
sb_mean=mean(Spinup_SBent.bio(:,lyr),2);
mb_mean=mean(Spinup_MBent.bio(:,lyr),2);

sd_prod=mean(Spinup_Sml_d.prod(:,lyr),2);
md_prod=mean(Spinup_Med_d.prod(:,lyr),2);
ld_prod=mean(Spinup_Lrg_d.prod(:,lyr),2);

save([fpath 'Means_Climatol_' harv '.mat'],'time','lyr',...
    'sd_tmean','md_tmean','ld_tmean','sb_tmean','mb_tmean',...
    'sd_mean','md_mean','ld_mean','sb_mean','mb_mean',...
    'sd_prod','md_prod','ld_prod');

%% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

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
   
set(groot,'defaultAxesColorOrder',cm9);

%% Plots in time
y = time;
nt = length(time);

% All size classes of all
figure(1)
plot(y,log10(sd_tmean),'Linewidth',1); hold on;
plot(y,log10(md_tmean),'Linewidth',1); hold on;
plot(y,log10(ld_tmean),'Linewidth',1); hold on;
plot(y,log10(sb_tmean),'Linewidth',1); hold on;
plot(y,log10(mb_tmean),'Linewidth',1); hold on;
legend('SD','MD','LD','SB','MB')
legend('location','eastoutside')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
title(['Climatol ' harv])
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_all_sizes.png'])

%%
figure(2)
D = sd_tmean+md_tmean+ld_tmean;
B = sb_tmean+mb_tmean;

plot(y,log10(D),'color',[0 0.7 0.2],'Linewidth',2); hold on;
plot(y,log10(B),'k','Linewidth',2); hold on;
legend('D','B')
legend('location','eastoutside')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
title(['Climatol ' harv])
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_all_types.png'])

%% Plots in space
Zsd=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zsb=NaN*ones(ni,nj);
Zmb=NaN*ones(ni,nj);

Zsd(ID)=sd_mean;
Zmd(ID)=md_mean;
Zld(ID)=ld_mean;
Zsb(ID)=sb_mean;
Zmb(ID)=mb_mean;

AllD = Zsd+Zmd+Zld;
AllB = Zsb+Zmb;

%% bent
figure(3)
subplot('position',[0.05 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zsb+eps))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 2]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Climatology log10 mean Small benthos (g m^-^2)')

%
subplot('position',[0.05 0.0 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zmb+eps))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 2]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Climatology log10 mean Medium benthos (g m^-^2)')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_global_BENT.png'])


%% sd
figure(4)
subplot('position',[0.0 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zsd))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 -5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('SD')

% md
subplot('position',[0.5 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zmd))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 -5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('MD')

% ld
subplot('position',[0.0 0.0 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zld))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 -5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('LD')

%all
subplot('position',[0.5 0.0 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 -5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('All D')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_global_AllD.png'])

