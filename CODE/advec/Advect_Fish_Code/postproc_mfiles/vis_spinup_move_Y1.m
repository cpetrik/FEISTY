% Visualize output of Spinup Y1
% Saved as mat files

clear 
close all

%%
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/CORE/'];
ppath = [pp cfile '/CORE/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%%
mod = 'Spinup1988_move_prey_v13_zerovel_All_fish03_Y1';
load([fpath mod '.mat']);
%load([fpath 'Means_' exper cfile '.mat']);

save([fpath mod '.mat'],'S_Bent_bio','S_Sml_f','S_Sml_p','S_Sml_d','S_Med_f',...
    'S_Med_p','S_Med_d','S_Lrg_p','S_Lrg_d','GRD1','GRD2','exper');

%%
load('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');
load('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

[ni,nj]=size(geolon_t);
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; 
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

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

%% Take means
time = 1:365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];


%Time
sp_tmean=mean(S_Sml_p,1);
sf_tmean=mean(S_Sml_f,1);
sd_tmean=mean(S_Sml_d,1);
mp_tmean=mean(S_Med_p,1);
mf_tmean=mean(S_Med_f,1);
md_tmean=mean(S_Med_d,1);
lp_tmean=mean(S_Lrg_p,1);
ld_tmean=mean(S_Lrg_d,1);
b_tmean=mean(S_Bent_bio,1);

% sp_tmean=mean(S_Sml_p,1,'omitnan');
% sf_tmean=mean(S_Sml_f,1,'omitnan');
% sd_tmean=mean(S_Sml_d,1,'omitnan');
% mp_tmean=mean(S_Med_p,1,'omitnan');
% mf_tmean=mean(S_Med_f,1,'omitnan');
% md_tmean=mean(S_Med_d,1,'omitnan');
% lp_tmean=mean(S_Lrg_p,1,'omitnan');
% ld_tmean=mean(S_Lrg_d,1,'omitnan');
% b_tmean=mean(S_Bent_bio,1,'omitnan');


% Space
%start 
sp_mean1=mean(S_Sml_p(:,1),2,'omitnan');
sf_mean1=mean(S_Sml_f(:,1),2,'omitnan');
sd_mean1=mean(S_Sml_d(:,1),2,'omitnan');
mp_mean1=mean(S_Med_p(:,1),2,'omitnan');
mf_mean1=mean(S_Med_f(:,1),2,'omitnan');
md_mean1=mean(S_Med_d(:,1),2,'omitnan');
lp_mean1=mean(S_Lrg_p(:,1),2,'omitnan');
ld_mean1=mean(S_Lrg_d(:,1),2,'omitnan');
b_mean1=mean(S_Bent_bio(:,1),2,'omitnan');

%middle 
% sp_mean2=mean(S_Sml_p(:,10),2);
% sf_mean2=mean(S_Sml_f(:,10),2);
% sd_mean2=mean(S_Sml_d(:,10),2);
% mp_mean2=mean(S_Med_p(:,10),2);
% mf_mean2=mean(S_Med_f(:,10),2);
% md_mean2=mean(S_Med_d(:,10),2);
% lp_mean2=mean(S_Lrg_p(:,10),2);
% ld_mean2=mean(S_Lrg_d(:,10),2);
% b_mean2=mean(S_Bent_bio(:,10),2);

sp_mean2=mean(S_Sml_p(:,183),2);
sf_mean2=mean(S_Sml_f(:,183),2);
sd_mean2=mean(S_Sml_d(:,183),2);
mp_mean2=mean(S_Med_p(:,183),2);
mf_mean2=mean(S_Med_f(:,183),2);
md_mean2=mean(S_Med_d(:,183),2);
lp_mean2=mean(S_Lrg_p(:,183),2);
ld_mean2=mean(S_Lrg_d(:,183),2);
b_mean2=mean(S_Bent_bio(:,183),2);

%end
sp_mean3=mean(S_Sml_p(:,365),2);
sf_mean3=mean(S_Sml_f(:,365),2);
sd_mean3=mean(S_Sml_d(:,365),2);
mp_mean3=mean(S_Med_p(:,365),2);
mf_mean3=mean(S_Med_f(:,365),2);
md_mean3=mean(S_Med_d(:,365),2);
lp_mean3=mean(S_Lrg_p(:,365),2);
ld_mean3=mean(S_Lrg_d(:,365),2);
b_mean3=mean(S_Bent_bio(:,365),2);

%whole yr
sp_mean=mean(S_Sml_p,2);
sf_mean=mean(S_Sml_f,2);
sd_mean=mean(S_Sml_d,2);
mp_mean=mean(S_Med_p,2);
mf_mean=mean(S_Med_f,2);
md_mean=mean(S_Med_d,2);
lp_mean=mean(S_Lrg_p,2);
ld_mean=mean(S_Lrg_d,2);
b_mean=mean(S_Bent_bio,2);

%% Plots in time
y = time;
F = sf_tmean+mf_tmean;
P = sp_tmean+mp_tmean+lp_tmean;
D = sd_tmean+md_tmean+ld_tmean;
B = b_tmean;

% All size classes of all
figure(1)
plot(y,log10(b_tmean),'Linewidth',1); hold on;
plot(y,log10(sf_tmean),'Linewidth',1); hold on;
plot(y,log10(mf_tmean),'Linewidth',1); hold on;
plot(y,log10(sp_tmean),'Linewidth',1); hold on;
plot(y,log10(mp_tmean),'Linewidth',1); hold on;
plot(y,log10(lp_tmean),'Linewidth',1); hold on;
plot(y,log10(sd_tmean),'Linewidth',1); hold on;
plot(y,log10(md_tmean),'Linewidth',1); hold on;
plot(y,log10(ld_tmean),'Linewidth',1); hold on;
legend('B','SF','MF','SP','MP','LP','SD','MD','LD')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-5 2])
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
title('Spinup')
stamp(exper)
print('-dpng',[ppath mod '_all_sizes.png'])

figure(2)
plot(y,log10(B),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,log10(F),'r','Linewidth',2); hold on;
plot(y,log10(P),'b','Linewidth',2); hold on;
plot(y,log10(D),'k','Linewidth',2); hold on;
legend('B','F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-5 2])
xlabel('Time (y)')
ylabel('log10 Biomass (g m^-^2)')
title(['Spinup'])
print('-dpng',[ppath mod '_all_types.png'])

%% Plots in space
Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zb=NaN*ones(ni,nj);

Zsf(GRD.ID)=sf_mean1;
Zsp(GRD.ID)=sp_mean1;
Zsd(GRD.ID)=sd_mean1;
Zmf(GRD.ID)=mf_mean1;
Zmp(GRD.ID)=mp_mean1;
Zmd(GRD.ID)=md_mean1;
Zlp(GRD.ID)=lp_mean1;
Zld(GRD.ID)=ld_mean1;
Zb(GRD.ID)=b_mean1;

All1 = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
AllF1 = Zsf+Zmf;
AllP1 = Zsp+Zmp+Zlp;
AllD1 = Zsd+Zmd+Zld;

% middle
Msf=NaN*ones(ni,nj);
Msp=NaN*ones(ni,nj);
Msd=NaN*ones(ni,nj);
Mmf=NaN*ones(ni,nj);
Mmp=NaN*ones(ni,nj);
Mmd=NaN*ones(ni,nj);
Mlp=NaN*ones(ni,nj);
Mld=NaN*ones(ni,nj);
Mb=NaN*ones(ni,nj);

Msf(GRD.ID)=sf_mean2;
Msp(GRD.ID)=sp_mean2;
Msd(GRD.ID)=sd_mean2;
Mmf(GRD.ID)=mf_mean2;
Mmp(GRD.ID)=mp_mean2;
Mmd(GRD.ID)=md_mean2;
Mlp(GRD.ID)=lp_mean2;
Mld(GRD.ID)=ld_mean2;
Mb(GRD.ID)=b_mean2;

All2 = Msp+Msf+Msd+Mmp+Mmf+Mmd+Mlp+Mld;
AllF2 = Msf+Mmf;
AllP2 = Msp+Mmp+Mlp;
AllD2 = Msd+Mmd+Mld;

%end
Esf=NaN*ones(ni,nj);
Esp=NaN*ones(ni,nj);
Esd=NaN*ones(ni,nj);
Emf=NaN*ones(ni,nj);
Emp=NaN*ones(ni,nj);
Emd=NaN*ones(ni,nj);
Elp=NaN*ones(ni,nj);
Eld=NaN*ones(ni,nj);
Eb=NaN*ones(ni,nj);

Esf(GRD.ID)=sf_mean3;
Esp(GRD.ID)=sp_mean3;
Esd(GRD.ID)=sd_mean3;
Emf(GRD.ID)=mf_mean3;
Emp(GRD.ID)=mp_mean3;
Emd(GRD.ID)=md_mean3;
Elp(GRD.ID)=lp_mean3;
Eld(GRD.ID)=ld_mean3;
Eb(GRD.ID)=b_mean3;

All3 = Esp+Esf+Esd+Emp+Emf+Emd+Elp+Eld;
AllF3 = Esf+Emf;
AllP3 = Esp+Emp+Elp;
AllD3 = Esd+Emd+Eld;

%whole yr
Ysf=NaN*ones(ni,nj);
Ysp=NaN*ones(ni,nj);
Ysd=NaN*ones(ni,nj);
Ymf=NaN*ones(ni,nj);
Ymp=NaN*ones(ni,nj);
Ymd=NaN*ones(ni,nj);
Ylp=NaN*ones(ni,nj);
Yld=NaN*ones(ni,nj);
Yb=NaN*ones(ni,nj);

Ysf(GRD.ID)=sf_mean;
Ysp(GRD.ID)=sp_mean;
Ysd(GRD.ID)=sd_mean;
Ymf(GRD.ID)=mf_mean;
Ymp(GRD.ID)=mp_mean;
Ymd(GRD.ID)=md_mean;
Ylp(GRD.ID)=lp_mean;
Yld(GRD.ID)=ld_mean;
Yb(GRD.ID)=b_mean;

All4 = Ysp+Ysf+Ysd+Ymp+Ymf+Ymd+Ylp+Yld;
AllF4 = Ysf+Ymf;
AllP4 = Ysp+Ymp+Ylp;
AllD4 = Ysd+Ymd+Yld;

%% bent
figure(3)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 mean Benthos (g m^-^2) day 1')

%
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Mb))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean Benthos (g m^-^2) day 183')

% 
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Eb))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean Benthos (g m^-^2) day 365')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Yb))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean Benthos (g m^-^2) annual mean')

stamp(mod)
print('-dpng',[ppath mod 'global_BENT.png'])

%% ALL - day 1
figure(4)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllF1))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 mean All F (g m^-^2) Day 1')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllD1))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All D (g m^-^2)')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP1))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All P (g m^-^2)')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All1))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All fishes (g m^-^2)')
stamp(exper)
print('-dpng',[ppath mod 'All_subplot_day1.png'])

%% ALL - day 183
figure(5)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllF2))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 mean All F (g m^-^2) Day 183')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllD2))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All D (g m^-^2)')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP2))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All P (g m^-^2)')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All2))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All fishes (g m^-^2)')
stamp(exper)
print('-dpng',[ppath mod 'All_subplot_day183.png'])

%% ALL - day 365
figure(6)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllF3))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 mean All F (g m^-^2) Day 365')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllD3))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All D (g m^-^2)')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP3))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All P (g m^-^2)')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All3))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All fishes (g m^-^2)')
stamp(exper)
print('-dpng',[ppath mod 'All_subplot_day365.png'])

%% ALL - all year mean
figure(15)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllF4))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 mean All F (g m^-^2) Y1')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllD4))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All D (g m^-^2)')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP4))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All P (g m^-^2)')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All4))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All fishes (g m^-^2)')
stamp(exper)
print('-dpng',[ppath mod 'All_subplot_wholeYR.png'])

%% 8plot by fn type and size
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - sf
subplot('Position',[0.015 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Esf))
cmocean('dense')
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'SF','HorizontalAlignment','center')

%B - 
subplot('Position',[0.015 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Esp))
cmocean('dense')
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'SP','HorizontalAlignment','center')

%C - 
subplot('Position',[0.015 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Esd))
cmocean('dense')
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'SD','HorizontalAlignment','center')

%D - F
subplot('Position',[0.015 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Emf))
cmocean('dense')
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'MF','HorizontalAlignment','center')

%E - P
subplot('Position',[0.47 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Emp))
cmocean('dense')
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'MP','HorizontalAlignment','center')

%F - D
subplot('Position',[0.47 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Emd))
cmocean('dense')
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2])
colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'MD','HorizontalAlignment','center')

%G - B
subplot('Position',[0.47 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Elp))
cmocean('dense')
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'LP','HorizontalAlignment','center')

%H - all
subplot('Position',[0.47 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Eld))
cmocean('dense')
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'LD','HorizontalAlignment','center')
stamp(exper)
%print('-dpng',[ppath mod 'All_stages_subplot_day365.png'])



