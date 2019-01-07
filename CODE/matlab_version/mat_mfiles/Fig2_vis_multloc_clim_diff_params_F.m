% Visualize output of POEM
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/v4_kt85_BE75/';
cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_area_1deg.mat']);
AREA_OCN = max(area,1);

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

cfileA = 'Dc_enc70-b200_m4-b250-k086_c20-b250_D075_J100_A100_Sm025_nmort1_BE08_noCC_RE00100';
fpathA=['/Volumes/GFDL/NC/Matlab_new_size/' cfileA '/'];
load([fpathA 'Means_bio_prod_fish_Climatol_' harv '_' cfileA '.mat'],...
    'sf_mean','mf_mean','sp_mean','mp_mean','lp_mean');
Asf=NaN*ones(ni,nj);
Amf=NaN*ones(ni,nj);
Asf(ID)=sf_mean;
Amf(ID)=mf_mean;
Asp=NaN*ones(ni,nj);
Amp=NaN*ones(ni,nj);
Alp=NaN*ones(ni,nj);
Asp(ID)=sp_mean;
Amp(ID)=mp_mean;
Alp(ID)=lp_mean;
clear sf_mean mf_mean sp_mean mp_mean lp_mean

cfileB = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A100_Sm025_nmort1_BE08_noCC_RE00100';
fpathB=['/Volumes/GFDL/NC/Matlab_new_size/' cfileB '/'];
load([fpathB 'Means_bio_prod_fish_Climatol_' harv '_' cfileB '.mat'],...
    'sf_mean','mf_mean','sp_mean','mp_mean','lp_mean');
Bsf=NaN*ones(ni,nj);
Bmf=NaN*ones(ni,nj);
Bsf(ID)=sf_mean;
Bmf(ID)=mf_mean;
Bsp=NaN*ones(ni,nj);
Bmp=NaN*ones(ni,nj);
Blp=NaN*ones(ni,nj);
Bsp(ID)=sp_mean;
Bmp(ID)=mp_mean;
Blp(ID)=lp_mean;
clear sf_mean mf_mean sp_mean mp_mean lp_mean

cfileC = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpathC=['/Volumes/GFDL/NC/Matlab_new_size/' cfileC '/'];
ppath = [pp cfileC '/'];
load([fpathC 'Means_bio_prod_fish_Climatol_' harv '_' cfileC '.mat'],...
    'sf_mean','mf_mean','sp_mean','mp_mean','lp_mean');
Csf=NaN*ones(ni,nj);
Cmf=NaN*ones(ni,nj);
Csf(ID)=sf_mean;
Cmf(ID)=mf_mean;
Csp=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Csp(ID)=sp_mean;
Cmp(ID)=mp_mean;
Clp(ID)=lp_mean;
clear sf_mean mf_mean sp_mean mp_mean lp_mean

% Diff maps of all fish
AF = Asf+Amf;
BF = Bsf+Bmf;
CF = Csf+Cmf;
AP = Asp+Amp+Alp;
BP = Bsp+Bmp+Blp;
CP = Csp+Cmp+Clp;

%% Ratios on subplots red-white-blue
figure(1)
%A
subplot('Position',[0.1 0.68 0.8 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'A')

%B
subplot('Position',[0.1 0.37 0.8 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(BF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'B')

%C
subplot('Position',[0.1 0.06 0.8 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(CF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.33 0.04 0.34 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
text(-2.75,1.75,'C')
%stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_Fcomp_params.png'])

%%
figure(2)
%A
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'A')
text(-1.0,1.55,'Forage Fishes')

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-1.0,1.55,'Large Pelagics')

%B
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(BF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'B')

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(BP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')

%C
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(CF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'C')
%stamp([harv '_' cfile])

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(CP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.25 0.04 0.34 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
print('-dpng',[ppath 'Climatol_' harv '_FPcomp_params.png'])

%% Differences in biomass
AFtot = AF .* AREA_OCN;
BFtot = BF .* AREA_OCN;
CFtot = CF .* AREA_OCN;
APtot = AP .* AREA_OCN;
BPtot = BP .* AREA_OCN;
CPtot = CP .* AREA_OCN;

tot_bioAF = nansum(AF(:)) * 1e-6;
tot_bioBF = nansum(BF(:)) * 1e-6;
tot_bioCF = nansum(CF(:)) * 1e-6;
tot_bioAP = nansum(AP(:)) * 1e-6;
tot_bioBP = nansum(BP(:)) * 1e-6;
tot_bioCP = nansum(CP(:)) * 1e-6;

diffABF=(tot_bioAF-tot_bioBF)./tot_bioAF % 3930%
diffABP=(tot_bioAP-tot_bioBP)./tot_bioAP % 38.4%

diffBCF=(tot_bioBF-tot_bioCF)./tot_bioBF % 202%
diffBCP=(tot_bioBP-tot_bioCP)./tot_bioBP % 35.2%

multABF=(tot_bioBF)./tot_bioAF % 40.3034x greater
multABP=(tot_bioBP)./tot_bioAP % 0.6156x = 1/1.6156 = 1.6243x less

multBCF=(tot_bioCF)./tot_bioBF % 3.0194x
multBCP=(tot_bioCP)./tot_bioBP % 0.6480x = 1/1.6480 = 1.5433x less

tab(1,1) = tot_bioAF;
tab(2,1) = tot_bioBF;
tab(3,1) = tot_bioCF;
tab(1,2) = tot_bioAP;
tab(2,2) = tot_bioBP; 
tab(3,2) = tot_bioCP;
Tab=array2table(tab,'VariableNames',{'F','P'},...
    'RowNames',{'A','B','C'});
writetable(Tab,[spath 'Tot_biom_FPcomp_params.csv'],'Delimiter',',','WriteRowNames',true);
save([fpathC 'Tot_biom_FPcomp_params_Climatol_',harv,'_' cfile '.mat'],'tab','Tab');




