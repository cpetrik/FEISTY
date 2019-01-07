% Visualize output of POEM
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
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
cmBP=cbrewer('seq','BuPu',50);

%% Plots in space
ptext = {'base','h25','h100','gam25','gam100','amet2','amet8','lam056','lam084',...
    'bc100','bc320','be100','be320','bm100','bm320','BE0375','BE15','RE0005',...
    'RE002','fish015','fish06','kc0302','kc1208','ke0302','ke1208','kt0302',...
    'kt1208','kap25','kap75','unat05','unat20','A025','A100','Sm0125',...
    'Sm050','J025','J100','D025','D100'};

%GFDL/NC/Matlab_new_size/Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100/param_sens/
cfile = 'Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100';
% sfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,...
%         '/param_sens/Climatol_All_fish03_',pt,'.mat'];

cfileA = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,...
        '/param_sens/Climatol_All_fish03_kt0302.mat'];
load(cfileA,'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Afish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Aall=NaN*ones(ni,nj);
Aall(ID)=Afish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

cfileB = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,...
        '/param_sens/Climatol_All_fish03_base.mat'];
load(cfileB,'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Bfish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Ball=NaN*ones(ni,nj);
Ball(ID)=Bfish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

cfileC = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,...
        '/param_sens/Climatol_All_fish03_kt1208.mat'];
load(cfileC,'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Cfish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Call=NaN*ones(ni,nj);
Call(ID)=Cfish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

%%
figure(1)
%A
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Aall))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'A')
text(-0.75,1.75,'k_M=0.0302')
%B
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Ball))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'B')
text(-0.75,1.75,'k_M=0.0604')
%C
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Call))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'C')
text(-0.75,1.75,'k_M=0.1208')
%stamp('')
colorbar('Position',[0.04 0.04 0.339 0.025],'orientation','horizontal')
print('-dpng',[pp 'Fig3_Climatol_' harv '_All_comp_kt_param_sensv2_BP.png'])

