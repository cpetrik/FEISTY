% Calc survival rates
% nmort, die (pred), yield (fishing)
% nu - gamma?

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/'];

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
ppath = [pp cfile '/Climatol/'];
load('/Volumes/FEISTY/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_gridspec.mat');

%%
load([fpath 'Means_die_nmort_yield_Climatol_' harv '_' cfile '.mat'],...
    'sf_die','sp_die','sd_die',...
    'mf_die','mp_die','md_die',...
    'lp_die','ld_die','sf_mort','sp_mort','sd_mort',...
    'mf_mort','mp_mort','md_mort',...
    'lp_mort','ld_mort',...
    'mf_yield','mp_yield','md_yield',...
    'lp_yield','ld_yield');

load([fpath 'Means_nu_gam_die_clev_Climatol_' harv '_' cfile '.mat'],...
    'sf_gamma','sp_gamma','sd_gamma',...
    'mf_gamma','mp_gamma','md_gamma','lp_gamma','ld_gamma',...
    'sf_nu','sp_nu','sd_nu',...
    'mf_nu','mp_nu','md_nu','lp_nu','ld_nu');

load([fpath 'Means_Climatol_' harv '_' cfile '.mat'],...
    'sf_bio','sp_bio','sd_bio',...
    'mf_bio','mp_bio','md_bio','lp_bio','ld_bio');

load([fpath 'Means_con_rec_rep_Climatol_' harv '_' cfile '.mat'],...
    'mf_rep','lp_rep','ld_rep');

%% add losses
sf_out = sf_nu - sf_gamma;
sp_out = sp_nu - sp_gamma;
sd_out = sd_nu - sd_gamma;
mp_out = mf_nu - mp_gamma;
md_out = mp_nu - md_gamma;
mf_out = md_nu - mf_gamma - mf_rep;
lp_out = lp_nu - lp_gamma - lp_rep;
ld_out = ld_nu - ld_gamma - ld_rep;

sf_tmort = (0.1/365) + sf_die./sf_bio;
sp_tmort = (0.1/365) + sp_die./sp_bio;
sd_tmort = (0.1/365) + sd_die./sd_bio;
mp_tmort = (0.1/365) + mp_die./mp_bio + (0.1*0.3/365);
md_tmort = (0.1/365) + md_die./md_bio + (0.1*0.3/365);
mf_tmort = (0.1/365) + mf_die./mf_bio + (0.3/365);
lp_tmort = (0.1/365) + lp_die./lp_bio + (0.3/365);
ld_tmort = (0.1/365) + ld_die./ld_bio + (0.3/365);

%% v1
sf_res1 = 1 - sf_tmort;
sp_res1 = 1 - sp_tmort;
sd_res1 = 1 - sd_tmort;
mf_res1 = 1 - mf_tmort;
mp_res1 = 1 - mp_tmort;
md_res1 = 1 - md_tmort;
lp_res1 = 1 - lp_tmort;
ld_res1 = 1 - ld_tmort;

%% v2
sf_res2 = 1 - sf_out - sf_tmort;
sp_res2 = 1 - sp_out - sp_tmort;
sd_res2 = 1 - sd_out - sd_tmort;
mf_res2 = 1 - mf_out - mf_tmort;
mp_res2 = 1 - mp_out - mp_tmort;
md_res2 = 1 - md_out - md_tmort;
lp_res2 = 1 - lp_out - lp_tmort;
ld_res2 = 1 - ld_out - ld_tmort;

%% means
sf_min = mean(sf_tmort,2);
sp_min = mean(sp_tmort,2);
sd_min = mean(sd_tmort,2);
mf_min = mean(mf_tmort,2);
mp_min = mean(mp_tmort,2);
md_min = mean(md_tmort,2);
lp_min = mean(lp_tmort,2);
ld_min = mean(ld_tmort,2);

sf_mout = mean(sf_out,2);
sp_mout = mean(sp_out,2);
sd_mout = mean(sd_out,2);
mf_mout = mean(mf_out,2);
mp_mout = mean(mp_out,2);
md_mout = mean(md_out,2);
lp_mout = mean(lp_out,2);
ld_mout = mean(ld_out,2);

sf_mres1 = nanmean(sf_res1,2);
sp_mres1 = nanmean(sp_res1,2);
sd_mres1 = nanmean(sd_res1,2);
mf_mres1 = nanmean(mf_res1,2);
mp_mres1 = nanmean(mp_res1,2);
md_mres1 = nanmean(md_res1,2);
lp_mres1 = nanmean(lp_res1,2);
ld_mres1 = nanmean(ld_res1,2);

sf_mres2 = mean(sf_res2,2);
sp_mres2 = mean(sp_res2,2);
sd_mres2 = mean(sd_res2,2);
mf_mres2 = mean(mf_res2,2);
mp_mres2 = mean(mp_res2,2);
md_mres2 = mean(md_res2,2);
lp_mres2 = mean(lp_res2,2);
ld_mres2 = mean(ld_res2,2);

%% Save
save([fpath 'Mort_surv_means_Climatol_' harv '_' cfile '.mat'],...
  'sf_min','sp_min','sd_min','mf_min','mp_min','md_min','lp_min','ld_min',...
  'sf_mout','sp_mout','sd_mout','mf_mout','mp_mout','md_mout','lp_mout','ld_mout',...
  'sf_mres1','sp_mres1','sd_mres1','mf_mres1','mp_mres1','md_mres1','lp_mres1','ld_mres1',...
  'sf_mres2','sp_mres2','sd_mres2','mf_mres2','mp_mres2','md_mres2','lp_mres2','ld_mres2')

%% Histograms
edges = [0.5:0.05:1];
figure(1)
subplot(3,3,1)
histogram((sf_mres1),edges)
title('SF')

subplot(3,3,2)
histogram((sp_mres1),edges)
title('SP')

subplot(3,3,3)
histogram((sd_mres1),edges)
title('SD')

subplot(3,3,4)
histogram((mf_mres1),edges)
title('MF')

subplot(3,3,5)
histogram((mp_mres1),edges)
title('MP')

subplot(3,3,6)
histogram((md_mres1),edges)
title('MD')

subplot(3,3,8)
histogram((lp_mres1),edges)
title('LP')

subplot(3,3,9)
histogram((ld_mres1),edges)
title('LD')

figure(2)
subplot(3,3,1)
histogram((sf_mres2),edges)
title('SF')

subplot(3,3,2)
histogram((sp_mres2),edges)
title('SP')

subplot(3,3,3)
histogram((sd_mres2),edges)
title('SD')

subplot(3,3,4)
histogram((mf_mres2),edges)
title('MF')

subplot(3,3,5)
histogram((mp_mres2),edges)
title('MP')

subplot(3,3,6)
histogram((md_mres2),edges)
title('MD')

subplot(3,3,8)
histogram((lp_mres2),edges)
title('LP')

subplot(3,3,9)
histogram((ld_mres2),edges)
title('LD')

%% Maps
% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90;
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

%% Plots in space
%tmort
Isf=NaN*ones(ni,nj);
Isp=NaN*ones(ni,nj);
Isd=NaN*ones(ni,nj);
Imf=NaN*ones(ni,nj);
Imp=NaN*ones(ni,nj);
Imd=NaN*ones(ni,nj);
Ilp=NaN*ones(ni,nj);
Ild=NaN*ones(ni,nj);
%tmort + nu diff
Osf=NaN*ones(ni,nj);
Osp=NaN*ones(ni,nj);
Osd=NaN*ones(ni,nj);
Omf=NaN*ones(ni,nj);
Omp=NaN*ones(ni,nj);
Omd=NaN*ones(ni,nj);
Olp=NaN*ones(ni,nj);
Old=NaN*ones(ni,nj);
%res1
Rsf=NaN*ones(ni,nj);
Rsp=NaN*ones(ni,nj);
Rsd=NaN*ones(ni,nj);
Rmf=NaN*ones(ni,nj);
Rmp=NaN*ones(ni,nj);
Rmd=NaN*ones(ni,nj);
Rlp=NaN*ones(ni,nj);
Rld=NaN*ones(ni,nj);
%res2
Ssf=NaN*ones(ni,nj);
Ssp=NaN*ones(ni,nj);
Ssd=NaN*ones(ni,nj);
Smf=NaN*ones(ni,nj);
Smp=NaN*ones(ni,nj);
Smd=NaN*ones(ni,nj);
Slp=NaN*ones(ni,nj);
Sld=NaN*ones(ni,nj);

%in
Isf(ID)=sf_min;
Isp(ID)=sp_min;
Isd(ID)=sd_min;
Imf(ID)=mf_min;
Imp(ID)=mp_min;
Imd(ID)=md_min;
Ilp(ID)=lp_min;
Ild(ID)=ld_min;
%out
Osf(ID)=sf_mout;
Osp(ID)=sp_mout;
Osd(ID)=sd_mout;
Omf(ID)=mf_mout;
Omp(ID)=mp_mout;
Omd(ID)=md_mout;
Olp(ID)=lp_mout;
Old(ID)=ld_mout;
%res1
Rsf(ID)=sf_mres1;
Rsp(ID)=sp_mres1;
Rsd(ID)=sd_mres1;
Rmf(ID)=mf_mres1;
Rmp(ID)=mp_mres1;
Rmd(ID)=md_mres1;
Rlp(ID)=lp_mres1;
Rld(ID)=ld_mres1;
%res2
Ssf(ID)=sf_mres2;
Ssp(ID)=sp_mres2;
Ssd(ID)=sd_mres2;
Smf(ID)=mf_mres2;
Smp(ID)=mp_mres2;
Smd(ID)=md_mres2;
Slp(ID)=lp_mres2;
Sld(ID)=ld_mres2;

%% 8 plot of tmort
f4 = figure('Units','inches','Position',[1 3 6.5 8]);

%A
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Isf))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.25])
set(gcf,'renderer','painters')
text(0,1.75,'SF tmort','HorizontalAlignment','center')

%B
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Isp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.25])
set(gcf,'renderer','painters')
text(0,1.75,'SP tmort','HorizontalAlignment','center')

%C
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Isd))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.25])
set(gcf,'renderer','painters')
text(0,1.75,'SD tmort','HorizontalAlignment','center')

%D
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Imf))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.25])
set(gcf,'renderer','painters')
text(0,1.75,'MF tmort','HorizontalAlignment','center')

%E
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Imp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.25])
set(gcf,'renderer','painters')
text(0,1.75,'MP tmort','HorizontalAlignment','center')

%F
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Imd))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.25])
cb = colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
xlabel(cb,'d^-^1')
set(gcf,'renderer','painters')
text(0,1.75,'MD tmort','HorizontalAlignment','center')

%G
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Ilp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.25])
set(gcf,'renderer','painters')
text(0,1.75,'LP tmort','HorizontalAlignment','center')

%H
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Ild))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.25])
set(gcf,'renderer','painters')
text(0,1.75,'LD tmort','HorizontalAlignment','center')
print('-dpng',[ppath 'Climatol_map_mean_tmort_stages.png'])

%% 8 plot of tmort + nu diff
f5 = figure('Units','inches','Position',[1 3 6.5 8]);

%A
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Osf))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.25])
set(gcf,'renderer','painters')
text(0,1.75,'SF tmort + dnu','HorizontalAlignment','center')

%B
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Osp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.25])
set(gcf,'renderer','painters')
text(0,1.75,'SP tmort + dnu','HorizontalAlignment','center')

%C
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Osd))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.25])
set(gcf,'renderer','painters')
text(0,1.75,'SD tmort + dnu','HorizontalAlignment','center')

%D
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Omf))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.25])
set(gcf,'renderer','painters')
text(0,1.75,'MF tmort + dnu','HorizontalAlignment','center')

%E
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Omp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.25])
set(gcf,'renderer','painters')
text(0,1.75,'MP tmort + dnu','HorizontalAlignment','center')

%F
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Omd))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.25])
cb = colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
xlabel(cb,'d^-^1')
set(gcf,'renderer','painters')
text(0,1.75,'MD tmort + dnu','HorizontalAlignment','center')

%G
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Olp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.25])
set(gcf,'renderer','painters')
text(0,1.75,'LP tmort + dnu','HorizontalAlignment','center')

%H
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Old))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.25])
set(gcf,'renderer','painters')
text(0,1.75,'LD tmort + dnu','HorizontalAlignment','center')
print('-dpng',[ppath 'Climatol_map_mean_tmort_dnu_stages.png'])

%% 8 plot of res1
f9 = figure('Units','inches','Position',[1 3 6.5 8]);

%A
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'SF surv1','HorizontalAlignment','center')

%B
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'SP surv1','HorizontalAlignment','center')

%C
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsd))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'SD surv1','HorizontalAlignment','center')

%D
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rmf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'MF surv1','HorizontalAlignment','center')

%E
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rmp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'MP surv1','HorizontalAlignment','center')

%F
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rmd))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
cb = colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
xlabel(cb,'d')
set(gcf,'renderer','painters')
text(0,1.75,'MD surv1','HorizontalAlignment','center')

%G
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rlp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'LP surv1','HorizontalAlignment','center')

%H
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rld))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'LD surv1','HorizontalAlignment','center')
print('-dpng',[ppath 'Climatol_map_mean_surv1_stages.png'])

%% 8 plot of res2
f10 = figure('Units','inches','Position',[1 3 6.5 8]);

%A
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Ssf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'SF surv2','HorizontalAlignment','center')

%B
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Ssp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'SP surv2','HorizontalAlignment','center')

%C
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Ssd))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'SD surv2','HorizontalAlignment','center')

%D
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Smf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'MF surv2','HorizontalAlignment','center')

%E
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Smp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'MP surv2','HorizontalAlignment','center')

%F
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Smd))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
cb = colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
xlabel(cb,'d')
set(gcf,'renderer','painters')
text(0,1.75,'MD surv2','HorizontalAlignment','center')

%G
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Slp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'LP surv2','HorizontalAlignment','center')

%H
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Sld))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'LD surv2','HorizontalAlignment','center')
print('-dpng',[ppath 'Climatol_map_mean_surv2_stages.png'])

%% Just small
f11 = figure('Units','inches','Position',[1 3 6.5 7.25]);

%A
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'SF surv1','HorizontalAlignment','center')

%B
subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Ssf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'SF surv2','HorizontalAlignment','center')

%C
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'SP surv1','HorizontalAlignment','center')

%D
subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Ssp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'SP surv2','HorizontalAlignment','center')

%E
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsd))
cmocean('speed')
cb = colorbar('Position',[0.85 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
xlabel(cb,'days')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'SD surv1','HorizontalAlignment','center')

%F
subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Ssd))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.75 1])
set(gcf,'renderer','painters')
text(0,1.75,'SD surv2','HorizontalAlignment','center')
print('-dpng',[ppath 'Climatol_map_mean_surv_small.png'])
