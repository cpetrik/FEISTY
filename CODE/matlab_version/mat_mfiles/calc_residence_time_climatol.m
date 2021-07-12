% Calc residence time
% Residence time = 1 / input
% or             = 1 / output
% Total inputs: rec, nu
% Total outputs: gamma, rep, nmort, die (pred), yield (fishing)
% Use per biomass rates (g/g/m2/d)

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'SWhigh_All_fish03';
% cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm000_nmort1_BE08_noCC_RE00100';
% harv = 'SWmarth_MZ50_All_fish03';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/'];

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
ppath = [pp cfile '/Climatol/'];
load('/Volumes/MIP/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_gridspec.mat');

%%
load([fpath 'Means_die_nmort_yield_Climatol_' harv '_' cfile '.mat'],...
    'sf_die','sp_die','sd_die',...
    'mf_die','mp_die','md_die',...
    'lp_die','ld_die','sf_mort','sp_mort','sd_mort',...
    'mf_mort','mp_mort','md_mort',...
    'lp_mort','ld_mort',...
    'mf_yield','mp_yield','md_yield',...
    'lp_yield','ld_yield');

load([fpath 'Means_con_rec_rep_Climatol_' harv '_' cfile '.mat'],...
    'sf_rec','sp_rec','sd_rec',...
    'mf_rec','mp_rec','md_rec',...
    'lp_rec','ld_rec',...
    'mf_rep','lp_rep','ld_rep');

load([fpath 'Means_nu_gam_die_clev_Climatol_' harv '_' cfile '.mat'],...
    'sf_gamma','sp_gamma','sd_gamma',...
    'mf_gamma','mp_gamma','md_gamma','lp_gamma','ld_gamma',...
    'sf_nu','sp_nu','sd_nu',...
    'mf_nu','mp_nu','md_nu','lp_nu','ld_nu');

load([fpath 'Means_Climatol_' harv '_' cfile '.mat'],...
    'sf_bio','sp_bio','sd_bio',...
    'mf_bio','mp_bio','md_bio','lp_bio','ld_bio');

%% add gains and losses
% +nu is g/g/m2/d; rec is g/m2/d
sf_in = max(0,sf_nu) + sf_rec./sf_bio;
sp_in = max(0,sp_nu) + sp_rec./sp_bio;
sd_in = max(0,sd_nu) + sd_rec./sd_bio;
mp_in = max(0,mp_nu) + mp_rec./mf_bio;
md_in = max(0,md_nu) + md_rec./mp_bio;
mf_in = max(0,mf_nu) + mf_rec./md_bio;
lp_in = max(0,lp_nu) + lp_rec./lp_bio;
ld_in = max(0,ld_nu) + ld_rec./ld_bio;

% (-nu,) gamma, rep, nmort, pred, fmort are g/g/m2/d; mort, die, yield are g/m2/d
% sf_out = sf_gamma + sf_mort./sf_bio + sf_die./sf_bio;
% sp_out = sp_gamma + sp_mort./sp_bio + sp_die./sp_bio;
% sd_out = sd_gamma + sd_mort./sd_bio + sd_die./sd_bio;
% mp_out = mp_gamma + mp_mort./mp_bio + mp_die./mp_bio + mp_yield./mp_bio;
% md_out = md_gamma + md_mort./md_bio + md_die./md_bio + md_yield./md_bio;
% mf_out = mf_gamma + mf_rep + mf_mort./mf_bio + mf_die./mf_bio + mf_yield./mf_bio;
% lp_out = lp_gamma + lp_rep + lp_mort./lp_bio + lp_die./lp_bio + lp_yield./lp_bio;
% ld_out = ld_gamma + ld_rep + ld_mort./ld_bio + ld_die./ld_bio + ld_yield./ld_bio;

sf_out = abs(min(0,sf_nu)) + sf_gamma + (0.1/365) + sf_die./sf_bio;
sp_out = abs(min(0,sp_nu)) + sp_gamma + (0.1/365) + sp_die./sp_bio;
sd_out = abs(min(0,sd_nu)) + sd_gamma + (0.1/365) + sd_die./sd_bio;
mp_out = abs(min(0,mf_nu)) + mp_gamma + (0.1/365) + mp_die./mp_bio + (0.1*0.3/365);
md_out = abs(min(0,mp_nu)) + md_gamma + (0.1/365) + md_die./md_bio + (0.1*0.3/365);
mf_out = abs(min(0,md_nu)) + mf_gamma + mf_rep + (0.1/365) + mf_die./mf_bio + (0.3/365);
lp_out = abs(min(0,lp_nu)) + lp_gamma + lp_rep + (0.1/365) + lp_die./lp_bio + (0.3/365);
ld_out = abs(min(0,ld_nu)) + ld_gamma + ld_rep + (0.1/365) + ld_die./ld_bio + (0.3/365);

%% v1
sf_res1 = 1 ./ sf_in;
sp_res1 = 1 ./ sp_in;
sd_res1 = 1 ./ sd_in;
mf_res1 = 1 ./ mf_in;
mp_res1 = 1 ./ mp_in;
md_res1 = 1 ./ md_in;
lp_res1 = 1 ./ lp_in;
ld_res1 = 1 ./ ld_in;

% sf_res1(isinf(sf_res1(:))) = NaN;
% sp_res1(isinf(sp_res1(:))) = NaN;
% sd_res1(isinf(sd_res1(:))) = NaN;
% mf_res1(isinf(mf_res1(:))) = NaN;
% mp_res1(isinf(mp_res1(:))) = NaN;
% md_res1(isinf(md_res1(:))) = NaN;
% lp_res1(isinf(lp_res1(:))) = NaN;
% ld_res1(isinf(ld_res1(:))) = NaN;

%% v2
sf_res2 = 1 ./ sf_out;
sp_res2 = 1 ./ sp_out;
sd_res2 = 1 ./ sd_out;
mf_res2 = 1 ./ mf_out;
mp_res2 = 1 ./ mp_out;
md_res2 = 1 ./ md_out;
lp_res2 = 1 ./ lp_out;
ld_res2 = 1 ./ ld_out;

%% means
sf_mbio = mean(sf_bio,2);
sp_mbio = mean(sp_bio,2);
sd_mbio = mean(sd_bio,2);
mf_mbio = mean(mf_bio,2);
mp_mbio = mean(mp_bio,2);
md_mbio = mean(md_bio,2);
lp_mbio = mean(lp_bio,2);
ld_mbio = mean(ld_bio,2);

sf_min = mean(sf_in,2);
sp_min = mean(sp_in,2);
sd_min = mean(sd_in,2);
mf_min = mean(mf_in,2);
mp_min = mean(mp_in,2);
md_min = mean(md_in,2);
lp_min = mean(lp_in,2);
ld_min = mean(ld_in,2);

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
save([fpath 'Residence_time_means_Climatol_' harv '_' cfile '.mat'],...
  'sf_mbio','sp_mbio','sd_mbio','mf_mbio','mp_mbio','md_mbio','lp_mbio','ld_mbio',...
  'sf_min','sp_min','sd_min','mf_min','mp_min','md_min','lp_min','ld_min',...
  'sf_mout','sp_mout','sd_mout','mf_mout','mp_mout','md_mout','lp_mout','ld_mout',...
  'sf_mres1','sp_mres1','sd_mres1','mf_mres1','mp_mres1','md_mres1','lp_mres1','ld_mres1',...
  'sf_mres2','sp_mres2','sd_mres2','mf_mres2','mp_mres2','md_mres2','lp_mres2','ld_mres2')

%% Histograms
edges = [0:30:360 547 730 912 1095];
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
%bio
Bsf=NaN*ones(ni,nj);
Bsp=NaN*ones(ni,nj);
Bsd=NaN*ones(ni,nj);
Bmf=NaN*ones(ni,nj);
Bmp=NaN*ones(ni,nj);
Bmd=NaN*ones(ni,nj);
Blp=NaN*ones(ni,nj);
Bld=NaN*ones(ni,nj);
%in
Isf=NaN*ones(ni,nj);
Isp=NaN*ones(ni,nj);
Isd=NaN*ones(ni,nj);
Imf=NaN*ones(ni,nj);
Imp=NaN*ones(ni,nj);
Imd=NaN*ones(ni,nj);
Ilp=NaN*ones(ni,nj);
Ild=NaN*ones(ni,nj);
%out
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

%bio
Bsf(ID)=sf_mbio;
Bsp(ID)=sp_mbio;
Bsd(ID)=sd_mbio;
Bmf(ID)=mf_mbio;
Bmp(ID)=mp_mbio;
Bmd(ID)=md_mbio;
Blp(ID)=lp_mbio;
Bld(ID)=ld_mbio;
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

%% 8 plot of bio
f3 = figure('Units','inches','Position',[1 3 6.5 8]);

%A
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Bsf))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'SF bio','HorizontalAlignment','center')

%B
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Bsp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'SP bio','HorizontalAlignment','center')

%C
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Bsd))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'SD bio','HorizontalAlignment','center')

%D
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Bmf))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'MF bio','HorizontalAlignment','center')

%E
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Bmp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'MP bio','HorizontalAlignment','center')

%F
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Bmd))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2])
cb = colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
xlabel(cb,'log_1_0 g m^-^2')
set(gcf,'renderer','painters')
text(0,1.75,'MD bio','HorizontalAlignment','center')

%G
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Blp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'LP bio','HorizontalAlignment','center')

%H
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Bld))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'LD bio','HorizontalAlignment','center')
print('-dpng',[ppath 'Climatol_',harv,'_map_mean_bio_stages.png'])

%% 8 plot of rate in (in)
f4 = figure('Units','inches','Position',[1 3 6.5 8]);

%A
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Isf))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 1])
set(gcf,'renderer','painters')
text(0,1.75,'SF rate in','HorizontalAlignment','center')

%B
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Isp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 1])
set(gcf,'renderer','painters')
text(0,1.75,'SP rate in','HorizontalAlignment','center')

%C
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Isd))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 1])
set(gcf,'renderer','painters')
text(0,1.75,'SD rate in','HorizontalAlignment','center')

%D
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Imf))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 1])
set(gcf,'renderer','painters')
text(0,1.75,'MF rate in','HorizontalAlignment','center')

%E
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Imp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 1])
set(gcf,'renderer','painters')
text(0,1.75,'MP rate in','HorizontalAlignment','center')

%F
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Imd))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 1])
cb = colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
xlabel(cb,'log_1_0 d^-^1')
set(gcf,'renderer','painters')
text(0,1.75,'MD rate in','HorizontalAlignment','center')

%G
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Ilp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 1])
set(gcf,'renderer','painters')
text(0,1.75,'LP rate in','HorizontalAlignment','center')

%H
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Ild))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 1])
set(gcf,'renderer','painters')
text(0,1.75,'LD rate in','HorizontalAlignment','center')
print('-dpng',[ppath 'Climatol_',harv,'_map_mean_rateIn_stages.png'])

%% 8 plot of rate out
f5 = figure('Units','inches','Position',[1 3 6.5 8]);

%A
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Osf))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 1])
set(gcf,'renderer','painters')
text(0,1.75,'SF rate out','HorizontalAlignment','center')

%B
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Osp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 1])
set(gcf,'renderer','painters')
text(0,1.75,'SP rate out','HorizontalAlignment','center')

%C
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Osd))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 1])
set(gcf,'renderer','painters')
text(0,1.75,'SD rate out','HorizontalAlignment','center')

%D
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Omf))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 1])
set(gcf,'renderer','painters')
text(0,1.75,'MF rate out','HorizontalAlignment','center')

%E
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Omp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 1])
set(gcf,'renderer','painters')
text(0,1.75,'MP rate out','HorizontalAlignment','center')

%F
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Omd))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 1])
cb = colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
xlabel(cb,'log_1_0 d^-^1')
set(gcf,'renderer','painters')
text(0,1.75,'MD rate out','HorizontalAlignment','center')

%G
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Olp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 1])
set(gcf,'renderer','painters')
text(0,1.75,'LP rate out','HorizontalAlignment','center')

%H
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Old))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 1])
set(gcf,'renderer','painters')
text(0,1.75,'LD rate out','HorizontalAlignment','center')
print('-dpng',[ppath 'Climatol_',harv,'_map_mean_rateOut_stages.png'])

%% 8 plot of res1
f6 = figure('Units','inches','Position',[1 3 6.5 8]);

%A
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Rsf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 3.5])
set(gcf,'renderer','painters')
text(0,1.75,'SF res (in)','HorizontalAlignment','center')

%B
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Rsp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 3.5])
set(gcf,'renderer','painters')
text(0,1.75,'SP res (in)','HorizontalAlignment','center')

%C
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Rsd))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 3.5])
set(gcf,'renderer','painters')
text(0,1.75,'SD res (in)','HorizontalAlignment','center')

%D
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Rmf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 3.5])
set(gcf,'renderer','painters')
text(0,1.75,'MF res (in)','HorizontalAlignment','center')

%E
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Rmp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 3.5])
set(gcf,'renderer','painters')
text(0,1.75,'MP res (in)','HorizontalAlignment','center')

%F
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Rmd))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 3.5])
cb = colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
xlabel(cb,'log_1_0 d')
set(gcf,'renderer','painters')
text(0,1.75,'MD res (in)','HorizontalAlignment','center')

%G
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Rlp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 3.5])
set(gcf,'renderer','painters')
text(0,1.75,'LP res (in)','HorizontalAlignment','center')

%H
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Rld))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 3.5])
set(gcf,'renderer','painters')
text(0,1.75,'LD res (in)','HorizontalAlignment','center')
print('-dpng',[ppath 'Climatol_',harv,'_map_mean_resIn_stages_log.png'])

%% 8 plot of res2
f7 = figure('Units','inches','Position',[1 3 6.5 8]);

%A
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Ssf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 3.5])
set(gcf,'renderer','painters')
text(0,1.75,'SF res (out)','HorizontalAlignment','center')

%B
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Ssp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 3.5])
set(gcf,'renderer','painters')
text(0,1.75,'SP res (out)','HorizontalAlignment','center')

%C
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Ssd))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 3.5])
set(gcf,'renderer','painters')
text(0,1.75,'SD res (out)','HorizontalAlignment','center')

%D
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Smf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 3.5])
set(gcf,'renderer','painters')
text(0,1.75,'MF res (out)','HorizontalAlignment','center')

%E
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Smp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 3.5])
set(gcf,'renderer','painters')
text(0,1.75,'MP res (out)','HorizontalAlignment','center')

%F
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Smd))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 3.5])
cb = colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
xlabel(cb,'log_1_0 d')
set(gcf,'renderer','painters')
text(0,1.75,'MD res (out)','HorizontalAlignment','center')

%G
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Slp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 3.5])
set(gcf,'renderer','painters')
text(0,1.75,'LP res (out)','HorizontalAlignment','center')

%H
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Sld))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 3.5])
set(gcf,'renderer','painters')
text(0,1.75,'LD res (out)','HorizontalAlignment','center')
print('-dpng',[ppath 'Climatol_',harv,'_map_mean_resOut_stages_log.png'])

%% 8 plot of res1
f9 = figure('Units','inches','Position',[1 3 6.5 8]);

%A
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 350])
set(gcf,'renderer','painters')
text(0,1.75,'SF res (in)','HorizontalAlignment','center')

%B
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 350])
set(gcf,'renderer','painters')
text(0,1.75,'SP res (in)','HorizontalAlignment','center')

%C
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsd))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 350])
set(gcf,'renderer','painters')
text(0,1.75,'SD res (in)','HorizontalAlignment','center')

%D
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rmf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 350])
set(gcf,'renderer','painters')
text(0,1.75,'MF res (in)','HorizontalAlignment','center')

%E
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rmp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 350])
set(gcf,'renderer','painters')
text(0,1.75,'MP res (in)','HorizontalAlignment','center')

%F
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rmd))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 350])
cb = colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
xlabel(cb,'d')
set(gcf,'renderer','painters')
text(0,1.75,'MD res (in)','HorizontalAlignment','center')

%G
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rlp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 350])
set(gcf,'renderer','painters')
text(0,1.75,'LP res (in)','HorizontalAlignment','center')

%H
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rld))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 350])
set(gcf,'renderer','painters')
text(0,1.75,'LD res (in)','HorizontalAlignment','center')
print('-dpng',[ppath 'Climatol_',harv,'_map_mean_resIn_stages.png'])

%% 8 plot of res2
f10 = figure('Units','inches','Position',[1 3 6.5 8]);

%A
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Ssf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 350])
set(gcf,'renderer','painters')
text(0,1.75,'SF res (out)','HorizontalAlignment','center')

%B
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Ssp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 350])
set(gcf,'renderer','painters')
text(0,1.75,'SP res (out)','HorizontalAlignment','center')

%C
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Ssd))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 350])
set(gcf,'renderer','painters')
text(0,1.75,'SD res (out)','HorizontalAlignment','center')

%D
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Smf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 350])
set(gcf,'renderer','painters')
text(0,1.75,'MF res (out)','HorizontalAlignment','center')

%E
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Smp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 350])
set(gcf,'renderer','painters')
text(0,1.75,'MP res (out)','HorizontalAlignment','center')

%F
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Smd))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 350])
cb = colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
xlabel(cb,'d')
set(gcf,'renderer','painters')
text(0,1.75,'MD res (out)','HorizontalAlignment','center')

%G
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Slp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 350])
set(gcf,'renderer','painters')
text(0,1.75,'LP res (out)','HorizontalAlignment','center')

%H
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Sld))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 350])
set(gcf,'renderer','painters')
text(0,1.75,'LD res (out)','HorizontalAlignment','center')
print('-dpng',[ppath 'Climatol_',harv,'_map_mean_resOut_stages.png'])

%% Just small
f11 = figure('Units','inches','Position',[1 3 6.5 7.25]);

%A
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 60])
set(gcf,'renderer','painters')
text(0,1.75,'SF res (in)','HorizontalAlignment','center')

%B
subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Ssf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 60])
set(gcf,'renderer','painters')
text(0,1.75,'SF res (out)','HorizontalAlignment','center')

%C
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 60])
set(gcf,'renderer','painters')
text(0,1.75,'SP res (in)','HorizontalAlignment','center')

%D
subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Ssp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 60])
set(gcf,'renderer','painters')
text(0,1.75,'SP res (out)','HorizontalAlignment','center')

%E
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsd))
cmocean('speed')
cb = colorbar('Position',[0.85 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
xlabel(cb,'days')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 60])
set(gcf,'renderer','painters')
text(0,1.75,'SD res (in)','HorizontalAlignment','center')

%F
subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Ssd))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 60])
set(gcf,'renderer','painters')
text(0,1.75,'SD res (out)','HorizontalAlignment','center')
print('-dpng',[ppath 'Climatol_',harv,'_map_mean_resInOut_small.png'])

%% test how low in tropics/subtropics
f12 = figure('Units','inches','Position',[1 3 6.5 7.25]);

%A
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([10 30])
set(gcf,'renderer','painters')
text(0,1.75,'SF res (in)','HorizontalAlignment','center')

%B
subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Ssf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([10 30])
set(gcf,'renderer','painters')
text(0,1.75,'SF res (out)','HorizontalAlignment','center')

%C
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([10 30])
set(gcf,'renderer','painters')
text(0,1.75,'SP res (in)','HorizontalAlignment','center')

%D
subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Ssp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([10 30])
set(gcf,'renderer','painters')
text(0,1.75,'SP res (out)','HorizontalAlignment','center')

%E
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsd))
cmocean('speed')
cb = colorbar('Position',[0.85 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
xlabel(cb,'days')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([10 30])
set(gcf,'renderer','painters')
text(0,1.75,'SD res (in)','HorizontalAlignment','center')

%F
subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Ssd))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([10 30])
set(gcf,'renderer','painters')
text(0,1.75,'SD res (out)','HorizontalAlignment','center')

%%
quan(1,:) = quantile(Ssf(:),[0.1 0.25 0.5 0.75 0.9]);
quan(2,:) = quantile(Ssp(:),[0.1 0.25 0.5 0.75 0.9]);
quan(3,:) = quantile(Ssd(:),[0.1 0.25 0.5 0.75 0.9]);
quan(4,:) = quantile(Smf(:),[0.1 0.25 0.5 0.75 0.9]);
quan(5,:) = quantile(Smp(:),[0.1 0.25 0.5 0.75 0.9]);
quan(6,:) = quantile(Smd(:),[0.1 0.25 0.5 0.75 0.9]);
quan(7,:) = quantile(Slp(:),[0.1 0.25 0.5 0.75 0.9]);
quan(8,:) = quantile(Sld(:),[0.1 0.25 0.5 0.75 0.9]);

TQ = array2table(quan,'VariableNames',{'0.10','0.25','0.5','0.75','0.90'},...
    'RowNames',{'SF','SP','SD','MF','MP','MD','LP','LD'})

save([fpath 'Residence_time_means_Climatol_' harv '_' cfile '.mat'],...
    'quan','TQ','-append');
writetable(TQ,[fpath 'Residence_time_quant_Climatol_' harv '_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);



