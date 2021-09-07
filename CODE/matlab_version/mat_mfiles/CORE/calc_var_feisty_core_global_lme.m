% CORE forced FEISTY runs
% calc interann variability by grid cell and lme

clear all
close all

%% Fish data
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/CORE/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

load([fpath 'Means_core_fished_' cfile '.mat'],...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'b_mean','lp_mean','ld_mean');

load('/Volumes/MIP/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');
load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

[ni,nj]=size(geolon_t);
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);

%%
[nid,nt]=size(sf_mean);

%% Biomass on grid
Zsf=NaN*ones(ni*nj,nt);
Zsp=NaN*ones(ni*nj,nt);
Zsd=NaN*ones(ni*nj,nt);
Zmf=NaN*ones(ni*nj,nt);
Zmp=NaN*ones(ni*nj,nt);
Zmd=NaN*ones(ni*nj,nt);
Zlp=NaN*ones(ni*nj,nt);
Zld=NaN*ones(ni*nj,nt);
Zb=NaN*ones(ni*nj,nt);

Zsf(GRD.ID,:)=sf_mean;
Zsp(GRD.ID,:)=sp_mean;
Zsd(GRD.ID,:)=sd_mean;
Zmf(GRD.ID,:)=mf_mean;
Zmp(GRD.ID,:)=mp_mean;
Zmd(GRD.ID,:)=md_mean;
Zlp(GRD.ID,:)=lp_mean;
Zld(GRD.ID,:)=ld_mean;
Zb(GRD.ID,:)=b_mean;

All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
AllF = Zsf+Zmf;
AllP = Zsp+Zmp+Zlp;
AllD = Zsd+Zmd+Zld;
AllS = Zsp+Zsf+Zsd;
AllM = Zmp+Zmf+Zmd;
AllL = Zlp+Zld;
FracPD = AllP ./ (AllP+AllD);
FracPF = AllP ./ (AllP+AllF);
FracLM = AllL ./ (AllL+AllM);

%% anomalies
asf = Zsf - nanmean(Zsf,2);
asp = Zsp - nanmean(Zsp,2);
asd = Zsd - nanmean(Zsd,2);
amf = Zmf - nanmean(Zmf,2);
amp = Zmp - nanmean(Zmp,2);
amd = Zmd - nanmean(Zmd,2);
alp = Zlp - nanmean(Zlp,2);
ald = Zld - nanmean(Zld,2);
ab = Zb - nanmean(Zb,2);
aa = All - nanmean(All,2);
as = AllS - nanmean(AllS,2);
am = AllM - nanmean(AllM,2);
al = AllL - nanmean(AllL,2);
af = AllF - nanmean(AllF,2);
ap = AllP - nanmean(AllP,2);
ad = AllD - nanmean(AllD,2);
apf = FracPF - nanmean(FracPF,2);
apd = FracPD - nanmean(FracPD,2);
alm = FracLM - nanmean(FracLM,2);

%% var by grid cell
vsf = var(asf,0,2,'omitnan');
vsp = var(asp,0,2,'omitnan');
vsd = var(asd,0,2,'omitnan');
vmf = var(amf,0,2,'omitnan');
vmp = var(amp,0,2,'omitnan');
vmd = var(amd,0,2,'omitnan');
vlp = var(alp,0,2,'omitnan');
vld = var(ald,0,2,'omitnan');
vb = var(ab,0,2,'omitnan');
va = var(aa,0,2,'omitnan');
vs = var(as,0,2,'omitnan');
vm = var(am,0,2,'omitnan');
vl = var(al,0,2,'omitnan');
vf = var(af,0,2,'omitnan');
vp = var(ap,0,2,'omitnan');
vd = var(ad,0,2,'omitnan');
vpf = var(apf,0,2,'omitnan');
vpd = var(apd,0,2,'omitnan');
vlm = var(alm,0,2,'omitnan');

%% var by lme
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
load([cpath 'lme_mask_esm2m.mat']);
tlme = lme_mask_esm2m';

lme_sf_var1 = NaN*ones(66,1);
lme_sp_var1 = NaN*ones(66,1);
lme_sd_var1 = NaN*ones(66,1);
lme_mf_var1 = NaN*ones(66,1);
lme_mp_var1 = NaN*ones(66,1);
lme_md_var1 = NaN*ones(66,1);
lme_lp_var1 = NaN*ones(66,1);
lme_ld_var1 = NaN*ones(66,1);
lme_b_var1 = NaN*ones(66,1);
lme_a_var1 = NaN*ones(66,1);
lme_s_var1 = NaN*ones(66,1);
lme_m_var1 = NaN*ones(66,1);
lme_l_var1 = NaN*ones(66,1);
lme_f_var1 = NaN*ones(66,1);
lme_p_var1 = NaN*ones(66,1);
lme_d_var1 = NaN*ones(66,1);
lme_pf_var1 = NaN*ones(66,1);
lme_pd_var1 = NaN*ones(66,1);
lme_lm_var1 = NaN*ones(66,1);

for L=1:66
    lid = find(tlme==L);
    
    lme_sf_var1(L,1) = nanmean(var(asf(lid,:),0,2,'omitnan'));
    lme_sp_var1(L,1) = nanmean(var(asp(lid,:),0,2,'omitnan'));
    lme_sd_var1(L,1) = nanmean(var(asd(lid,:),0,2,'omitnan'));
    lme_mf_var1(L,1) = nanmean(var(amf(lid,:),0,2,'omitnan'));
    lme_mp_var1(L,1) = nanmean(var(amp(lid,:),0,2,'omitnan'));
    lme_md_var1(L,1) = nanmean(var(amd(lid,:),0,2,'omitnan'));
    lme_lp_var1(L,1) = nanmean(var(alp(lid,:),0,2,'omitnan'));
    lme_ld_var1(L,1) = nanmean(var(ald(lid,:),0,2,'omitnan'));
    lme_b_var1(L,1) = nanmean(var(ab(lid,:),0,2,'omitnan'));
    lme_a_var1(L,1) = nanmean(var(aa(lid,:),0,2,'omitnan'));
    lme_s_var1(L,1) = nanmean(var(as(lid,:),0,2,'omitnan'));
    lme_m_var1(L,1) = nanmean(var(am(lid,:),0,2,'omitnan'));
    lme_l_var1(L,1) = nanmean(var(al(lid,:),0,2,'omitnan'));
    lme_f_var1(L,1) = nanmean(var(af(lid,:),0,2,'omitnan'));
    lme_p_var1(L,1) = nanmean(var(ap(lid,:),0,2,'omitnan'));
    lme_d_var1(L,1) = nanmean(var(ad(lid,:),0,2,'omitnan'));
    lme_pf_var1(L,1) = nanmean(var(apf(lid,:),0,2,'omitnan'));
    lme_pd_var1(L,1) = nanmean(var(apd(lid,:),0,2,'omitnan'));
    lme_lm_var1(L,1) = nanmean(var(alm(lid,:),0,2,'omitnan'));
    
end

%% map info
latlim=[-90 90];
lonlim=[-280 80];
load coastlines

cmYR=cbrewer('seq','YlOrRd',66,'pchip');
cmYOB=cbrewer('seq','YlOrBr',66,'pchip');
cmOR=cbrewer('seq','OrRd',66,'pchip');

%% reshape
vsf = reshape(vsf,ni,nj);
vsp = reshape(vsp,ni,nj);
vsd = reshape(vsd,ni,nj);
vmf = reshape(vmf,ni,nj);
vmp = reshape(vmp,ni,nj);
vmd = reshape(vmd,ni,nj);
vlp = reshape(vlp,ni,nj);
vld = reshape(vld,ni,nj);
vb = reshape(vb,ni,nj);
va = reshape(va,ni,nj);
vs = reshape(vs,ni,nj);
vm = reshape(vm,ni,nj);
vl = reshape(vl,ni,nj);
vf = reshape(vf,ni,nj);
vp = reshape(vp,ni,nj);
vd = reshape(vd,ni,nj);
vpf = reshape(vpf,ni,nj);
vpd = reshape(vpd,ni,nj);
vlm = reshape(vlm,ni,nj);

%% map
% 8plot by stage 
f1 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - sf
subplot('Position',[0.025 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,vsf)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.005])
set(gcf,'renderer','painters')
text(0,1.75,'SF','HorizontalAlignment','center')

%B - sp
subplot('Position',[0.025 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,vsp)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.001])
set(gcf,'renderer','painters')
text(0,1.75,'SP','HorizontalAlignment','center')

%C - mp
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,vmp)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.05])
set(gcf,'renderer','painters')
text(0,1.75,'MP','HorizontalAlignment','center')

%D - lp
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,vlp)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 2])
set(gcf,'renderer','painters')
text(0,1.75,'LP','HorizontalAlignment','center')

%E - mf
subplot('Position',[0.5 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,vmf)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 5])
set(gcf,'renderer','painters')
text(0,1.75,'MF','HorizontalAlignment','center')

%F - sd
subplot('Position',[0.5 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,vsd)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.001])
%colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'SD','HorizontalAlignment','center')

%G - md
subplot('Position',[0.5 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,vmd)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.3])
set(gcf,'renderer','painters')
text(0,1.75,'MD','HorizontalAlignment','center')

%H - ld
subplot('Position',[0.5 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,vld)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 2])
set(gcf,'renderer','painters')
text(0,1.75,'LD','HorizontalAlignment','center')

print('-dpng',[ppath 'Map_FEISTY_CORE_interann_var_stages.png'])

% 8plot by type 
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - s
subplot('Position',[0.025 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,vs)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'Sm','HorizontalAlignment','center')

%B - m
subplot('Position',[0.025 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,vm)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 6])
set(gcf,'renderer','painters')
text(0,1.75,'Md','HorizontalAlignment','center')

%C - l
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,vl)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 6])
set(gcf,'renderer','painters')
text(0,1.75,'Lg','HorizontalAlignment','center')

%D - F
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,vf)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 6])
set(gcf,'renderer','painters')
text(0,1.75,'F','HorizontalAlignment','center')

%E - P
subplot('Position',[0.5 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,vp)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 6])
set(gcf,'renderer','painters')
text(0,1.75,'P','HorizontalAlignment','center')

%F - D
subplot('Position',[0.5 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,vd)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 6])
%colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'D','HorizontalAlignment','center')

%G - B
subplot('Position',[0.5 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,vb)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 100])
set(gcf,'renderer','painters')
text(0,1.75,'Bent','HorizontalAlignment','center')

%H - all
subplot('Position',[0.5 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,va)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 10])
set(gcf,'renderer','painters')
text(0,1.75,'All','HorizontalAlignment','center')

print('-dpng',[ppath 'Map_FEISTY_CORE_interann_var_types.png'])

%% save
save([fpath 'FEISTY_CORE_interann_var.mat'],...
    'vsf','vsp','vsd','vmf','vmp','vmd','vlp','vld','vb','va','vs','vm','vl',...
    'vf','vp','vd','vpf','vpd','vlm',...
    'lme_sf_var1','lme_sp_var1','lme_sd_var1',...
    'lme_mf_var1','lme_mp_var1','lme_md_var1',...
    'lme_lp_var1','lme_ld_var1',...
    'lme_b_var1','lme_a_var1',...
    'lme_s_var1','lme_m_var1','lme_l_var1',...
    'lme_f_var1','lme_p_var1','lme_d_var1',...
    'lme_pf_var1','lme_pd_var1','lme_lm_var1');
save([fpath 'FEISTY_CORE_ann_mean_anoms.mat'],...
    'asf','asp','asd','amf','amp','amd','alp','ald','ab','aa','as','am','al',...
    'af','ap','ad','apf','apd','alm');



