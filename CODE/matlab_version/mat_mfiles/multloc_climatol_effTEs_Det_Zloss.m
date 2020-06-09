% Visualize output of POEM Climatology globally
% 150 years, monthly means saved
% Transfer efficiency ("effective") 
% Use BE*det

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

% POEM
%cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_2B_BE08_noCC_RE00100';
BE = 0.075;
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
%load([fpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);
load([fpath 'Means_Climatol_' harv '_' cfile '.mat']);


%% Zoop and det and npp
gpath='/Volumes/FEISTY/GCM_DATA/ESM26_hist/';
load([gpath 'clim_det_biom_Dmeans_Ytot.mat'])
load([gpath 'clim_npp_Dmeans_Ytot.mat'])

%ESM2.6 in mg C m-2 or mg C m-2 d-1
%from mg C m-2 to g(WW) m-2
% 1e-3 g C in 1 mg C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)

mmz_mean = mz_mean_clim * 1e-3 * 9.0;
mlz_mean = lz_mean_clim * 1e-3 * 9.0;
mmz_loss = mzloss_mean_clim * 1e-3 * 9.0;
mlz_loss = lzloss_mean_clim * 1e-3 * 9.0;

tmz_mean = mz_tot_clim * 1e-3 * 9.0;
tlz_mean = lz_tot_clim * 1e-3 * 9.0;
tmz_loss = mzl_tot_clim * 1e-3 * 9.0;
tlz_loss = lzl_tot_clim * 1e-3 * 9.0;

mdet = det_mean_clim* 1e-3 * 9.0;
tdet = det_tot_clim* 1e-3 * 9.0;

mnpp = npp_mean_clim* 1e-3 * 9.0;
tnpp = npp_tot_clim* 1e-3 * 9.0;

%% plot info
[ni,nj]=size(lon);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

geolat_t=lat;
geolon_t=lon;

Psf=NaN*ones(ni,nj);
Psp=NaN*ones(ni,nj);
Psd=NaN*ones(ni,nj);
Pmf=NaN*ones(ni,nj);
Pmp=NaN*ones(ni,nj);
Pmd=NaN*ones(ni,nj);
Plp=NaN*ones(ni,nj);
Pld=NaN*ones(ni,nj);

Psf(ID)=sf_mprod;
Psp(ID)=sp_mprod;
Psd(ID)=sd_mprod;
Pmf(ID)=mf_mprod;
Pmp(ID)=mp_mprod;
Pmd(ID)=md_mprod;
Plp(ID)=lp_mprod;
Pld(ID)=ld_mprod;

Psf(Psf(:)<0) = 0;
Psp(Psp(:)<0) = 0;
Psd(Psd(:)<0) = 0;
Pmf(Pmf(:)<0) = 0;
Pmp(Pmp(:)<0) = 0;
Pmd(Pmd(:)<0) = 0;
Plp(Plp(:)<0) = 0;
Pld(Pld(:)<0) = 0;

All = Psp+Psf+Psd+Pmp+Pmf+Pmd+Plp+Pld;
AllF = Psf+Pmf;
AllP = Psp+Pmp+Plp;
AllD = Psd+Pmd+Pld;
AllS = Psp+Psf+Psd;
AllM = Pmp+Pmf+Pmd;
AllL = Plp+Pld;

%% Effective TEs
% With BE*det instead of Bent
TEeffM = AllM./(BE*mdet + mmz_loss + mlz_loss); 
%TEeff_L = production_L/NPP
TEeff_ATL = AllL./mnpp;
TEeff_ATL(TEeff_ATL==-Inf) = NaN;
TEeff_ATL(TEeff_ATL==Inf) = NaN;
TEeff_ATL(TEeff_ATL<0) = NaN;
%TEeff_LTL = (production_benthic_invert+mesozoo_prod_to_fish)/NPP
TEeff_LTL = (BE*mdet + mmz_loss + mlz_loss)./mnpp;
TEeff_LTL(TEeff_LTL==-Inf) = NaN;
TEeff_LTL(TEeff_LTL==Inf) = NaN;
TEeff_LTL(TEeff_LTL<0) = NaN;
%TEeff_HTL = production_L/(production_benthic_invert+mesozoo_prod_to_fish)
TEeff_HTL = AllL./(BE*mdet + mmz_loss + mlz_loss); 
TEeff_HTL(TEeff_HTL<0) = NaN;

TEM = real(TEeffM.^(1/2));
TEATL = real(TEeff_ATL.^(1/4));
TELTL = real(TEeff_LTL.^(1/1.3333));
TEHTL = real(TEeff_HTL.^(1/3));

q(:,1) = [0.01 0.05 0.25 0.5 0.75 0.95 0.99]';
q(:,2) = quantile((TEeff_LTL(:)),[0.01 0.05 0.25 0.5 0.75 0.95 0.99])';
q(:,3) = quantile((TEeff_HTL(:)),[0.01 0.05 0.25 0.5 0.75 0.95 0.99]);
q(:,4) = quantile((TEeff_ATL(:)),[0.01 0.05 0.25 0.5 0.75 0.95 0.99]);
q(:,5) = quantile((TELTL(:)),[0.01 0.05 0.25 0.5 0.75 0.95 0.99]);
q(:,6) = quantile((TEHTL(:)),[0.01 0.05 0.25 0.5 0.75 0.95 0.99]);
q(:,7) = quantile((TEATL(:)),[0.01 0.05 0.25 0.5 0.75 0.95 0.99]);


Q = array2table(q,'VariableNames',{'Quantile','TEeff_LTL','TEeff_HTL',...
    'TEeff_ATL','TE_LTL','TE_HTL','TE_ATL'});

%% save
writetable(Q,[fpath 'TEeff_quant_Climatol_All_fish03_' cfile '.csv'],'Delimiter',',');

save([fpath 'TEeffDet_Climatol_All_fish03_' cfile '.mat'],'TEeffM',...
    'Pmf','Pmp','Pmd','Plp','Pld','mmz_loss','mlz_loss','mnpp',...
    'TEeff_ATL','TEeff_LTL','TEeff_HTL');

%% Figures

% All 4 on subplots
%Detritus----------------------
figure(2)
subplot('Position',[0 0.53 0.5 0.5])
%LTL
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(TEeff_LTL))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 -0.5]);
colorbar('Position',[0.025 0.555 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('log10 TEeff LTL')
text(-2.75,1.25,'A')

%HTL
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(TEeff_HTL))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3.5 -1]);
colorbar('Position',[0.525 0.555 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('log10 TEeff HTL')
text(-2.75,1.25,'B')

%ATL
subplot('Position',[0.25 0.025 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(TEeff_ATL))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5.5 -1.5]);
colorbar('Position',[0.275 0.05 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('log10 TEeff ATL')
text(-2.75,1.25,'C')
%stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_global_DeffTEs_subplot.png'])

%% All 3 converted on subplots
%Detritus----------------------
figure(3)
subplot('Position',[0 0.53 0.5 0.5])
%LTL
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(TELTL))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.25]);
colorbar('Position',[0.025 0.555 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('A. TE LTL')

%ATL
subplot('Position',[0.25 0.025 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,TEATL)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.3]);
colorbar('Position',[0.275 0.05 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('C. TE ATL')

subplot('Position',[0.5 0.53 0.5 0.5])
%HTL
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,TEHTL)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.4]);
colorbar('Position',[0.525 0.555 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('B. TE HTL')
%stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_global_DTEs_subplot.png'])


