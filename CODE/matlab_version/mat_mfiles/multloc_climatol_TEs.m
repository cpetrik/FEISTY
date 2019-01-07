% Visualize output of POEM Climatology globally
% 150 years, monthly means saved
% Transfer efficiency ("effective") 

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

% POEM
cfile = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);


%% Zoop and det
gpath='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([gpath 'clim_det_biom_Dmeans_Ytot.mat'])

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
Plb=NaN*ones(ni,nj);

Psf(ID)=sf_mprod;
Psp(ID)=sp_mprod;
Psd(ID)=sd_mprod;
Pmf(ID)=mf_mprod;
Pmp(ID)=mp_mprod;
Pmd(ID)=md_mprod;
Plp(ID)=lp_mprod;
Pld(ID)=ld_mprod;
Plb(ID)=b_mean;

All = Psp+Psf+Psd+Pmp+Pmf+Pmd+Plp+Pld;
AllF = Psf+Pmf;
AllP = Psp+Pmp+Plp;
AllD = Psd+Pmd+Pld;
AllS = Psp+Psf+Psd;
AllM = Pmp+Pmf+Pmd;
AllL = Plp+Pld;

%% Production (= nu * biom) efficiency
%Try using losses vs. standing biomass
TEprodS = AllS./mmz_loss;   %0.0020    0.0264    0.0522    0.0855    0.2533
TEprodM = AllM./AllS;       %-0.0484    1.2105    2.0928    3.8484   21.4620
TEprodL = AllL./AllM;       %-0.0748    0.0116    0.0337    0.1292    0.6602

AllM1=(Pmp+Pmf);
AllM2=(Pmd);
AllL1=(Plp);
AllL2=(Pld);

TEprodMpel = AllM1./(AllS+mlz_loss);    %-0.0223    0.0661    0.1704    0.2792    0.8204
TEprodMdem = AllM2./(Plb);              %0.0001    0.0005    0.0011    0.0020    0.0155
TEprodLP = AllL1./AllM1;                %-0.0021   -0.0000    0.0000    0.0403    0.2857
TEprodLD = AllL2./(AllM2+Plb);          %0.0000    0.0001    0.0003    0.0006    0.0087


%% Effective TEs
TEeffM = AllM./(Plb + mmz_loss + mlz_loss); %-0.0004    0.0034    0.0132    0.0347    0.1338
TEeffL = AllL./(Plb + mmz_loss + mlz_loss); %0.0000    0.0001    0.0003    0.0020    0.0179

%TEeff_L = production_L/NPP
%TEeff_LTL = (production_benthic_invert+mesozoo_prod_to_fish)/NPP
%TEeff_HTL = production_L/(production_benthic_invert+mesozoo_prod_to_fish)

%% save
% TEprods don't seem right b/c go over 1

TEM = real(TEeffM.^(1/3));
TEL = real(TEeffL.^(1/4));

save([fpath 'TEeff_Climatol_All_fish03_' cfile '.mat'],'TEeffM','TEeffL',...
    'Pmf','Pmp','Pmd','Plp','Pld','Plb','mmz_loss','mlz_loss');

%% Figures
% all S
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,TEprodS)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.3]);
hcb = colorbar('h');
ylim(hcb,[0 0.3])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology TEprod S')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_TEprodS.png'])

% all M
figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,TEprodM)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.3]);
hcb = colorbar('h');
ylim(hcb,[0 0.3])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology TEprod M')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_TEprodM.png'])

% all L
figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,TEprodL)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.3]);
hcb = colorbar('h');
ylim(hcb,[0 0.3])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology TEprod L')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_TEprodL.png'])

% Mpel
figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,TEprodMpel)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.3]);
hcb = colorbar('h');
ylim(hcb,[0 0.3])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology TEprod MF & MP')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_TEprodMpel.png'])

% MD
figure(5)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,TEprodMdem)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.3]);
hcb = colorbar('h');
ylim(hcb,[0 0.3])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology TEprod MD')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_TEprodMD.png'])

% LP
figure(6)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,TEprodLP)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.3]);
hcb = colorbar('h');
ylim(hcb,[0 0.3])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology TEprod LP')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_TEprodLP.png'])

% LD
figure(7)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,TEprodLD)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.3]);
hcb = colorbar('h');
ylim(hcb,[0 0.3])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology TEprod LD')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_TEprodLD.png'])

%% Effective
% all M
figure(8)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,TEeffM)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.005 0.03]);
hcb = colorbar('h');
ylim(hcb,[0.005 0.03])                   
set(gcf,'renderer','painters')
title('Climatology TEeff M')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_TEeffM.png'])

% all L
figure(9)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,TEeffL)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.0005 0.003]);
hcb = colorbar('h');
ylim(hcb,[0.0005 0.003])                   
set(gcf,'renderer','painters')
title('Climatology TEeff L')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_TEeffL.png'])

%% all M
figure(10)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,TEM)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.40]);
hcb = colorbar('h');
ylim(hcb,[0 0.40])                   
set(gcf,'renderer','painters')
title('Climatology TEeff M')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_TEeffM_converted.png'])

% all L
figure(11)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,TEL)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.40]);
hcb = colorbar('h');
ylim(hcb,[0 0.40])                   
set(gcf,'renderer','painters')
title('Climatology TEeff L')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_TEeffL_converted.png'])

%% All 4 on subplots
% figure(10)
% % all F
% subplot('Position',[0 0.51 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(AllF))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0 0.3]);
% %     hcb = colorbar('h');
% %     ylim(hcb,[0 0.3])                   %Set color axis if needed
% colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
% set(gcf,'renderer','painters')
% title('log10 mean All F (g m^-^2)')
% 
% % all D
% subplot('Position',[0 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(AllD))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0 0.3]);
% %     hcb = colorbar('h');
% %     ylim(hcb,[0 0.3])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('log10 mean All D (g m^-^2)')
% 
% % All P
% subplot('Position',[0.5 0.51 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(AllP))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0 0.3]);
% %     hcb = colorbar('h');
% %     ylim(hcb,[0 0.3])
% set(gcf,'renderer','painters')
% title('log10 mean All P (g m^-^2)')
% 
% % All
% subplot('Position',[0.5 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(All))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0 0.3]);
% %     hcb = colorbar('h');
% %     ylim(hcb,[0 0.3])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('log10 mean All fishes (g m^-^2)')
% %     stamp([harv '_' cfile])
% %print('-dpng',[ppath 'Climatol_' harv '_global_All_subplot.png'])
% 
