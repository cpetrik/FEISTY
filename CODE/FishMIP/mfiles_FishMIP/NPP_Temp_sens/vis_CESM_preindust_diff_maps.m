% Visualize output of FEISTY forced with CESM
% Preindustrial 1800-2100 initialized with spinup biomass
% Global maps of 2090-2100 vs 1860-1870

clear all
close all

%% Fish data
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP_CESM/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'pristine';
tharv = 'F=0';
fpath=['/Volumes/FEISTY/NC/FishMIP/CESM1-BEC/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'Means_Preindust_' cfile '.mat'],'sf_meanP','sp_meanP','sd_meanP',...
    'mf_meanP','mp_meanP','md_meanP',...
    'lp_meanP','ld_meanP','b_meanP',...
    'sf_meanF','sp_meanF','sd_meanF',...
    'mf_meanF','mp_meanF','md_meanF',...
    'lp_meanF','ld_meanF','b_meanF');

%% Map data
cpath = '/Volumes/FEISTY/Fish-MIP/CESM/';
load([cpath 'gridspec_cesm.mat']);
load([cpath 'Data_grid_cesm.mat']);
[ni,nj]=size(LON);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

% colors
cmBP=cbrewer('seq','BuPu',50,'PCHIP');
cmGy=cbrewer('seq','Greys',50,'PCHIP');


%% Plots in space
% 1860-1870
Psf=NaN*ones(ni,nj);
Psp=NaN*ones(ni,nj);
Psd=NaN*ones(ni,nj);
Pmf=NaN*ones(ni,nj);
Pmp=NaN*ones(ni,nj);
Pmd=NaN*ones(ni,nj);
Plp=NaN*ones(ni,nj);
Pld=NaN*ones(ni,nj);
Pb=NaN*ones(ni,nj);

Psf(GRD.ID)=sf_meanP;
Psp(GRD.ID)=sp_meanP;
Psd(GRD.ID)=sd_meanP;
Pmf(GRD.ID)=mf_meanP;
Pmp(GRD.ID)=mp_meanP;
Pmd(GRD.ID)=md_meanP;
Plp(GRD.ID)=lp_meanP;
Pld(GRD.ID)=ld_meanP;
Pb(GRD.ID)=b_meanP;

% 2090-2100
Fsf=NaN*ones(ni,nj);
Fsp=NaN*ones(ni,nj);
Fsd=NaN*ones(ni,nj);
Fmf=NaN*ones(ni,nj);
Fmp=NaN*ones(ni,nj);
Fmd=NaN*ones(ni,nj);
Flp=NaN*ones(ni,nj);
Fld=NaN*ones(ni,nj);
Fb=NaN*ones(ni,nj);

Fsf(GRD.ID)=sf_meanF;
Fsp(GRD.ID)=sp_meanF;
Fsd(GRD.ID)=sd_meanF;
Fmf(GRD.ID)=mf_meanF;
Fmp(GRD.ID)=mp_meanF;
Fmd(GRD.ID)=md_meanF;
Flp(GRD.ID)=lp_meanF;
Fld(GRD.ID)=ld_meanF;
Fb(GRD.ID)=b_meanF;

ocean=NaN*ones(ni,nj);
ocean(GRD.ID)=ones(size(sf_meanP));

%% Diff maps of all fish
PAll = Psp+Psf+Psd+Pmp+Pmf+Pmd+Plp+Pld;
PAllF = Psf+Pmf;
PAllP = Psp+Pmp+Plp;
PAllD = Psd+Pmd+Pld;
PAllS = Psp+Psf+Psd;
PAllM = Pmp+Pmf+Pmd;
PAllL = Plp+Pld;

FAll = Fsp+Fsf+Fsd+Fmp+Fmf+Fmd+Flp+Fld;
FAllF = Fsf+Fmf;
FAllP = Fsp+Fmp+Flp;
FAllD = Fsd+Fmd+Fld;
FAllS = Fsp+Fsf+Fsd;
FAllM = Fmp+Fmf+Fmd;
FAllL = Flp+Fld;

pdiffA = (FAll - PAll) ./ PAll;
pdiffF = (FAllF - PAllF) ./ PAllF;
pdiffP = (FAllP - PAllP) ./ PAllP;
pdiffD = (FAllD - PAllD) ./ PAllD;
pdiffM = (FAllM - PAllM) ./ PAllM;
pdiffL = (FAllL - PAllL) ./ PAllL;

%% All 4 on subplots 
figure(1)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,100*pdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar('Position',[0.25 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
text(2.6,1.75,'Preindustrial 2090-2100 vs 1860-1870 % change biomass','HorizontalAlign','center')
text(-2.75,1.25,'F')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,100*pdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(-2.75,1.25,'P')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,100*pdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(-2.75,1.25,'D')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,100*pdiffA)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(-2.75,1.25,'All')
%stamp(cfile)
print('-dpng',[ppath 'Preindustrial_',harv,'_global_All_subplot_pdiffs.png'])

%% Large & medium
figure(2)
subplot('Position',[0 0.5 0.49 0.49])
%M
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,100*pdiffM)
cmocean('balance')
caxis([-100 100]);
colorbar('Position',[0.25 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
text(-2.75,1.25,'Med')
text(2.6,1.75,'Preindustrial 2090-2100 vs 1860-1870 % change biomass','HorizontalAlign','center')

%L
subplot('Position',[0.5 0.5 0.49 0.49])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,100*pdiffL)
cmocean('balance')
caxis([-100 100]);
set(gcf,'renderer','painters')
text(-2.75,1.25,'Lrg')

% %L:M
% subplot('Position',[0.25 0.0 0.49 0.49])
% axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','on','FLineWidth',1)
% surfm(LAT,LON,FracLM)
% cmocean('balance')
% caxis([0 1]);
% set(gcf,'renderer','painters')
% text(0,0.835,'\bf Fraction Large vs. Medium','HorizontalAlignment','center')
% stamp(cfile)
print('-dpng',[ppath 'Preindustrial_',harv,'_global_sizes_subplot_pdiffs.png'])
 
