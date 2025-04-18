% Visualize output of FEISTY forced with CESM
% Temperature control 1800-2100 initialized with spinup biomass
% Global maps of 2090-2100 vs 1860-1870

clear all
close all

%% Fish data
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP_CESM/';
cfile = 'NoDc_enc70-b200_m4-b175-k086_c20-b250_noD_J100_A050_Sm025_nmort1_BE00_noCC_RE00100';
harv = 'pristine';
tharv = 'F=0';
fpath=['/Volumes/FEISTY/NC/FishMIP/CESM1-BEC/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'Means_Temp_cont_' cfile '.mat'],'sf_meanP','sp_meanP',...
    'mf_meanP','mp_meanP',...
    'lp_meanP',...
    'sf_meanF','sp_meanF',...
    'mf_meanF','mp_meanF',...
    'lp_meanF');

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

%% Plots in space
% 1860-1870
Psf=NaN*ones(ni,nj);
Psp=NaN*ones(ni,nj);
Pmf=NaN*ones(ni,nj);
Pmp=NaN*ones(ni,nj);
Plp=NaN*ones(ni,nj);

Psf(GRD.ID)=sf_meanP;
Psp(GRD.ID)=sp_meanP;
Pmf(GRD.ID)=mf_meanP;
Pmp(GRD.ID)=mp_meanP;
Plp(GRD.ID)=lp_meanP;

% 2090-2100
Fsf=NaN*ones(ni,nj);
Fsp=NaN*ones(ni,nj);
Fmf=NaN*ones(ni,nj);
Fmp=NaN*ones(ni,nj);
Flp=NaN*ones(ni,nj);

Fsf(GRD.ID)=sf_meanF;
Fsp(GRD.ID)=sp_meanF;
Fmf(GRD.ID)=mf_meanF;
Fmp(GRD.ID)=mp_meanF;
Flp(GRD.ID)=lp_meanF;

ocean=NaN*ones(ni,nj);
ocean(GRD.ID)=ones(size(sf_meanP));

%% Diff maps of all fish
PAll = Psp+Psf+Pmp+Pmf+Plp;
PAllF = Psf+Pmf;
PAllP = Psp+Pmp+Plp;
PAllS = Psp+Psf;
PAllM = Pmp+Pmf;
PAllL = Plp;

FAll = Fsp+Fsf+Fmp+Fmf+Flp;
FAllF = Fsf+Fmf;
FAllP = Fsp+Fmp+Flp;
FAllS = Fsp+Fsf;
FAllM = Fmp+Fmf;
FAllL = Flp;

pdiffA = (FAll - PAll) ./ PAll;
pdiffF = (FAllF - PAllF) ./ PAllF;
pdiffP = (FAllP - PAllP) ./ PAllP;
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
text(2.6,1.75,'Temp control 2090-2100 vs 1860-1870 % change biomass','HorizontalAlign','center')
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

% All
subplot('Position',[0.25 0 0.5 0.5])
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
print('-dpng',[ppath 'Temp_cont_',harv,'_global_All_subplot_pdiffs.png'])

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
text(2.6,1.75,'Temp control 2090-2100 vs 1860-1870 % change biomass','HorizontalAlign','center')

%L
subplot('Position',[0.5 0.5 0.49 0.49])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,100*pdiffL)
cmocean('balance')
caxis([-100 100]);
set(gcf,'renderer','painters')
text(-2.75,1.25,'Lrg')

print('-dpng',[ppath 'Temp_cont_',harv,'_global_sizes_subplot_pdiffs.png'])
 