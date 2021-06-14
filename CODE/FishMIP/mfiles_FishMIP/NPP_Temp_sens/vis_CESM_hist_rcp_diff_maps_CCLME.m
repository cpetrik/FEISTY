% Visualize output of FEISTY forced with CESM
% Historic 1850-2005 initialized with spinup biomass
% CCLME maps of 2090-2100 vs 1860-1870

clear all
close all

%% Map data
cpath = '/Volumes/FEISTY/Fish-MIP/CMIP5/CESM/';
load([cpath 'gridspec_cesm.mat']);
load([cpath 'Data_grid_cesm.mat']);
[ni,nj]=size(LON);
% plotminlat=-90; 
% plotmaxlat=90;
% plotminlon=-280;
% plotmaxlon=80;
plotminlat=25; 
plotmaxlat=55;
plotminlon=-134;
plotmaxlon=-110;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines

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

%% Hist
load([fpath 'Means_Historic_' harv '_' cfile '.mat'],'time','mo',...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean');

% Comparing 2090-2100 to 1990-2000
Hyr = mo(1:12:end);
yrP=find(Hyr>1990 & Hyr<=2000);  

sp_meanP=mean(sp_mean(:,yrP),2);
sf_meanP=mean(sf_mean(:,yrP),2);
sd_meanP=mean(sd_mean(:,yrP),2);
mp_meanP=mean(mp_mean(:,yrP),2);
mf_meanP=mean(mf_mean(:,yrP),2);
md_meanP=mean(md_mean(:,yrP),2);
lp_meanP=mean(lp_mean(:,yrP),2);
ld_meanP=mean(ld_mean(:,yrP),2);
b_meanP = mean(b_mean(:,yrP),2);

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

%% Fore
load([fpath 'Means_Forecast_' harv '_' cfile '.mat'],'time','mo',...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean');

% Comparing 2090-2100 to 1990-2000
Fyr = mo(1:12:end);
yrF=find(Fyr>2090 & Fyr<=2100); 

sp_meanF=mean(sp_mean(:,yrF),2);
sf_meanF=mean(sf_mean(:,yrF),2);
sd_meanF=mean(sd_mean(:,yrF),2);
mp_meanF=mean(mp_mean(:,yrF),2);
mf_meanF=mean(mf_mean(:,yrF),2);
md_meanF=mean(md_mean(:,yrF),2);
lp_meanF=mean(lp_mean(:,yrF),2);
ld_meanF=mean(ld_mean(:,yrF),2);
b_meanF = mean(b_mean(:,yrF),2);

%2090-2100
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

PFracPD = PAllP ./ (PAllP+PAllD);
PFracPF = PAllP ./ (PAllP+PAllF);
PFracLM = PAllL ./ (PAllL+PAllM);
PFracD = PAllD ./ (PAll);
PFracF = PAllF ./ (PAll);
PFracP = PAllP ./ (PAll);
PFracB = Pb ./ (PAll + Pb);

FAll = Fsp+Fsf+Fsd+Fmp+Fmf+Fmd+Flp+Fld;
FAllF = Fsf+Fmf;
FAllP = Fsp+Fmp+Flp;
FAllD = Fsd+Fmd+Fld;
FAllS = Fsp+Fsf+Fsd;
FAllM = Fmp+Fmf+Fmd;
FAllL = Flp+Fld;

FFracPD = FAllP ./ (FAllP+FAllD);
FFracPF = FAllP ./ (FAllP+FAllF);
FFracLM = FAllL ./ (FAllL+FAllM);
FFracD = FAllD ./ (FAll);
FFracF = FAllF ./ (FAll);
FFracP = FAllP ./ (FAll);
FFracB = Fb ./ (FAll + Fb);

pdiffA = (FAll - PAll) ./ PAll;
pdiffF = (FAllF - PAllF) ./ PAllF;
pdiffP = (FAllP - PAllP) ./ PAllP;
pdiffD = (FAllD - PAllD) ./ PAllD;
pdiffM = (FAllM - PAllM) ./ PAllM;
pdiffL = (FAllL - PAllL) ./ PAllL;
pdiffB = (Fb - Pb) ./ Pb;

diffPD = (FFracPD - PFracPD);% ./ PFracPD;
diffPF = (FFracPF - PFracPF);% ./ PFracPF;
diffLM = (FFracLM - PFracLM);% ./ PFracLM;
diffD = (FFracD - PFracD);% ./ PFracD;
diffP = (FFracP - PFracP);% ./ PFracP;
diffF = (FFracF - PFracF);% ./ PFracF;
diffB = (FFracB - PFracB);% ./ PFracF;

%% All 4 on subplots 
figure(1)
% all F
subplot('Position',[0 0.51 0.3 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,100*pdiffF)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 50]);
set(gcf,'renderer','painters')
text(-0.30,0.935,'\bf Forage Fish','HorizontalAlignment','center')

% All P
subplot('Position',[0.3 0.51 0.3 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,100*pdiffP)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 50]);
colorbar('Position',[0.6 0.25 0.03 0.5],'orientation','vertical')
set(gcf,'renderer','painters')
text(-0.30,0.935,'\bf Large Pelagics','HorizontalAlignment','center')
text(-0.075,0.80,'RCP 8.5 2090-2100 vs Historic 1990-2000','HorizontalAlign','left')
text(-0.075,0.75,'% change biomass','HorizontalAlign','left')

% all D
subplot('Position',[0 0 0.3 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,100*pdiffD)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 50]);
set(gcf,'renderer','painters')
text(-0.30,0.935,'\bf Groundfish','HorizontalAlignment','center')

% all B
subplot('Position',[0.3 0 0.3 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,100*pdiffB)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 50]);
set(gcf,'renderer','painters')
text(-0.30,0.935,'\bf Benthic Inverts','HorizontalAlignment','center')
print('-dpng',[ppath 'Hist_RCP_',harv,'_CCLME_All_subplot_pdiffs.png'])

%% Large & medium
figure(2)
subplot('Position',[0 0.5 0.49 0.49])
%P:D
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,diffPD)
cmocean('balance')
caxis([-0.25 0.25]);
colorbar('Position',[0.475 0.55 0.03 0.4],'orientation','vertical')
set(gcf,'renderer','painters')
text(-0.315,0.935,'\bf Fraction Lg Pel vs. Dem','HorizontalAlignment','center')
text(2.6,1.75,'RCP 8.5 2090-2100 vs Historic 1860-1870 % change biomass','HorizontalAlign','center')

%P:F
subplot('Position',[0.5 0.5 0.49 0.49])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,diffPF)
cmocean('balance')
caxis([-0.25 0.25]);
set(gcf,'renderer','painters')
text(-0.315,0.935,'\bf Fraction Lg Pel vs. Forage','HorizontalAlignment','center')

%L:M
subplot('Position',[0.25 0.0 0.49 0.49])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','on','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,diffLM)
cmocean('balance')
caxis([-0.25 0.25]);
set(gcf,'renderer','painters')
text(-0.315,0.935,'\bf Fraction Large vs. Medium','HorizontalAlignment','center')
stamp(cfile)
print('-dpng',[ppath 'Hist_RCP_',harv,'_CCLME_fracs_subplot_pdiffs.png'])

%% All 4 fractions on subplots 
figure(3)
% all F
subplot('Position',[0 0.51 0.3 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,diffF)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.2 0.2]);
set(gcf,'renderer','painters')
text(-0.30,0.935,'\bf Forage Fish','HorizontalAlignment','center')

% All P
subplot('Position',[0.3 0.51 0.3 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,diffP)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.2 0.2]);
colorbar('Position',[0.6 0.25 0.03 0.5],'orientation','vertical')
set(gcf,'renderer','painters')
text(-0.30,0.935,'\bf Large Pelagics','HorizontalAlignment','center')
text(-0.075,0.80,'RCP 8.5 2090-2100 vs Historic 1990-2000','HorizontalAlign','left')
text(-0.075,0.75,'change in relative fraction','HorizontalAlign','left')

% all D
subplot('Position',[0 0 0.3 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,diffD)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.2 0.2]);
set(gcf,'renderer','painters')
text(-0.30,0.935,'\bf Groundfish','HorizontalAlignment','center')

% all B
subplot('Position',[0.3 0 0.3 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,diffB)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.2 0.2]);
set(gcf,'renderer','painters')
text(-0.30,0.935,'\bf Benthos','HorizontalAlignment','center')
print('-dpng',[ppath 'Hist_RCP_',harv,'_CCLME_fracs_types_subplot.png'])
