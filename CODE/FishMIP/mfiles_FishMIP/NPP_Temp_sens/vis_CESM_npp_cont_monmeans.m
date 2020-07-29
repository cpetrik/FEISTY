% Visualize output of FEISTY forced with CESM
% NPP control 1800-2100 initialized with spinup biomass
% Global maps

clear all
close all

%% Fish data
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP_CESM/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'pristine';
tharv = 'F=0';
fpath=['/Volumes/GFDL/NC/FishMIP/CESM1-BEC/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'Means_NPP_cont_' cfile '.mat']);

%% Map data
cpath = '/Volumes/GFDL/Fish-MIP/CESM/';
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
cmBP=cbrewer('seq','BuPu',50);
cmGy=cbrewer('seq','Greys',50);


%% Plots in space
% 1851-1900
Psf=NaN*ones(ni,nj);
Psp=NaN*ones(ni,nj);
Psd=NaN*ones(ni,nj);
Pmf=NaN*ones(ni,nj);
Pmp=NaN*ones(ni,nj);
Pmd=NaN*ones(ni,nj);
Plp=NaN*ones(ni,nj);
Pld=NaN*ones(ni,nj);
Pb=NaN*ones(ni,nj);

Psf(GRD.ID)=sf_mean18;
Psp(GRD.ID)=sp_mean18;
Psd(GRD.ID)=sd_mean18;
Pmf(GRD.ID)=mf_mean18;
Pmp(GRD.ID)=mp_mean18;
Pmd(GRD.ID)=md_mean18;
Plp(GRD.ID)=lp_mean18;
Pld(GRD.ID)=ld_mean18;
Pb(GRD.ID)=b_mean18;

% 1951-2000
Hsf=NaN*ones(ni,nj);
Hsp=NaN*ones(ni,nj);
Hsd=NaN*ones(ni,nj);
Hmf=NaN*ones(ni,nj);
Hmp=NaN*ones(ni,nj);
Hmd=NaN*ones(ni,nj);
Hlp=NaN*ones(ni,nj);
Hld=NaN*ones(ni,nj);
Hb=NaN*ones(ni,nj);

Hsf(GRD.ID)=sf_mean19;
Hsp(GRD.ID)=sp_mean19;
Hsd(GRD.ID)=sd_mean19;
Hmf(GRD.ID)=mf_mean19;
Hmp(GRD.ID)=mp_mean19;
Hmd(GRD.ID)=md_mean19;
Hlp(GRD.ID)=lp_mean19;
Hld(GRD.ID)=ld_mean19;
Hb(GRD.ID)=b_mean19;

% 2051-2100
Fsf=NaN*ones(ni,nj);
Fsp=NaN*ones(ni,nj);
Fsd=NaN*ones(ni,nj);
Fmf=NaN*ones(ni,nj);
Fmp=NaN*ones(ni,nj);
Fmd=NaN*ones(ni,nj);
Flp=NaN*ones(ni,nj);
Fld=NaN*ones(ni,nj);
Fb=NaN*ones(ni,nj);

Fsf(GRD.ID)=sf_mean20;
Fsp(GRD.ID)=sp_mean20;
Fsd(GRD.ID)=sd_mean20;
Fmf(GRD.ID)=mf_mean20;
Fmp(GRD.ID)=mp_mean20;
Fmd(GRD.ID)=md_mean20;
Flp(GRD.ID)=lp_mean20;
Fld(GRD.ID)=ld_mean20;
Fb(GRD.ID)=b_mean20;

ocean=NaN*ones(ni,nj);
ocean(GRD.ID)=ones(size(sf_mean18));

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

HAll = Hsp+Hsf+Hsd+Hmp+Hmf+Hmd+Hlp+Hld;
HAllF = Hsf+Hmf;
HAllP = Hsp+Hmp+Hlp;
HAllD = Hsd+Hmd+Hld;
HAllS = Hsp+Hsf+Hsd;
HAllM = Hmp+Hmf+Hmd;
HAllL = Hlp+Hld;
HFracPD = HAllP ./ (HAllP+HAllD);
HFracPF = HAllP ./ (HAllP+HAllF);
HFracLM = HAllL ./ (HAllL+HAllM);

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

%% All 4 on subplots 1851-1900
figure(1)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(PAllF))
%colormap(cmBP)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar('Position',[0.25 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
%title('log_1_0 mean All F (g m^-^2)')
text(2.6,1.75,'NPP control 1851-1900 log_1_0 mean biomass (g m^-^2)','HorizontalAlign','center')
text(-2.75,1.25,'F')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(PAllP))
%colormap(cmBP)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
%title('log_1_0 mean All P (g m^-^2)')
text(-2.75,1.25,'P')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(PAllD))
%colormap(cmBP)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
%title('log_1_0 mean All D (g m^-^2)')
text(-2.75,1.25,'D')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(PAll))
%colormap(cmBP)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
%title('log_1_0 mean All fishes (g m^-^2)')
text(-2.75,1.25,'All')
stamp(cfile)
print('-dpng',[ppath 'NPP_cont_',harv,'_global_All_subplot_1800s.png'])

%% All 4 on subplots 1951-2000
figure(2)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(HAllF))
%colormap(cmBP)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar('Position',[0.25 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
%title('log_1_0 mean All F (g m^-^2)')
text(2.6,1.75,'NPP control 1951-2000 log_1_0 mean biomass (g m^-^2)','HorizontalAlign','center')
text(-2.75,1.25,'F')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(HAllP))
%colormap(cmBP)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
%title('log_1_0 mean All P (g m^-^2)')
text(-2.75,1.25,'P')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(HAllD))
%colormap(cmBP)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
%title('log_1_0 mean All D (g m^-^2)')
text(-2.75,1.25,'D')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(HAll))
%colormap(cmBP)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
%title('log_1_0 mean All fishes (g m^-^2)')
text(-2.75,1.25,'All')
stamp(cfile)
print('-dpng',[ppath 'NPP_cont_',harv,'_global_All_subplot_1900s.png'])

%% All 4 on subplots 2051-2100
figure(3)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(FAllF))
%colormap(cmBP)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar('Position',[0.25 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
%title('log_1_0 mean All F (g m^-^2)')
text(2.6,1.75,'NPP control 2051-2100 log_1_0 mean biomass (g m^-^2)','HorizontalAlign','center')
text(-2.75,1.25,'F')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(FAllP))
%colormap(cmBP)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
%title('log_1_0 mean All P (g m^-^2)')
text(-2.75,1.25,'P')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(FAllD))
%colormap(cmBP)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
%title('log_1_0 mean All D (g m^-^2)')
text(-2.75,1.25,'D')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(FAll))
%colormap(cmBP)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
%title('log_1_0 mean All fishes (g m^-^2)')
text(-2.75,1.25,'All')
stamp(cfile)
print('-dpng',[ppath 'NPP_cont_',harv,'_global_All_subplot_2000s.png'])

%% Ratios on subplots red-white-blue
% 3 figure subplot P:D, P:F, M:L
% figure(19)
% subplot('Position',[0 0.5 0.49 0.49])
% %P:D
% axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','on','FLineWidth',1)
% surfm(LAT,LON,FracPD)
% cmocean('balance')
% caxis([0 1]);
% set(gcf,'renderer','painters')
% text(0,0.835,'\bf Fraction Large Pelagics vs. Demersals','HorizontalAlignment','center')
% 
% %P:F
% subplot('Position',[0.5 0.5 0.49 0.49])
% axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','on','FLineWidth',1)
% surfm(LAT,LON,FracPF)
% cmocean('balance')
% caxis([0 1]);
% colorbar('Position',[0.475 0.55 0.035 0.4],'orientation','vertical')
% set(gcf,'renderer','painters')
% text(0,0.835,'\bf Fraction Large Pelagics vs. Forage Fishes','HorizontalAlignment','center')
% 
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
% print('-dpng',[ppath 'NPP control_',harv,'_CC_ratios_subplot.png'])
% 
