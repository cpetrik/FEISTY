% Visualize output of FEISTY forced with CESM
% Forecast 2006-2100 initialized with spinup biomass
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
load([fpath 'Means_Forecast_' harv '_' cfile '.mat'],'sf_meanP','sp_meanP','sd_meanP',...
    'mf_meanP','mp_meanP','md_meanP',...
    'lp_meanP','ld_meanP','b_meanP',...
    'sf_meanF','sp_meanF','sd_meanF',...
    'mf_meanF','mp_meanF','md_meanF',...
    'lp_meanF','ld_meanF','b_meanF');

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

%% Plots in space
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
text(2.6,1.75,'Forecast 2051-2100 log_1_0 mean biomass (g m^-^2)','HorizontalAlign','center')
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
print('-dpng',[ppath 'Forecast_',harv,'_global_All_subplot_2000s.png'])

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
% print('-dpng',[ppath 'Forecast_',harv,'_CC_ratios_subplot.png'])
% 
