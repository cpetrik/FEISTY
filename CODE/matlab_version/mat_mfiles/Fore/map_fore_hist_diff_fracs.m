% Visualize difference between
% ESM2M Hindcast of 1951-2000 
% and Forecast of 2051-2100

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

% colors
cmG=cbrewer('seq','Greens',50,'PCHIP');
cmB=cbrewer('seq','Blues',50,'PCHIP');
cmP=cbrewer('seq','Purples',50,'PCHIP');
cmR=cbrewer('seq','Reds',50,'PCHIP');
cmBW=cbrewer('seq','Greys',50,'PCHIP');

%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
%fpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

% save hist and fore together
load([fpath 'Means_hist_fore_',harv,'_cobalt_' cfile '.mat']);

load([fpath 'Hist_Fore_',harv,'_FPDfracs_5yr_ts.mat']);

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

%%
cPel = cF+cP;
cAll = cF+cP+cD;
cFracPD = cP ./ (cP+cD);
cFracPF = cP ./ (cP+cF);
cFracDF = cD ./ (cD+cF);
cFracFD = cF ./ (cF+cD);
cFracFP = cF ./ (cF+cP);
cFracLM = cL ./ (cL+cM);
cFracML = cM ./ (cM+cL);

hPel = hF+hP;
hAll = hF+hP+hD;
hFracPD = hP ./ (hP+hD);
hFracPF = hP ./ (hP+hF);
hFracDF = hD ./ (hD+hF);
hFracFD = hF ./ (hF+hD);
hFracFP = hF ./ (hF+hP);
hFracLM = hL ./ (hL+hM);
hFracML = hM ./ (hM+hL);


pdiffL = (cL-hL) ./ hL;
pdiffM = (cM-hM) ./ hM;
pdiffF = (cF-hF) ./ hF;
pdiffP = (cP-hP) ./ hP;
pdiffD = (cD-hD) ./ hD;
pdiffB = (cB-hB) ./ hB;
pdiffAll = (cAll-hAll) ./ hAll;
pdiffPD = (cFracPD-hFracPD) ./ hFracPD;
pdiffPF = (cFracPF-hFracPF) ./ hFracPF;
pdiffDF = (cFracDF-hFracDF) ./ hFracDF;
pdiffFD = (cFracFD-hFracFD) ./ hFracFD;
pdiffFP = (cFracFP-hFracFP) ./ hFracFP;
pdiffLM = (cFracLM-hFracLM) ./ hFracLM;

diffPD = (cFracPD-hFracPD);
diffPF = (cFracPF-hFracPF);
diffDF = (cFracDF-hFracDF);
diffFD = (cFracFD-hFracFD);
diffFP = (cFracFP-hFracFP);
diffLM = (cFracLM-hFracLM);
diffML = (cFracML-hFracML);

hFracPelD = hPel ./ (hPel+hD);
cFracPelD = cPel ./ (cPel+cD);
diffPelD = (cFracPelD-hFracPelD);

%% time series of 5-yr means as difference from 1951
%'tPD','tPF','tFD','tFP','tPelD','tLM','y','tZD','mtemp')
dtFD = tFD - tFD(19);
dtPD = tPD - tPD(19);
dtPelD = tPelD - tPelD(19);
dtZD = tZD - tZD(19);
dtemp = mtemp - mtemp(19);

pdtFD = (tFD - tFD(19)) / tFD(19);
pdtPD = (tPD - tPD(19)) / tPD(19);
pdtPelD = (tPelD - tPelD(19)) / tPelD(19);

%% ts colors
cm = [...
     0.57255      0.58431      0.56863   %grey
           1       0.8431            0   %yellow
     0.97647         0.19            0   %red
           0            0         0.65   %blue 
         0.4          0.2            0   %brown
         0.1         0.65      0.10196]; %green 

%set(groot,'defaultAxesColorOrder',cm);

cmap = [...
     0       0         0.65   %blue
     0       0            0   %black
     0.1     0.65      0.10196 %green
     0.97647 0.19           0]; %red

set(groot,'defaultAxesColorOrder',cmap);


%% Individual Zprod, Det w/ ts as subplot pdiff Z, D
% USE RED AND BLUE

figure(1)
subplot('Position',[0 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.25 0.55 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Change P/(P+D)');
text(-2.75,1.75,'A')

subplot('Position',[0 0.05 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffPelD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('Change Pel/(Pel+D)');
text(-2.75,2,'B')

% Zp:Det diff
subplot('Position',[0.5 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffFD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
%colorbar('Position',[0.55 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Change in F/(F+D)');
text(-2.75,1.75,'C')

%ts 
subplot('Position',[0.6 0.075 0.375 0.4])
line(y(19:end),dtPD(19:end),'Linewidth',2); hold on;
line(y(19:end),dtFD(19:end),'color',[0.97647 0.19 0],'Linewidth',2); hold on;
line(y(19:end),dtPelD(19:end),'color',[0.1 0.65 0.10196],'Linewidth',2); hold on;
legend('P','F','Pel')
legend('location','northwest')
xlabel('Year')
ylabel('Change in fractions');
text(1918,0.03,'D')
print('-dpng',[pp 'Hist_Fore_' harv '_global_diff_fracs_4plot.png'])

