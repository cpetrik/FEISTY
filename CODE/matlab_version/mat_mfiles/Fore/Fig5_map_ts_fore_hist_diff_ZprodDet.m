% Visualize difference between
% ESM2M Hindcast of 1951-2000 
% and Forecast of 2051-2100

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ffold=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
%ffold=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

%% Zoop & Det
load([bpath 'cobalt_det_temp_zoop_npp_means.mat'],'mzprod_mean_hist',...
    'mzprod_mean_fore','lzprod_mean_hist','lzprod_mean_fore','det_mean_hist',...
    'det_mean_fore');
load([ffold 'ESM2M_Hist_Fore/ts_Hist_Fore_Zp_D_ZpDet_intA.mat']);

% molN/m2/s --> g/m2/d
mzprod_hist = mzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_hist = lzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_hist = det_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

mzprod_fore = mzprod_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_fore = lzprod_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_fore = det_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

zprod_hist = mzprod_hist + lzprod_hist;
zprod_fore = mzprod_fore + lzprod_fore;

ZpDet_hist = (zprod_hist./det_hist);
ZpDet_fore = (zprod_fore./det_fore);

%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

% Hindcast
load([ffold 'Historic_ESM2M/Means_Historic_' harv '_' cfile '.mat'],...
    'b_mean50');

[hi,hj]=size(geolon_t);
Hb =NaN*ones(hi,hj);
Hb(grid(:,1)) =b_mean50;

clear b_mean50


% Forecast
load([ffold 'Forecast_RCP85_ESM2M/Means_fore_',harv,'_' cfile '.mat'],...
    'b_mean50');

[ni,nj]=size(geolon_t);
Cb =NaN*ones(ni,nj);
Cb(ID) =b_mean50;

clear b_mean50


%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

%%
diffZD = (ZpDet_fore-ZpDet_hist);
pdiffZD = (ZpDet_fore-ZpDet_hist) ./ ZpDet_hist;

pdiffDet = (det_fore-det_hist) ./ det_hist;
pdiffZ = (zprod_fore-zprod_hist) ./ zprod_hist;

%time series of 5-yr means as difference from 1951
dtD = tD - tD(19);
dtZ = tZ - tZ(19);
dtZD = tZD - tZD(19);

pdtD = (tD - tD(19)) / tD(19);
pdtZ = (tZ - tZ(19)) / tZ(19);
pdtZD = (tZD - tZD(19)) / tZD(19);

%% ts colors
% USE RED AND BLUE
% cmap = [...
%     0            0         0.65   %blue
%     0       0            0   %black
%     0.97647         0.19            0]; %red

 cmap = [...
     0.5            0.5         0.5   %blue
           0       0            0]; %black

set(groot,'defaultAxesColorOrder',cmap);

%% Map Zprod:Det w/ ts as subplot 
figure(1)
%subplot('Position',[0 0.05 0.5 0.4])
subplot('Position',[0.4 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffZD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-20 20]);
colorbar('Position',[0.45 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Change in Zoo:Det');
text(-2.9,1.75,'A')

%ts 
subplot('Position',[0.46 0.075 0.385 0.4])
yyaxis left
line(y(19:end),100*pdtD(19:end),'Linewidth',2); hold on;
line(y(19:end),100*pdtZ(19:end),'LineStyle','-.','Linewidth',2); hold on;
ylabel('Percent change in production');

yyaxis right
line(y(19:end),dtZD(19:end),'Linewidth',2); hold on;
xlabel('Year')
ylabel('Change in Zoo:Det');
text(1929,0.42,'B')
print('-dpng',[pp 'Hist_Fore_' harv '_global_diff_ZpDet_2plot_ms.png'])

%% Mollweide Map Zprod:Det w/ ts as subplot 
figure(2)
%subplot('Position',[0 0.05 0.5 0.4])
subplot('Position',[0.4 0.575 0.5 0.4])
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffZD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-20 20]);
colorbar('Position',[0.45 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Change in Zoo:Det');
text(-2.9,1.75,'A')

%ts 
subplot('Position',[0.46 0.075 0.385 0.4])
yyaxis left
line(y(19:end),100*pdtD(19:end),'Linewidth',2); hold on;
line(y(19:end),100*pdtZ(19:end),'LineStyle','-.','Linewidth',2); hold on;
ylabel('Percent change in production');

yyaxis right
line(y(19:end),dtZD(19:end),'Linewidth',2); hold on;
xlabel('Year')
ylabel('Change in Zoo:Det');
text(1929,0.42,'B')
print('-dpng',[pp 'Hist_Fore_' harv '_global_diff_ZpDet_2plot_ms_mollw.png'])

