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

cDark=cbrewer('qual','Dark2',3);

%% Find regions of Tyle I, II, and III
%I: an incr in zooplankton prod and decr in detritus flux
%II: an incr in zooplankton prod greater than the incr in detritus flux
%III: a decr in zooplankton production less than the decr in detritus flux

typeZD = nan(size(diffZD));
Zincr = nan(size(diffZD));
Dincr = nan(size(diffZD));
Zdecr = nan(size(diffZD));
Ddecr = nan(size(diffZD));

zincr = find(zprod_fore > zprod_hist);
zdecr = find(zprod_fore < zprod_hist);
dincr = find(det_fore > det_hist);
ddecr = find(det_fore < det_hist);

t1 = intersect(zincr,ddecr);

Zincr(zincr) = pdiffZ(zincr);
Dincr(dincr) = pdiffDet(dincr);
t2 = find(Zincr > Dincr);

Zdecr(zdecr) = pdiffZ(zdecr);
Ddecr(ddecr) = pdiffDet(ddecr);
t3 = find(Zdecr > Ddecr);

typeZD(t1) = ones;
typeZD(t2) = 2*ones;
typeZD(t3) = 3*ones;

%individual ones for testing
type1 = nan(size(diffZD));
type2 = nan(size(diffZD));
type3 = nan(size(diffZD));

type1(t1) = ones;
type2(t2) = 2*ones;
type3(t3) = 3*ones;

%% Mollweide Map of 3 types subplot 
figure(1)
subplot(2,2,1)
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,type1)
%cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([-20 20]);
%colorbar('Position',[0.45 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Type I');

subplot(2,2,2)
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,type2)
%cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([-20 20]);
%colorbar('Position',[0.45 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Type II');

subplot(2,2,3)
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,type3)
%cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([-20 20]);
%colorbar('Position',[0.45 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Type III');

subplot(2,2,4)
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,typeZD)
colormap(cDark)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([-20 20]);
colorbar
%colorbar('Position',[0.45 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Types');
print('-dpng',[pp 'Hist_Fore_' harv '_global_diff_ZpDet_types.png'])

%% Mollweide Map Zprod:Det w/ ts as subplot 
figure(2)
subplot('Position',[0.5 0.575 0.5 0.4])
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,typeZD)
colormap(cDark)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([-20 20]);
colorbar('Position',[0.55 0.57 0.4 0.03],'orientation','horizontal',...
    'Ticks',[1.33,2,2.67],'TickLabels',1:3)
set(gcf,'renderer','painters')
title('Types of increase in Zoo:Det');
text(-2.9,1.75,'B')

%subplot('Position',[0 0.05 0.5 0.4])
subplot('Position',[0.0 0.575 0.5 0.4])
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffZD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-20 20]);
colorbar('Position',[0.05 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Change in Zoo:Det');
text(-2.9,1.75,'A')

%ts 
subplot('Position',[0.1 0.075 0.8 0.4])
yyaxis left
line(y(19:end),100*pdtD(19:end),'Linewidth',2); hold on;
line(y(19:end),100*pdtZ(19:end),'LineStyle','-.','Linewidth',2); hold on;
ylabel('Percent change in production');

yyaxis right
line(y(19:end),dtZD(19:end),'Linewidth',2); hold on;
xlabel('Year')
ylabel('Change in Zoo:Det');
text(1934,0.42,'C')
print('-dpng',[pp 'Hist_Fore_' harv '_global_diff_ZpDet_2plot_mollw_types.png'])

%% Map Zprod:Det w/ ts as subplot - percent change map
figure(2)
%subplot('Position',[0 0.05 0.5 0.4])
subplot('Position',[0.4 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,100*pdiffZD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar('Position',[0.45 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Percent change in Zoo:Det');
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
print('-dpng',[pp 'Hist_Fore_' harv '_global_pdiff_ZpDet_2plot_ms.png'])

%% Comp diff & pdiff
figure(3)
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

subplot('Position',[0.4 0.075 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,100*pdiffZD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
%colorbar('Position',[0.45 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Percent change in Zoo:Det');
print('-dpng',[pp 'Hist_Fore_' harv '_global_diff_pdiff_ZpDet_map_comp.png'])


%% Try plot only in places with det above a threshold
did = find(det_hist(:)>1e-2);
ldet = nan*ones(size(det_hist));
ldet(did) = ones;

figure(4)
%subplot('Position',[0 0.05 0.5 0.4])
subplot('Position',[0.4 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffZD.*ldet)
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
print('-dpng',[pp 'Hist_Fore_' harv '_global_diff_ZpDet_2plot_1e-2.png'])


%% Try plot places with Z:D below threshold
zid = find(ZpDet_hist(:)<21.5654);
ldet = nan*ones(size(det_hist));
ldet(zid) = ones;

figure(5)
%subplot('Position',[0 0.05 0.5 0.4])
subplot('Position',[0.4 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffZD.*ldet)
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
print('-dpng',[pp 'Hist_Fore_' harv '_global_diff_ZpDet_2plot_21.png'])

%% Where are big changes
hid=find(diffZD(:)>15);
dep = nan*ones(size(diffZD));
dep(ID) = grid(:,4);

% figure
% subplot(2,2,1)
% hist(det_hist(hid))
% subplot(2,2,2)
% hist(det_fore(hid))
% subplot(2,2,3)
% hist(dep(hid))
% 
% quantile(det_hist(hid),[0.05 0.25 0.5 0.75 0.95])
% quantile(det_fore(hid),[0.05 0.25 0.5 0.75 0.95])
% mean(det_hist(hid))
% mean(det_fore(hid))
% mean(dep(hid))
% 
% max(det_hist(hid))
% max(det_fore(hid))
% 
% %%
% quantile(ZpDet_hist(hid),[0.05 0.25 0.5 0.75 0.95])
% quantile(ZpDet_fore(hid),[0.05 0.25 0.5 0.75 0.95])
% mean(ZpDet_hist(hid))
% mean(ZpDet_fore(hid))
% max(ZpDet_hist(hid))
% max(ZpDet_fore(hid))
% 
% %%
% quantile(ZpDet_hist(:),[0.05 0.25 0.5 0.75 0.95])
% quantile(ZpDet_fore(:),[0.05 0.25 0.5 0.75 0.95])
% nanmean(ZpDet_hist(:))
% nanmean(ZpDet_fore(:))
% max(ZpDet_hist(:))
% max(ZpDet_fore(:))

%
['In areas where diffZD > 15, the mean depth is '...
num2str(mean(dep(hid))) ...
', the mean historic ZD was ' ...
num2str(mean(ZpDet_hist(hid))) ...
' and det flux was '...
num2str(mean(det_hist(hid))) ]

sid=find(dep(:)<=200);
['In shelf areas <= 200 m, the mean change in ZD is '...
num2str(mean(diffZD(sid))) ...
' and the projected ZD is ' ...
num2str(mean(ZpDet_fore(sid)))]

lid=find(dep(:)<=1000);
['In areas <= 1000 m, the mean change in ZD is '...
num2str(mean(diffZD(lid))) ...
' and the projected ZD is ' ...
num2str(mean(ZpDet_fore(lid)))]

