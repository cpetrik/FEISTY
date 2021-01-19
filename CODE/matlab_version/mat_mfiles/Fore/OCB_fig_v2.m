% Figure for OCB Highlight
% Visualize difference between ESM2M Hindcast (1951-2000) and Forecast (2051-2100)
% fractions of F and P production instead of biomass

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); 
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/ESM2M_Hist_Fore/'];
epath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';

%% hist and fore saved together
load([fpath 'Means_hist_fore_',harv,'_prod_' cfile '.mat'],...
    'cF','cP','cD','cS','cM','cL','hF',...
    'hP','hD','hS','hM','hL','Cb','Hb',...
    'HtP','HtF','HtPel','CtP','CtF','CtPel','Hyr','Cyr');
% load([epath 'Means_hist_fore_',harv,'_prod_' cfile '.mat'],...
%     'cF','cP','cD','cS','cM','cL','hF',...
%     'hP','hD','hS','hM','hL','Cb','Hb',...
%     'HtP','HtF','HtPel','CtP','CtF','CtPel','Hyr','Cyr');


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

hPel = hF+hP;
hAll = hF+hP+hD;

hFracPel = hPel ./ (hAll);
cFracPel = cPel ./ (cAll);
diffPel = (cFracPel-hFracPel);

hFracF = hF ./ (hAll);
cFracF = cF ./ (cAll);
diffF = (cFracF-hFracF);

hFracP = hP ./ (hAll);
cFracP = cP ./ (cAll);
diffP = (cFracP-hFracP);

tF = [HtF CtF];
tP = [HtP CtP];
tPel = [HtPel CtPel];

mtF = movmean(tF,12,2);
mtP = movmean(tP,12,2);
mtPel = movmean(tPel,12,2);

%% time series of 12-mo means as difference from 1951
yr = [Hyr Cyr];
id = find(yr>1950);
yid = id(1);

dtF = mtF - mtF(yid);
dtP = mtP - mtP(yid);
dtPel = mtPel - mtPel(yid);

pdtF = (tF - tF(yid)) / tF(yid);
pdtP = (tP - tP(yid)) / tP(yid);
pdtPel = (tPel - tPel(yid)) / tPel(yid);

%% Zoop:Det
% Visualize difference between ESM2M Hindcast (1951-2000) and Forecast (2051-2100)

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
 cmap = [...
     0.5            0.5         0.5   %grey
           0       0            0]; %black

set(groot,'defaultAxesColorOrder',cmap);

%% Bar plot of trophic amplification
% with std dev from ensemble sims var temp and equal temp
% NO benthic production = benthic biomass * detritus flux * benthic efficiency
% calc benthic production = detritus flux * benthic efficiency


cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

%% Varying temp-dep
% epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
%     'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
% load([epath 'Prod_diff_50yr_ensem6_mid_kt2_bestAIC_multFup_multPneg.mat'])
nfile = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
load([nfile 'Prod_diff_50yr_ensem6_mid_kt3_bestAIC_multFup_multPneg.mat'])

% take mean and error bars
%functional groups and sizes
fmean = mean(fbar);
fstd = std(fbar);
fmax = max(fbar);
fmin = min(fbar);

%size
sbar = pbar(1:2);
sbar(3:4) = fmean(1:2);
smean = fmean(1:2);
sstd = fstd(1:2);
smax = fmax(1:2);
smin = fmin(1:2);

%type
tbar = fmean(3:5);
tstd = fstd(3:5);
tmax = fmax(3:5);
tmin = fmin(3:5);

%type with bent
bbar = fmean([3:5 8]);
bstd = fstd([3:5 8]);

clear fbar 

%% bar graphs
% BW
%figure(5)
figure('Units','inches','Position',[1 3 7 7.5]);
%A  Panel
%subplot('Position',[0.08 0.7 0.4 0.275])
subplot('Position',[0.08 0.7 0.375 0.275])
b=bar(sbar*100,'k'); hold on;
b.FaceColor = [0.65 0.65 0.65];
ylim([-27 0])
er = errorbar(3:4,100*sbar(3:4),100*sstd,'LineWidth',2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
ylabel('Percent change')
set(gca,'XTickLabel','')
text(1,-28,'NPP','HorizontalAlignment','right','Rotation',45)
text(2,-28,{'Meso', 'Zoo'},'HorizontalAlignment','right','Rotation',45)
text(3,-28,{'Med','Fish'},'HorizontalAlignment','right','Rotation',45)
text(4,-28,{'Lrg','Fish'},'HorizontalAlignment','right','Rotation',45)

text(0,-1,'a','FontWeight','bold')

% Type & size subplot
%B  Panel
subplot('Position',[0.555 0.7 0.375 0.275])
b=bar(bbar*100,'k'); hold on;
b.FaceColor = [0.65 0.65 0.65];
ylim([-27 0])
er = errorbar(1:3,100*bbar(1:3),100*bstd(1:3),'LineWidth',1.5);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
set(gca,'XTickLabel','')
text(1,-28,'Forage','HorizontalAlignment','right','Rotation',45)
text(2,-28,{'Large', 'Pelagic'},'HorizontalAlignment','right','Rotation',45)
text(3,-28,'Demersal','HorizontalAlignment','right','Rotation',45)
text(4,-28,{'Benthic', 'invert'},'HorizontalAlignment','right','Rotation',45)
ylabel('Percent change')
text(0,-1,'b','FontWeight','bold')

%% Mollweide Map Zprod:Det w/ ts as subplot 
%C  Panel
subplot('Position',[0.01 0.28 0.42 0.39])
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffZD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-20 20]);
colorbar('Position',[0.02 0.32 0.38 0.02],'orientation','horizontal')
set(gcf,'renderer','painters')
%title('Change in Zoo:Det');
text(0,1.75,'\Delta Zoo prod:Det flux','HorizontalAlignment','center')
text(-2.75,1.75,'c','FontWeight','bold')

%ts 
%D  Panel
subplot('Position',[0.53 0.34 0.39 0.26])
yyaxis left
line(y(19:end),100*pdtD(19:end),'Linewidth',2); hold on;
line(y(19:end),100*pdtZ(19:end),'LineStyle','-.','Linewidth',2); hold on;
ylabel('Percent change in production');

yyaxis right
line(y(19:end),dtZD(19:end),'Linewidth',2); hold on;
xlabel('Year')
ylabel('\Delta Zoo prod:Det flux');
text(1955,0.38,'d','FontWeight','bold')



%E  Panel
% subplot('Position',[0.01 0.04 0.4 0.3])
% 
% %F Panel
% subplot('Position',[0.525 0.06 0.4 0.275])

print('-dpng',[pp 'Hist_Fore_' harv '_OCB_bar_map_ts_v2.png'])


