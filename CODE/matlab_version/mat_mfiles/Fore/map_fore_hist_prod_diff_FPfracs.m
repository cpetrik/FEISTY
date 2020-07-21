% Visualize difference between
% ESM2M Hindcast of 1951-2000 
% and Forecast of 2051-2100
% fractions of F and P
% production instead of biomass

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
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

% dtF = tF - tF(yid);
% dtP = tP - tP(yid);
% dtPel = tPel - tPel(yid);

dtF = mtF - mtF(yid);
dtP = mtP - mtP(yid);
dtPel = mtPel - mtPel(yid);

pdtF = (tF - tF(yid)) / tF(yid);
pdtP = (tP - tP(yid)) / tP(yid);
pdtPel = (tPel - tPel(yid)) / tPel(yid);


%%
figure(1)
%P
subplot('Position',[0 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.25 0.55 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Change P/All');
text(-2.75,1.75,'A')
%Pel
subplot('Position',[0 0.05 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffPel)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('Change Pel/All');
text(-2.75,2,'B')
%F
subplot('Position',[0.5 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
%colorbar('Position',[0.55 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Change in F/All');
text(-2.75,1.75,'C')

%ts 
subplot('Position',[0.6 0.075 0.375 0.4])
line(yr(yid:end),dtP(yid:end),'color',[0 0.1 0.9],'Linewidth',2); hold on;
line(yr(yid:end),dtF(yid:end),'color',[0.97647 0.19 0],'Linewidth',2); hold on;
line(yr(yid:end),dtPel(yid:end),'color',[0 0 0],'Linewidth',2); hold on;
set(gca,'Box','on')
ylim([-0.05 0.05])
legend('P','F','Pel')
legend('location','northwest')
xlabel('Year')
ylabel('Change in fractions');
text(1918,0.048,'D')
print('-dpng',[pp 'Hist_Fore_' harv '_global_prod_diff_fracs_4plot_v1.png'])

%% zero line
x0=1951:2100;
y0=zeros(size(x0));
figure(2)
% P
subplot('Position',[0 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar('Position',[0.25 0.55 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Change P/All');
text(-2.75,1.75,'A')
% Pel
subplot('Position',[0 0.05 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffPel)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Change Pel/All');
text(-2.75,2,'B')
% F
subplot('Position',[0.5 0.575 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
%colorbar('Position',[0.55 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Change in F/All');
text(-2.75,1.75,'C')

%ts 
subplot('Position',[0.6 0.075 0.375 0.4])
line(yr(yid:end),dtP(yid:end),'color',[0 0.1 0.9],'Linewidth',2); hold on;
line(yr(yid:end),dtF(yid:end),'color',[0.97647 0.19 0],'Linewidth',2); hold on;
line(yr(yid:end),dtPel(yid:end),'color',[0 0 0],'Linewidth',2); hold on;
legend('P','F','Pel')
legend('location','northwest')
legend('AutoUpdate','off')
line(x0,y0,'LineStyle',':','color',[0 0 0]);
set(gca,'Box','on')
ylim([-0.05 0.05])
xlabel('Year')
ylabel('Change in fractions');
text(1918,0.048,'D')
print('-dpng',[pp 'Hist_Fore_' harv '_global_prod_diff_fracs_4plot_v2.png'])

%% zero line, Mollweid projection
x0=1951:2100;
y0=zeros(size(x0));
figure(3)
% P
subplot('Position',[0 0.575 0.5 0.4])
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar('Position',[0.25 0.55 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Change P/All');
text(-2.75,1.75,'A')
% Pel
subplot('Position',[0 0.05 0.5 0.4])
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffPel)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Change Pel/All');
text(-2.75,2,'B')
% F
subplot('Position',[0.5 0.575 0.5 0.4])
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
%colorbar('Position',[0.55 0.57 0.4 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Change in F/All');
text(-2.75,1.75,'C')

%ts 
subplot('Position',[0.6 0.075 0.375 0.4])
line(yr(yid:end),dtP(yid:end),'color',[0 0.1 0.9],'Linewidth',2); hold on;
line(yr(yid:end),dtF(yid:end),'color',[0.97647 0.19 0],'Linewidth',2); hold on;
line(yr(yid:end),dtPel(yid:end),'color',[0 0 0],'Linewidth',2); hold on;
legend('P','F','Pel')
legend('location','northwest')
legend('AutoUpdate','off')
line(x0,y0,'LineStyle',':','color',[0 0 0]);
set(gca,'Box','on')
ylim([-0.05 0.05])
xlabel('Year')
ylabel('Change in fractions');
text(1918,0.048,'D')
print('-dpng',[pp 'Hist_Fore_' harv '_global_prod_diff_fracs_4plot_v3.png'])


