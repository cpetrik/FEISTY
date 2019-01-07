% Visualize difference between
% ESM2M Hindcast of 1951-2000 and Forecasts of 2051-2100
% Forecasts start from fished state, then test w/ and w/o fishing

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/GoMex/'];

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

%Hindcast fished
load([fpath 'Map_vals_Historic_',harv,'_' cfile '.mat'])
HBtime = b_tmean;
HSFtime = sf_tmean;
HMFtime = mf_tmean;
HSPtime = sp_tmean;
HMPtime = mp_tmean;
HLPtime = lp_tmean;
HSDtime = sd_tmean;
HMDtime = md_tmean;
HLDtime = ld_tmean;
HFtime = F;
HPtime = P;
HDtime = D;
HB = Zb;
HF = AllF;
HP = AllP;
HD = AllD;
HS = AllS;
HM = AllM;
HL = AllL;
HPD = FracPD;
HPF = FracPF;
HLM = FracLM;

%Forecast fished
load([fpath 'Map_vals_Forecast_',harv,'_' cfile '.mat'])
FfBtime = b_tmean;
FfSFtime = sf_tmean;
FfMFtime = mf_tmean;
FfSPtime = sp_tmean;
FfMPtime = mp_tmean;
FfLPtime = lp_tmean;
FfSDtime = sd_tmean;
FfMDtime = md_tmean;
FfLDtime = ld_tmean;
FfFtime = F;
FfPtime = P;
FfDtime = D;
FfB = Zb;
FfF = AllF;
FfP = AllP;
FfD = AllD;
FfS = AllS;
FfM = AllM;
FfL = AllL;
FfPD = FracPD;
FfPF = FracPF;
FfLM = FracLM;

%Forecast pristine
load([fpath 'Map_vals_Forecast_pristine_' cfile '.mat'])
FpBtime = b_tmean;
FpSFtime = sf_tmean;
FpMFtime = mf_tmean;
FpSPtime = sp_tmean;
FpMPtime = mp_tmean;
FpLPtime = lp_tmean;
FpSDtime = sd_tmean;
FpMDtime = md_tmean;
FpLDtime = ld_tmean;
FpFtime = F;
FpPtime = P;
FpDtime = D;
FpB = Zb;
FpF = AllF;
FpP = AllP;
FpD = AllD;
FpS = AllS;
FpM = AllM;
FpL = AllL;
FpPD = FracPD;
FpPF = FracPF;
FpLM = FracLM;

HAll = HF+HP+HD; 
FfAll = FfF+FfP+FfD; 
FpAll = FpF+FpP+FpD;

HAllT = HFtime+HPtime+HDtime; 
FfAllT = FfFtime+FfPtime+FfDtime; 
FpAllT = FpFtime+FpPtime+FpDtime;

%% Year-avg time means
st=1:12:(length(HBtime));
en=12:12:(length(HBtime));

for n=1:length(st)
    HBt(:,n)=nanmean(HBtime(st(n):en(n)));
    HSFt(:,n)=nanmean(HSFtime(st(n):en(n)));
    HSPt(:,n)=nanmean(HSPtime(st(n):en(n)));
    HSDt(:,n)=nanmean(HSDtime(st(n):en(n)));
    HMFt(:,n)=nanmean(HMFtime(st(n):en(n)));
    HMPt(:,n)=nanmean(HMPtime(st(n):en(n)));
    HMDt(:,n)=nanmean(HMDtime(st(n):en(n)));
    HLPt(:,n)=nanmean(HLPtime(st(n):en(n)));
    HLDt(:,n)=nanmean(HLDtime(st(n):en(n)));
end

st=1:12:(length(FpBtime));
en=12:12:(length(FpBtime));
for n=1:length(st)
    FpBt(:,n)=nanmean(FpBtime(st(n):en(n)));
    FpSFt(:,n)=nanmean(FpSFtime(st(n):en(n)));
    FpSPt(:,n)=nanmean(FpSPtime(st(n):en(n)));
    FpSDt(:,n)=nanmean(FpSDtime(st(n):en(n)));
    FpMFt(:,n)=nanmean(FpMFtime(st(n):en(n)));
    FpMPt(:,n)=nanmean(FpMPtime(st(n):en(n)));
    FpMDt(:,n)=nanmean(FpMDtime(st(n):en(n)));
    FpLPt(:,n)=nanmean(FpLPtime(st(n):en(n)));
    FpLDt(:,n)=nanmean(FpLDtime(st(n):en(n)));
    
    FfBt(:,n)=nanmean(FfBtime(st(n):en(n)));
    FfSFt(:,n)=nanmean(FfSFtime(st(n):en(n)));
    FfSPt(:,n)=nanmean(FfSPtime(st(n):en(n)));
    FfSDt(:,n)=nanmean(FfSDtime(st(n):en(n)));
    FfMFt(:,n)=nanmean(FfMFtime(st(n):en(n)));
    FfMPt(:,n)=nanmean(FfMPtime(st(n):en(n)));
    FfMDt(:,n)=nanmean(FfMDtime(st(n):en(n)));
    FfLPt(:,n)=nanmean(FfLPtime(st(n):en(n)));
    FfLDt(:,n)=nanmean(FfLDtime(st(n):en(n)));
end

%% Plots over time
htime = 1861:2005;
ftime = 2006:2100;

HFt = HSFt+HMFt;
HPt = HSPt+HMPt+HLPt;
HDt = HSDt+HMDt+HLDt;

FpFt = FpSFt+FpMFt;
FpPt = FpSPt+FpMPt+FpLPt;
FpDt = FpSDt+FpMDt+FpLDt;

FfFt = FfSFt+FfMFt;
FfPt = FfSPt+FfMPt+FfLPt;
FfDt = FfSDt+FfMDt+FfLDt;
%%
figure(1)
plot(htime,log10(HFt),'r','Linewidth',2); hold on;
plot(htime,log10(HPt),'b','Linewidth',2); hold on;
plot(htime,log10(HDt),'k','Linewidth',2); hold on;
plot(ftime,log10(FpFt),'m','Linewidth',2); hold on;
plot(ftime,log10(FpPt),'c','Linewidth',2); hold on;
plot(ftime,log10(FpDt),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(ftime,log10(FfFt),'r','Linewidth',2); hold on;
plot(ftime,log10(FfPt),'b','Linewidth',2); hold on;
plot(ftime,log10(FfDt),'k','Linewidth',2); hold on;
xlim([1860 2100])
%ylim([-0.4 0.6])
xlabel('Year')
ylabel('log_1_0 Biomass (g m^-^2)')
%title(['Forecast fished'])
print('-dpng',[ppath 'Timeseries_GoMex_all_types.png'])

%% since 1980
figure(2)
plot(htime,log10(HFt),'r','Linewidth',2); hold on;
plot(htime,log10(HPt),'b','Linewidth',2); hold on;
plot(htime,log10(HDt),'k','Linewidth',2); hold on;
plot(ftime,log10(FpFt),'m','Linewidth',2); hold on;
plot(ftime,log10(FpPt),'c','Linewidth',2); hold on;
plot(ftime,log10(FpDt),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(ftime,log10(FfFt),'r','Linewidth',2); hold on;
plot(ftime,log10(FfPt),'b','Linewidth',2); hold on;
plot(ftime,log10(FfDt),'k','Linewidth',2); hold on;
xlim([1980 2100])
xlabel('Year')
ylabel('log_1_0 Biomass (g m^-^2)')
%title(['Forecast fished'])
print('-dpng',[ppath 'Timeseries_1980_GoMex_all_types.png'])

%% since 1980
FallFt = [HFt FfFt];
FallPt = [HPt FfPt];
FallDt = [HDt FfDt];

PallFt = [HFt FpFt];
PallPt = [HPt FpPt];
PallDt = [HDt FpDt];

yr = [htime ftime];

figure(3)
plot(yr,log10(PallFt),'m','Linewidth',2); hold on;
plot(yr,log10(PallPt),'c','Linewidth',2); hold on;
plot(yr,log10(PallDt),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(yr,log10(FallFt),'r','Linewidth',2); hold on;
plot(yr,log10(FallPt),'b','Linewidth',2); hold on;
plot(yr,log10(FallDt),'k','Linewidth',2); hold on;
legend('F unfished','P unfished','D unfished','F fished','P fished','D fished')
legend('location','northeast')
xlim([1980 2100])
ylim([-0.2 0.6])
xlabel('Year')
ylabel('log_1_0 Biomass (g m^-^2)')
%title(['Forecast fished'])
print('-dpng',[ppath 'Timeseries_1980_GoMex_all_types_legend.png'])

%%
figure(4)
plot(yr,log10(FallFt),'r','Linewidth',2); hold on;
plot(yr,log10(FallPt),'b','Linewidth',2); hold on;
plot(yr,log10(FallDt),'k','Linewidth',2); hold on;
legend('F','P','D')
legend('location','northeast')
xlim([1980 2100])
ylim([-0.2 0.6])
xlabel('Year')
ylabel('log_1_0 Biomass (g m^-^2)')
%title(['Forecast fished'])
print('-dpng',[ppath 'Timeseries_fished_1980_GoMex_all_types_legend.png'])

%% Diffs

diffF = FfF-HF;
diffP = FfP-HP;
diffD = FfD-HD;
diffAll = FfAll-HAll;

pdiffF = (FfF-HF) ./ HF;
pdiffP = (FfP-HP) ./ HP;
pdiffD = (FfD-HD) ./ HD;
pdiffAll = (FfAll-HAll) ./ HAll;

%% plot info
plotminlat=15; %Set these bounds for your data
plotmaxlat=35;
plotminlon=-100;
plotmaxlon=-75;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];


%% Maps
% TOTAL CHANGE
figure(5)
subplot('Position',[0 0.5 0.48 0.45])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('F');

%
subplot('Position',[0 0 0.48 0.45])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.47 0.25 0.03 0.5],'orientation','vertical')
set(gcf,'renderer','painters')
title('D');
%
subplot('Position',[0.51 0.5 0.48 0.45])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('P');

%
subplot('Position',[0.51 0 0.48 0.45])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('All');
print('-dpng',[ppath 'Fore_fished_Hist_' harv '_global_diffAll.png'])

%% 
% PERCENT CHANGE
figure(6)
subplot('Position',[0 0.5 0.48 0.45])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('F');

%
subplot('Position',[0 0 0.48 0.45])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.47 0.25 0.03 0.5],'orientation','vertical')
set(gcf,'renderer','painters')
title('D');
%
subplot('Position',[0.51 0.5 0.48 0.45])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('P');

%
subplot('Position',[0.51 0 0.48 0.45])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('All');
print('-dpng',[ppath 'Fore_fished_Hist_' harv '_global_percent_diffAll.png'])


