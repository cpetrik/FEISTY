% Visualize time series output of FEISTY ensembles
% Historic (1860-2005) and Forecast time period (2006-2100) at all locations
% Saved as mat files
% Ensemble mid6, temp3
% Production instead of biomass, annual means

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

%% Original parameters
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath = ['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
%epath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

%% Historic
% load([fpath 'Means_Historic_',harv,'_prod_' cfile '.mat'],...
%     'sf_tamean','sp_tamean','sd_tamean',...
%     'mf_tamean','mp_tamean','md_tamean',...
%     'lp_tamean','ld_tamean','b_tamean');
% 
% HF_tamean = sf_tamean + mf_tamean;
% HP_tamean = sp_tamean + mp_tamean + lp_tamean;
% HD_tamean = sd_tamean + md_tamean + ld_tamean;
% HB_tamean = b_tamean;
% HA_tamean = HF_tamean + HP_tamean + HD_tamean;
% 
% clear sf_tamean sp_tamean sd_tamean mf_tamean mp_tamean md_tamean lp_tamean ld_tamean b_tamean

%% Forecast
% load([fpath 'Means_fore_',harv,'_' cfile '.mat'],...
%     'sf_tamean','sp_tamean','sd_tamean',...
%     'mf_tamean','mp_tamean','md_tamean',...
%     'lp_tamean','ld_tamean','b_tamean');
% 
% FF_tamean = sf_tamean + mf_tamean;
% FP_tamean = sp_tamean + mp_tamean + lp_tamean;
% FD_tamean = sd_tamean + md_tamean + ld_tamean;
% FB_tamean = b_tamean;
% FA_tamean = FF_tamean + FP_tamean + FD_tamean;
% 
% clear sf_tamean sp_tamean sd_tamean mf_tamean mp_tamean md_tamean lp_tamean ld_tamean b_tamean
% 
% save([fpath 'Time_Means_Historic_Forecast_',harv,'_' cfile '.mat'],...
%     'HF_tamean','HP_tamean','HD_tamean','HB_tamean',...
%     'FF_tamean','FP_tamean','FD_tamean','FB_tamean',...
%     'HA_tamean','FA_tamean','-append');
load([fpath 'ESM2M_Hist_Fore/Time_Means_Historic_Forecast_',harv,'_' cfile '.mat'],...
    'HF_tamean','HP_tamean','HD_tamean','HB_tamean',...
    'FF_tamean','FP_tamean','FD_tamean','FB_tamean',...
    'HA_tamean','FA_tamean');

%% Ensemble parameter sets
% epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
%     'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
epath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
load([epath 'Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'haTsF','haTsP','haTsD','haTmF','haTmP','haTmD','haTB','haTlP','haTlD');
load([epath 'Forecast_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'faTsF','faTsP','faTsD','faTmF','faTmP','faTmD','faTB','faTlP','faTlD');

%aT=tamean=nansum(prod.*area,1);

%% ts
%In original saved file
y1 = 1860+(1/12):(1/12):2005;
y2 = 2005+(1/12):(1/12):2100;
y = [y1 y2];

HF = haTsF + haTmF;
HP = haTsP + haTmP + haTlP;
HD = haTsD + haTmD + haTlD;
HA = HF + HP + HD;

FF = faTsF + faTmF;
FP = faTsP + faTmP + faTlP;
FD = faTsD + faTmD + faTlD;
FA = FF + FP + FD;

tForig = [HF_tamean FF_tamean];
tPorig = [HP_tamean FP_tamean];
tDorig = [HD_tamean FD_tamean];
tAorig = [HA_tamean FA_tamean];

%% Prod in g/d --> g/yr?
tF = [HF FF];
tP = [HP FP];
tD = [HD FD];
tA = [HA FA];

tF = [tF; tForig];
tP = [tP; tPorig];
tD = [tD; tDorig];
tA = [tA; tAorig];

%% Calc annual means before difference
st = 1:12:length(y);
en = 12:12:length(y);

mtF = NaN*ones(44,length(st));
mtP = mtF;
mtD = mtF;
mtA = mtF;
for i=1:length(st)
    mtF(:,i) = mean(tF(:,st(i):en(i)),2);
    mtP(:,i) = mean(tP(:,st(i):en(i)),2);
    mtD(:,i) = mean(tD(:,st(i):en(i)),2);
    mtA(:,i) = mean(tA(:,st(i):en(i)),2);
end

%% percent difference from 1951
yr = 1860 + [1:length(st)];
test=find(yr==1951);
yid=test(1);

mmdF = (mtF - mtF(:,yid)) ./ mtF(:,yid);
mmdP = (mtP - mtP(:,yid)) ./ mtP(:,yid);
mmdD = (mtD - mtD(:,yid)) ./ mtD(:,yid);
mmdA = (mtA - mtA(:,yid)) ./ mtA(:,yid);

hyr = (yr>=1951 & yr<=2000);
fyr = (yr>=2051 & yr<=2100);

vt(1,1) = mean(var(mmdF(:,hyr),0,2)); %0.0010         
vt(1,2) = mean(var(mmdF(:,fyr),0,2)); %0.0004 F decr
vt(2,1) = mean(var(mmdP(:,hyr),0,2)); %0.0029
vt(2,2) = mean(var(mmdP(:,fyr),0,2)); %0.0024 P decr
vt(3,1) = mean(var(mmdD(:,hyr),0,2)); %0.0003
vt(3,2) = mean(var(mmdD(:,fyr),0,2))  %0.0009 D incr

%% CONE OF UNCERTAINTY ann means pdiff
mF = mean(mmdF);
mP = mean(mmdP);
mD = mean(mmdD);
mA = mean(mmdA);

sF = std(mmdF);
sP = std(mmdP);
sD = std(mmdD);
sA = std(mmdA);

%create continuous x value array for plotting
X=[yr fliplr(yr)]; 
%create y values for out and then back
%+/- 1 stdev
Sa=[mA+sA fliplr(mA-sA)]; 
Sf=[mF+sF fliplr(mF-sF)]; 
Sp=[mP+sP fliplr(mP-sP)]; 
Sd=[mD+sD fliplr(mD-sD)]; 


%%
% set(groot,'defaultFontName','TimesNewRoman');
% set(groot,'defaultFontSize',16);


%% types - percent diff
figure(1)
fill(X,(Sf),'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Sp),'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Sd),'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(yr,(mF),'r','LineWidth',2); hold on;
plot(yr,(mP),'b','LineWidth',2); hold on;
plot(yr,(mD),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([1951 2100])
%ylim([-0.2 0.3])
%title('All functional types')
xlabel('Year')
ylabel('Production (g d^-^1) relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_pdiff1951_prod_annmean_types_ensem_mid6_temp3_cone_1std_yr.png'])


%% save for multipanel plot w/biom changes
dpath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
save([epath 'Hist_Fore_All_fish03_ensem6_mid_temp3_ts_prod_annmean_pdiff.mat'],...
    'Sf','Sp','Sd','Sa','mF','mP','mD','mA','mtF','mtP','mtD','mtA',...
    'X','yr');
save([dpath 'Hist_Fore_All_fish03_ensem6_mid_temp3_ts_prod_annmean_pdiff.mat'],...
    'Sf','Sp','Sd','Sa','mF','mP','mD','mA','mtF','mtP','mtD','mtA',...
    'X','yr');



