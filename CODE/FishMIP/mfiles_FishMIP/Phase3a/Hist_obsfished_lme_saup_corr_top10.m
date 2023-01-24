% FEISTY catch vs. SAUP catch by LME
% Use same methods as Stock et al. 2017 to reduce SAUP dataset
% 1961-2010 ctrlclim & obsclim
% Observed effort

clear 
close all

%%
load(['/Users/cpetrik/Dropbox/Princeton/FEISTY_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

cdir = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/grid_cobalt/';
load([cdir 'LME_clim_temp_zoop_det.mat'],'lme_ptemp');

cpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';
load([cpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat']);
load([cpath 'lme_gfdl-mom6-cobalt2_onedeg.mat'],'tlme');
[ni,nj]=size(LON);

%% FEISTY LME biomass in MT/km2
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];

mod = 'obsclim_All_fishobs_v3_';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP/Phase3a/';
ppath = [pp cfile '/OneDeg/'];
if (~isfolder(ppath))
    mkdir(ppath)
end
load([fpath 'LME_Hist_',mod,'Catch_top10.mat']);

plme_mcatch = alme_mcatch10;
plme_Fmcatch = Flme_mcatch10;
plme_Pmcatch = Plme_mcatch10;
plme_Dmcatch = Dlme_mcatch10;

pFracPD = sFracPD;
clear sFracPD

l10p=log10(plme_mcatch+eps);
l10pF=log10(plme_Fmcatch+eps);
l10pP=log10(plme_Pmcatch+eps);
l10pD=log10(plme_Dmcatch+eps);

clear Flme_mcatch10 Plme_mcatch10 Dlme_mcatch10 sFracPD

%% SAUP in MT/km2
spath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/SAUP/';
load([spath 'SAUP_LME_Catch_top10_Stock.mat'])

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

%on grid
sFracPD_grid = NaN*ones(size(tlme));
pFracPD_grid = NaN*ones(size(tlme));
for L=1:66
    lid = find(tlme==L);
    sFracPD_grid(lid) = sFracPD(L);
    pFracPD_grid(lid) = pFracPD(L);
end

%% DvD on grid
load('/Users/cpetrik/Dropbox/Princeton/FEISTY_other/DanielVD_PelDem/Colleen_modeledfish_LME.mat')
dFracPD = FracLP;
clear FracLP
dFracPD_grid = NaN*ones(size(tlme));
for L=1:63
    lid = find(tlme==L);
    dFracPD_grid(lid) = dFracPD(L);
end

%% Drop Arctic, Antarctic, Hawaii, Australia -------------------------
did2 = notLELC(notLELC<=63);

diffD = pFracPD_grid - dFracPD_grid;
diffS = pFracPD_grid - sFracPD_grid;


% Stats
%r
[rall,pall]=corr(l10s(keep),l10p(keep));
[rF,pF]=corr(l10sF(keep),l10pF(keep));
[rP,pP]=corr(l10sP(keep),l10pP(keep));
[rD,pD]=corr(l10sD(keep),l10pD(keep));
[rPD,pPD]=corr(sFracPD(keep),pFracPD(keep));
[rdPD2,pdPD2]=corr(dFracPD(did2),pFracPD(did2));

%root mean square error
o=l10s(keep);
p=l10p(keep);
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);

o=l10sF(keep);
p=l10pF(keep);
n = length(o);
num=nansum((p-o).^2);
rmseF = sqrt(num/n);

o=l10sP(keep);
p=l10pP(keep);
n = length(o);
num=nansum((p-o).^2);
rmseP = sqrt(num/n);

o=l10sD(keep);
p=l10pD(keep);
n = length(o);
num=nansum((p-o).^2);
rmseD = sqrt(num/n);

o=sFracPD(keep);
p=pFracPD(keep);
n = length(o);
num=nansum((p-o).^2);
rmsePD = sqrt(num/n);

o=dFracPD(did2);
p=pFracPD(did2);
n = length(o);
num=nansum((p-o).^2);
rmsePD2 = sqrt(num/n);

%% Plots
cmYOR=cbrewer('seq','YlOrRd',28,'PCHIP');
cmRP=cbrewer('seq','RdPu',28,'PCHIP');
cmPR=cbrewer('seq','PuRd',28,'PCHIP');

x=-8:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

% For ms
figure(1)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10sF(keep),l10pF(keep),20,lme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.375 0.5 0.3 0.025],'orientation','horizontal')
text(-5.5,1.5,['r = ' sprintf('%2.2f',rF) ' (p = ' sprintf('%2.2f',pF) ')'])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseF)])
axis([-7 2 -7 2])
xlabel('SAU')
ylabel('FEISTY')
title('Forage Fishes')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10sP(keep),l10pP(keep),20,lme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-5.5,1.5,['r = ' sprintf('%2.2f',rP) ' (p = ' sprintf('%2.2f',pP) ')'])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseP)])
axis([-7 2 -7 2])
xlabel('SAU')
ylabel('FEISTY')
title('Large Pelagics')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10sD(keep),l10pD(keep),20,lme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-3.75,1.5,['r = ' sprintf('%2.2f',rD) ' (p = ' sprintf('%2.2f',pD) ')'])
text(-3.75,1.0,['RMSE = ' sprintf('%2.2f',rmseD)])
axis([-4 2 -4 2])
xlabel('SAU')
ylabel('FEISTY')
title('Demersals')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10s(keep),l10p(keep),20,lme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-3.75,1.5,['r = ' sprintf('%2.2f',rall) ' (p = ' sprintf('%2.2f',pall) ')'])
text(-3.75,1.0,['RMSE = ' sprintf('%2.2f',rmse)])
axis([-4 2 -4 2])
xlabel('SAU')
ylabel('FEISTY')
title('All fishes')
stamp(mod)
print('-dpng',[ppath 'LME_Hist_',mod,'onedeg_Catch_top10_SAUP_comp.png'])

%% LgPel vs. Dem
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

x=0:0.1:1;

load coastlines

%% Subplot with maps and corr and pval
figure(2)
%SAU
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,diffS)
cmocean('balance')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.25 0.56 0.5 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
title('FEISTY - SAU difference')

%DvD
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,diffD)
cmocean('balance')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('FEISTY - vanD difference')

%SAU corr
subplot('Position',[0.1 0.16 0.35 0.35])
plot(x,x,'--k');hold on;
scatter(sFracPD(notLELC),pFracPD(notLELC),20,lme_ptemp(notLELC,1),'filled'); hold on;
cmocean('thermal');
text(0.725,0.55,['r = ' sprintf('%2.2f',rPD)])
text(0.725,0.49,['(p = ' sprintf('%2.2f',pPD) ')'])
text(0.725,0.41,['RMSE = ' sprintf('%2.2f',rmsePD)])
axis([0 1.05 0 1.05])
xlabel('SAU')
ylabel('FEISTY')
%title('Fraction Large Pelagics')

%DvD Corr
subplot('Position',[0.575 0.16 0.35 0.35])
plot(x,x,'--k');hold on;
scatter(dFracPD(did2),pFracPD(did2),20,lme_ptemp(did2,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.25 0.05 0.5 0.025],'orientation','horizontal')
text(0.725,0.55,['r = ' sprintf('%2.2f',rdPD2)])
text(0.725,0.49,['(p = ' sprintf('%2.2f',pdPD2) ')'])
text(0.725,0.41,['RMSE = ' sprintf('%2.2f',rmsePD2)])
axis([0 1.05 0 1.05])
xlabel('vanD')
ylabel('FEISTY')
%title('Fraction Large Pelagics')
stamp(mod)
print('-dpng',[ppath 'LME_Hist_',mod,'onedeg_fracPD_catch_SAUP_DvD_comp.png'])
