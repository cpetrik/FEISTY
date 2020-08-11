function [r,rmse,ss,mis] = hist_1meso_lme_saup_corr_loop(lme_mcatch,lme_area,ppath)

%FEISTY catch vs. SAUP catch by LME
%Use same methods as Stock et al. 2017 to reduce SAUP dataset

gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

load([gpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t');
load([gpath 'lme_mask_esm2m.mat'],'lme_mask_esm2m');
load([cpath 'LME_hist_temp_zoop_det.mat'],'lme_ptemp');
tlme = lme_mask_esm2m';
lme_area_km2 = lme_area * 1e-6;

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/SAUP_Stock_top10.mat');
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

harv = 'All_fish03';

%% Plot info
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; 

x=-8:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% SAUP
% MT/km2
sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

%% FEISTY  LME biomass in MT
plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
% MT/km2
plme_mcatch = plme_mcatch ./ lme_area_km2;
plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;

pFracPD = plme_Pmcatch ./ (plme_Pmcatch + plme_Dmcatch);

l10p=log10(plme_mcatch+eps);
l10pF=log10(plme_Fmcatch+eps);
l10pP=log10(plme_Pmcatch+eps);
l10pD=log10(plme_Dmcatch+eps);

%% Drop Arctic, Antarctic, Hawaii, Australia -------------------------
% Stats
mis = nan*ones(length(keep),5);
%r
r(1)=corr(l10s(keep),l10p(keep));
r(2)=corr(l10sF(keep),l10pF(keep));
r(3)=corr(l10sP(keep),l10pP(keep));
r(4)=corr(l10sD(keep),l10pD(keep));
r(5)=corr(sFracPD(keep),pFracPD(keep));

%root mean square error
o=l10s(keep);
p=l10p(keep);
n = length(o);
num=nansum((p-o).^2);
rmse(1) = sqrt(num/n);
ss(1) = num;
mis(:,1) = (p-o);

o=l10sF(keep);
p=l10pF(keep);
n = length(o);
num=nansum((p-o).^2);
rmse(2) = sqrt(num/n);
ss(2) = num;
mis(:,2) = (p-o);

o=l10sP(keep);
p=l10pP(keep);
n = length(o);
num=nansum((p-o).^2);
rmse(3) = sqrt(num/n);
ss(3) = num;
mis(:,3) = (p-o);

o=l10sD(keep);
p=l10pD(keep);
n = length(o);
num=nansum((p-o).^2);
rmse(4) = sqrt(num/n);
ss(4) = num;
mis(:,4) = (p-o);

o=sFracPD(keep);
p=pFracPD(keep);
n = length(o);
num=nansum((p-o).^2);
rmse(5) = sqrt(num/n);
ss(5) = num;
mis(:,5) = (p-o);

%% P:D ratios
% DvD on grid
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/DanielVD_PelDem/Colleen_modeledfish_LME.mat','FracLP')
dlme_Pfrac = NaN*ones(360,200);
for L=1:63
    lid = find(tlme==L);
    dlme_Pfrac(lid) = FracLP(L);
end

% SAUP & FEISTY on grid
sFracPD_grid = NaN*ones(360,200);
rPD_catch = NaN*ones(360,200);
for L=1:66
    lid = find(tlme==L);
    sFracPD_grid(lid) = sFracPD(L);
    rPD_catch(lid) = pFracPD(L);
end

% Comparison stats
did=[1:61,63];

diffD = rPD_catch - dlme_Pfrac;
diffS = rPD_catch - sFracPD_grid;

%r
%r(5)=corr(sFracPD(keep),pFracPD(keep));
r(6)=corr(FracLP(did),pFracPD(did));

%root mean square error
o=FracLP(did);
p=pFracPD(did);
n = length(o);
num=nansum((p-o).^2);
rmse(6) = sqrt(num/n);

%% Plot of catch vs. SAUP
close all
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
text(-5.5,1.5,['r = ' sprintf('%2.2f',r(2))])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmse(2))])
axis([-6 2 -6 2])
xlabel('SAU')
ylabel('FEISTY ')
title('Forage Fishes')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10sP(keep),l10pP(keep),20,lme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-5.5,1.5,['r = ' sprintf('%2.2f',r(3))])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmse(3))])
axis([-6 3 -6 3])
xlabel('SAU')
ylabel('FEISTY ')
title('Large Pelagics')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10sD(keep),l10pD(keep),20,lme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-1.75,1.7,['r = ' sprintf('%2.2f',r(4))])
text(-1.75,1.4,['RMSE = ' sprintf('%2.2f',rmse(4))])
axis([-3 3 -3 3])
xlabel('SAU')
ylabel('FEISTY ')
title('Demersals')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10s(keep),l10p(keep),20,lme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-1.75,1.7,['r = ' sprintf('%2.2f',r(1))])
text(-1.75,1.4,['RMSE = ' sprintf('%2.2f',rmse(1))])
axis([-3 3 -3 3])
xlabel('SAU')
ylabel('FEISTY ')
title('All fishes')
print('-dpng',[ppath 'Hist_1meso_',harv,'_SAUP_comp_types_temp_Stock_LELC_ms.png'])


%% Subplot with maps and corr no Fmed
x=0:0.1:1;

figure(2)
%SAU
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffS)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.25 0.56 0.5 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
title('FEISTY - SAU difference')

%DvD
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('FEISTY - vanD difference')

%SAU corr
subplot('Position',[0.1 0.16 0.35 0.35])
plot(x,x,'--k');hold on;
scatter(sFracPD(notLELC),pFracPD(notLELC),20,lme_ptemp(notLELC,1),'filled'); hold on;
cmocean('thermal');
text(0.725,0.55,['r = ' sprintf('%2.2f',r(5))])
text(0.725,0.49,['RMSE = ' sprintf('%2.2f',rmse(5))])
axis([0 1.05 0 1.05])
xlabel('SAU')
ylabel('FEISTY')
%title('Fraction Large Pelagics')

%DvD Corr
subplot('Position',[0.575 0.16 0.35 0.35])
plot(x,x,'--k');hold on;
scatter(FracLP(did),pFracPD(did),20,lme_ptemp(did,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.25 0.05 0.5 0.025],'orientation','horizontal')
text(0.725,0.55,['r = ' sprintf('%2.2f',r(6))])
text(0.725,0.49,['RMSE = ' sprintf('%2.2f',rmse(6))])
axis([0 1.05 0 1.05])
xlabel('vanD')
ylabel('FEISTY')
%title('Fraction Large Pelagics')
print('-dpng',[ppath 'Hist_1meso_' harv '_LME_fracPD_catch_SAUP_DvD_comp_subplot.png'])



