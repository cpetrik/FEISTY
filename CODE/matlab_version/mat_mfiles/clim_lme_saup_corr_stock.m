%POEM catch vs. SAUP catch by LME
%Use same methods as Stock et al. 2017 to reduce SAUP dataset

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'LME_clim_temp_zoop_det.mat']);

%use weighted catches
load([spath 'SAUP_LME_Catch_annual.mat'],'yr','totcatch','lme_catch',...
    'Flme_wcatch','Dlme_wcatch','Plme_wcatch');

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')

AREA_OCN = max(area,1);

% POEM file info
frate = 0.3;
tfish = num2str(100+int64(10*frate));

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat']);
lme_area_km2 = lme_area * 1e-6;


%% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);
close all

load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

x=-8:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% SAUP
load([spath 'SAUP_LME_Catch_top10_Stock.mat']);
sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

%% POEM LME biomass in MT
plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
plme_Bmbio = lme_mbio(:,9) * 1e-6;
plme_Bsbio = lme_sbio(:,9) * 1e-6;
% MT/km2
plme_mcatch = plme_mcatch ./ lme_area_km2;
plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;
plme_Bmbio = plme_Bmbio ./ lme_area_km2;
plme_Bsbio = plme_Bsbio ./ lme_area_km2;

pFracPD = plme_Pmcatch ./ (plme_Pmcatch + plme_Dmcatch);

l10p=log10(plme_mcatch);
l10pF=log10(plme_Fmcatch);
l10pP=log10(plme_Pmcatch);
l10pD=log10(plme_Dmcatch);

%% on grid
tlme = lme_mask_onedeg;
pFracPD_grid = NaN*ones(180,360);
sFracPD_grid = NaN*ones(180,360);
l10sF_grid = NaN*ones(180,360);
l10pF_grid = NaN*ones(180,360);
l10sP_grid = NaN*ones(180,360);
l10pP_grid = NaN*ones(180,360);
for L=1:66
    lid = find(tlme==L);
    pFracPD_grid(lid) = pFracPD(L);
    sFracPD_grid(lid) = sFracPD(L);
%     l10sF_grid(lid) = l10sF(L);
%     l10pF_grid(lid) = l10pF(L);
%     l10sP_grid(lid) = l10sP(L);
%     l10pP_grid(lid) = l10pP(L);
end

for L=1:length(keep)
    lme=keep(L);
    lid = find(tlme==lme);
    l10sF_grid(lid) = l10sF(lme);
    l10pF_grid(lid) = l10pF(lme);
    l10sP_grid(lid) = l10sP(lme);
    l10pP_grid(lid) = l10pP(lme);
end

diffPD = pFracPD_grid - sFracPD_grid;
diffF = (l10pF_grid - l10sF_grid);
diffP = (l10pP_grid - l10sP_grid);

%% Drop Arctic, Antarctic, Hawaii, Australia -------------------------

% Table of Zoop, Det, Bent, Temp, modcatch, SAUcatch
%mean zoop (lme_az) in g(WW) from mgC/m2
%mean zoop loss (lme_azl) in g(WW) from mgC/m2/d
%mean det (lme_adet) in g(WW) from mgC/m2/d 

% MT/km2
%means of all grid cells in an LME
lme_az_mtkm2 = lme_az * 1e-6 ./ lme_area_km2;
lme_azl_mtkm2 = lme_azl * 1e-6 ./ lme_area_km2;
lme_adet_mtkm2 = lme_adet * 1e-6 ./ lme_area_km2;
%sum of all grid cells in an LME
lme_asz_mtkm2 = lme_asz * 1e-6 ./ lme_area_km2;
lme_aszl_mtkm2 = lme_aszl * 1e-6 ./ lme_area_km2;
lme_asdet_mtkm2 = lme_asdet * 1e-6 ./ lme_area_km2;

% MT/km2
tab(:,1)=keep;
tab(:,2)=lme_aszl_mtkm2(keep);
tab(:,3)=lme_asdet_mtkm2(keep);
tab(:,4)=plme_Bsbio(keep);
tab(:,5)=lme_ptemp(keep);
tab(:,6)=plme_mcatch(keep);
tab(:,7)=slme_mcatch10(keep);

% gC/m2
%MTWW -> gWW = 1e6
%gWW -> gC = 1/9
%1/km2 -> 1/m2 = 1e-6
tab2 = tab;
tab2(:,2:4) = tab(:,2:4)*(1/9);
tab2(:,6:7) = tab(:,6:7)*(1/9);

T1 = array2table(tab,'VariableNames',{'LME','ZP','Det','Bent','T','Pcatch','Scatch'});
T2 = array2table(tab2,'VariableNames',{'LME','ZP','Det','Bent','T','Pcatch','Scatch'});

writetable(T1,[dpath 'LME_prod_catch_SAUP_mtkm2_' cfile '.csv'],'Delimiter',',')
writetable(T2,[dpath 'LME_prod_catch_SAUP_gCm2_' cfile '.csv'],'Delimiter',',')
save([dpath 'LME_prod_catch_SAUP_' cfile '.mat'],'tab','tab2')

%% Stats
%r
rall=corr(l10s(keep),l10p(keep));
rF=corr(l10sF(keep),l10pF(keep));
rP=corr(l10sP(keep),l10pP(keep));
rD=corr(l10sD(keep),l10pD(keep));
rPD=corr(sFracPD(keep),pFracPD(keep));

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

%Fmed
Fall=10^(median(l10s(keep)-l10p(keep)));
FF=10^(median(l10sF(keep)-l10pF(keep)));
FP=10^(median(l10sP(keep)-l10pP(keep)));
FD=10^(median(l10sD(keep)-l10pD(keep)));
FPD=10^(median(sFracPD(keep)-pFracPD(keep)));

% Table
fish_stat(1,1) = rall;
fish_stat(2,1) = rF;
fish_stat(3,1) = rP;
fish_stat(4,1) = rD;
fish_stat(1,2) = rmse;
fish_stat(2,2) = rmseF;
fish_stat(3,2) = rmseP;
fish_stat(4,2) = rmseD;
fish_stat(1,3) = Fall;
fish_stat(2,3) = FF;
fish_stat(3,3) = FP;
fish_stat(4,3) = FD;
fish_stat(5,1) = rPD;
fish_stat(5,2) = rmsePD;
fish_stat(5,3) = FPD;

Fstat = array2table(fish_stat,'VariableNames',{'r','RMSE','Fmed'},...
    'RowNames',{'All','F','P','D','FracP'});
writetable(Fstat,[dpath 'LME_SAUP_stats_' cfile '.csv'],'Delimiter',',',...
    'WriteRowNames',true)
save([dpath 'LME_SAUP_stats_' cfile '.mat'],'fish_stat')

%% Plots
figure(1)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10s(lme),l10p(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-3.5,1.5,['r = ' num2str(rall)])
text(-3.5,1.0,['RMSE = ' num2str(rmse)])
text(-3.5,0.5,['Fmed = ' num2str(Fall)])
axis([-4 2 -4 2])
xlabel('SAUP mean of top 10 years')
ylabel('POEM total mean of Climatology')
title('Mean catch without Polar and Australia')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10sF(lme),l10pF(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' num2str(rF)])
text(-5.5,1.0,['RMSE = ' num2str(rmseF)])
text(-5.5,0.5,['Fmed = ' num2str(FF)])
axis([-6 2 -6 2])
xlabel('SAUP F catch (log10 MT km^-^2)')
ylabel('POEM F catch (log10 MT km^-^2)')
title('Mean F catch')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10sP(lme),l10pP(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' num2str(rP)])
text(-5.5,1.0,['RMSE = ' num2str(rmseP)])
text(-5.5,0.5,['Fmed = ' num2str(FP)])
axis([-6 2 -6 2])
xlabel('SAUP P catch (log10 MT km^-^2)')
ylabel('POEM P catch (log10 MT km^-^2)')
title('Mean P catch')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10sD(lme),l10pD(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' num2str(rD)])
text(-5.5,1.0,['RMSE = ' num2str(rmseD)])
text(-5.5,0.5,['Fmed = ' num2str(FD)])
axis([-6 2 -6 2])
xlabel('SAUP D catch (log10 MT km^-^2)')
ylabel('POEM D catch (log10 MT km^-^2)')
title('Mean D catch')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_SAUP_comp_types_temp_Stock_LELC.png'])

%% For ms
figure(2)
subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10s(lme),l10p(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-1.75,1.7,['r = ' sprintf('%2.2f',rall)])
text(-1.75,1.4,['RMSE = ' sprintf('%2.2f',rmse)])
axis([-2 2 -2 2])
xlabel('SAUP')
ylabel('POEM')
title('D. All fishes')

subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10sF(lme),l10pF(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' sprintf('%2.2f',rF)])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseF)])
axis([-6 2 -6 2])
xlabel('SAUP')
ylabel('POEM')
title('A. Forage Fishes')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10sP(lme),l10pP(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' sprintf('%2.2f',rP)])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseP)])
axis([-6 2 -6 2])
xlabel('SAUP')
ylabel('POEM')
title('B. Large Pelagics')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10sD(lme),l10pD(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-1.75,1.7,['r = ' sprintf('%2.2f',rD)])
text(-1.75,1.4,['RMSE = ' sprintf('%2.2f',rmseD)])
axis([-2 2 -2 2])
xlabel('SAUP')
ylabel('POEM')
title('C. Demersals')
% stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_SAUP_comp_types_temp_Stock_LELC_ms.png'])

%% For ppt
% figure(3)
% subplot(2,2,1)
% plot(x,x,'--k'); hold on;
% % plot(x,x2h,'--b'); hold on;
% % plot(x,x2l,'--b'); hold on;
% % plot(x,x5h,'--r'); hold on;
% % plot(x,x5l,'--r'); hold on;
% for i=1:length(keep)
%     lme=keep(i);
%     plot(l10s(lme),l10p(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
% end
% text(-1.75,1.75,['r = ' sprintf('%2.2f',rall)])
% text(-1.75,1.50,['RMSE = ' sprintf('%2.2f',rmse)])
% axis([-2 2 -2 2])
% xlabel('SAU')
% ylabel('POEM')
% title('All fishes')
% 
% subplot(2,2,2)
% plot(x,x,'--k'); hold on;
% % plot(x,x2h,'--b'); hold on;
% % plot(x,x2l,'--b'); hold on;
% % plot(x,x5h,'--r'); hold on;
% % plot(x,x5l,'--r'); hold on;
% for i=1:length(keep)
%     lme=keep(i);
%     plot(l10sF(lme),l10pF(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
% end
% text(-5.5,1.5,['r = ' sprintf('%2.2f',rF)])
% text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseF)])
% axis([-6 2 -6 2])
% xlabel('SAU')
% ylabel('POEM')
% title('Forage Fishes')
% 
% subplot(2,2,3)
% plot(x,x,'--k'); hold on;
% % plot(x,x2h,'--b'); hold on;
% % plot(x,x2l,'--b'); hold on;
% % plot(x,x5h,'--r'); hold on;
% % plot(x,x5l,'--r'); hold on;
% for i=1:length(keep)
%     lme=keep(i);
%     plot(l10sP(lme),l10pP(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
% end
% text(-5.5,1.5,['r = ' sprintf('%2.2f',rP)])
% text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseP)])
% axis([-6 2 -6 2])
% xlabel('SAU')
% ylabel('POEM')
% title('Large Pelagics')
% 
% subplot(2,2,4)
% plot(x,x,'--k'); hold on;
% % plot(x,x2h,'--b'); hold on;
% % plot(x,x2l,'--b'); hold on;
% % plot(x,x5h,'--r'); hold on;
% % plot(x,x5l,'--r'); hold on;
% for i=1:length(keep)
%     lme=keep(i);
%     plot(l10sD(lme),l10pD(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
% end
% text(-1.75,1.75,['r = ' sprintf('%2.2f',rD)])
% text(-1.75,1.50,['RMSE = ' sprintf('%2.2f',rmseD)])
% axis([-2 2 -2 2])
% xlabel('SAU')
% ylabel('POEM')
% title('Demersals')
% % stamp([harv '_' cfile])
% print('-dpng',[ppath 'Clim_',harv,'_SAUP_comp_types_temp_Stock_LELC_ppt.png'])
% 
% %
% rlelc=corr(l10s,l10p);
% figure(4)
% plot(x,x,'--k'); hold on;
% % plot(x,x2h,'--b'); hold on;
% % plot(x,x2l,'--b'); hold on;
% % plot(x,x5h,'--r'); hold on;
% % plot(x,x5l,'--r'); hold on;
% for i=1:66
%     lme=i;
%     plot(l10s(lme),l10p(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
% end
% axis([-3.5 1.5 -3.5 1.5])
% xlabel('SAU')
% ylabel('POEM')
% title('All fishes')
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Clim_',harv,'_SAUP_comp_temp_Stock_ppt.png'])
% 
% figure(5)
% plot(x,x,'--k'); hold on;
% % plot(x,x2h,'--b'); hold on;
% % plot(x,x2l,'--b'); hold on;
% % plot(x,x5h,'--r'); hold on;
% % plot(x,x5l,'--r'); hold on;
% for i=1:length(keep)
%     lme=keep(i);
%     plot(l10s(lme),l10p(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
% end
% text(-5.5,1.5,['r = ' sprintf('%2.2f',rall)])
% axis([-3.5 1.5 -3.5 1.5])
% xlabel('SAU')
% ylabel('POEM')
% title('All fishes')
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Clim_',harv,'_SAUP_comp_temp_Stock_LELC_ppt.png'])

%% Just P
%Corr
figure(6)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10sP(lme),l10pP(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sP(lme),l10pP(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center'); hold on;
end
text(-5.5,1.5,['r = ' sprintf('%2.2f',rP)])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseP)])
axis([-7 2 -7 2])
xlabel('SAUP')
ylabel('POEM')
title('Large Pelagics')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_SAUP_comp_Stock_LELC_Pcorr.png'])

%Map
figure(7)
% Diff
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.4 1.4]);
colorbar('Ticks',[-1.4 -0.7 -0.3 0 0.3 0.7 1.4])
set(gcf,'renderer','painters')
title('POEM - SAU P difference')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_SAUP_comp_Stock_LELC_Pmap.png'])

%% Just F
%Corr
figure(8)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10sF(lme),l10pF(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sF(lme),l10pF(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center'); hold on;
end
text(-5.5,1.5,['r = ' sprintf('%2.2f',rF)])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseF)])
axis([-7 2 -7 2])
xlabel('SAUP')
ylabel('POEM')
title('Forage fishes')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_SAUP_comp_Stock_LELC_Fcorr.png'])

%Map
figure(9)
% Diff
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.4 1.4]);
colorbar('Ticks',[-1.4 -0.7 -0.3 0 0.3 0.7 1.4])
set(gcf,'renderer','painters')
title('POEM - SAU F difference')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_SAUP_comp_Stock_LELC_Fmap.png'])

%% FRACTION LARGE PELAGIC ----------------------------------
% Correlation
% figure(10)
% plot(x,x,'--k');hold on;
% for i=1:length(keep)
%     lme=keep(i);
%     plot(sFracPD(lme),pFracPD(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
%     text(sFracPD(lme),pFracPD(lme),num2str(lme),...
%         'Color','k','HorizontalAlignment','center'); hold on;
% end
% text(0.05,0.95,['r = ' sprintf('%2.2f',rPD)])
% text(0.05,0.9,['RMSE = ' sprintf('%2.2f',rmsePD)])
% text(0.05,0.85,['Fmed = ' sprintf('%2.2f',FPD)])
% axis([0 1 0 1])
% xlabel('SAUP')
% ylabel('POEM')
% title('Fraction Large Pelagics')
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Clim_',harv,'_SAUP_comp_fracP.png'])

%% Catch map
figure(11)
% POEM
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,pFracPD_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology Fraction LP catch')

% Diff
subplot('Position',[0.25 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('POEM - SAU difference')

% SAUP
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,sFracPD_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('SAU Fraction LP catch')
stamp(cfile)
print('-dpng',[ppath 'Clim_' harv '_LME_fracPD_catch_SAUP_comp.png'])

%% Subplot with maps and corr
figure(12)
% POEM
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,pFracPD_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.01 0.5 0.5 0.05],'orientation','horizontal')                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology Fraction LP catch')

% Diff
subplot('Position',[0.0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('POEM - SAU difference')

% SAUP
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,sFracPD_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('SAU Fraction LP catch')

%Corr
subplot('Position',[0.55 0.075 0.4 0.4])
plot(x,x,'--k');hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(sFracPD(lme),pFracPD(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.75,0.65,['r = ' sprintf('%2.2f',rPD)])
text(0.75,0.6,['RMSE = ' sprintf('%2.2f',rmsePD)])
text(0.75,0.55,['Fmed = ' sprintf('%2.2f',FPD)])
axis([0 1.05 0 1.05])
xlabel('SAUP')
ylabel('POEM')
title('Fraction Large Pelagics')
stamp(cfile)
print('-dpng',[ppath 'Clim_' harv '_LME_fracPD_catch_SAUP_comp_subplot.png'])


