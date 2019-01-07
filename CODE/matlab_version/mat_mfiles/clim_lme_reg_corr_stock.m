%POEM catch vs. Reg Watson catch by LME
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
load([spath 'Reg_LME_Catch_annual.mat']);

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')

AREA_OCN = max(area,1);

% "lfile" never changes, has lme areas
lfile = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100';
lpath = ['/Volumes/GFDL/NC/Matlab_new_size/' lfile '/'];
load([lpath 'LME_clim_fished03_' lfile '.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;

% POEM file info
frate = 0.3;
tfish = num2str(100+int64(10*frate));

cfile = 'Dc_enc70-b200_cm20_m-b175-k08_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat'],'lme_mcatch','lme_mbio','lme_sbio');
%load([dpath 'LME_clim_',harv,'_loop_' cfile '.mat'],'lme_mcatch');


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
%Reg does not have 65 or 66
kid = (notLELC<65);
keep = keep(kid==1);

x=-6:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% Reg
% sum of functional types in each
Flme_land_all = nansum(Flme_wland,3);
Plme_land_all = nansum(Plme_wland,3);
Dlme_land_all = nansum(Dlme_wland,3);
Flme_catch_all = nansum(Flme_wall,3);
Plme_catch_all = nansum(Plme_wall,3);
Dlme_catch_all = nansum(Dlme_wall,3);

%1950-2006 Reg average
id = find(yr>1950 & yr<=2006);

slme_mland = nanmean(lme_land(id,:));
slme_mland = slme_mland';
Fslme_mland = nanmean(Flme_land_all(id,:));
Fslme_mland = Fslme_mland';
Pslme_mland = nanmean(Plme_land_all(id,:));
Pslme_mland = Pslme_mland';
Dslme_mland = nanmean(Dlme_land_all(id,:));
Dslme_mland = Dslme_mland';

slme_mcatch = nanmean(lme_all(id,:));
slme_mcatch = slme_mcatch';
Fslme_mcatch = nanmean(Flme_catch_all(id,:));
Fslme_mcatch = Fslme_mcatch';
Pslme_mcatch = nanmean(Plme_catch_all(id,:));
Pslme_mcatch = Pslme_mcatch';
Dslme_mcatch = nanmean(Dlme_catch_all(id,:));
Dslme_mcatch = Dslme_mcatch';

%Top 10 yrs by total
slme_mland10 = NaN*ones(size(slme_mland));
Flme_mland10 = NaN*ones(size(slme_mland));
Plme_mland10 = NaN*ones(size(slme_mland));
Dlme_mland10 = NaN*ones(size(slme_mland));
slme_mcatch10 = NaN*ones(size(slme_mcatch));
Flme_mcatch10 = NaN*ones(size(slme_mcatch));
Plme_mcatch10 = NaN*ones(size(slme_mcatch));
Dlme_mcatch10 = NaN*ones(size(slme_mcatch));
for i=1:66
    [sort_lme_land,ix] = sort(lme_land(:,i),'descend');
    sort_Flme_land = Flme_land_all(ix,i);
    sort_Plme_land = Plme_land_all(ix,i);
    sort_Dlme_land = Dlme_land_all(ix,i);
    slme_mland10(i) = nanmean(sort_lme_land(1:10));
    Flme_mland10(i) = nanmean(sort_Flme_land(1:10));
    Plme_mland10(i) = nanmean(sort_Plme_land(1:10));
    Dlme_mland10(i) = nanmean(sort_Dlme_land(1:10));
    
    [sort_lme_all,ix] = sort(lme_all(:,i),'descend');
    sort_Flme_catch = Flme_catch_all(ix,i);
    sort_Plme_catch = Plme_catch_all(ix,i);
    sort_Dlme_catch = Dlme_catch_all(ix,i);
    slme_mcatch10(i) = nanmean(sort_lme_all(1:10));
    Flme_mcatch10(i) = nanmean(sort_Flme_catch(1:10));
    Plme_mcatch10(i) = nanmean(sort_Plme_catch(1:10));
    Dlme_mcatch10(i) = nanmean(sort_Dlme_catch(1:10));
end

% MT/km2
slme_mland10 = slme_mland10 ./ lme_area_km2;
Flme_mland10 = Flme_mland10 ./ lme_area_km2;
Plme_mland10 = Plme_mland10 ./ lme_area_km2;
Dlme_mland10 = Dlme_mland10 ./ lme_area_km2;

slme_mcatch10 = slme_mcatch10 ./ lme_area_km2;
Flme_mcatch10 = Flme_mcatch10 ./ lme_area_km2;
Plme_mcatch10 = Plme_mcatch10 ./ lme_area_km2;
Dlme_mcatch10 = Dlme_mcatch10 ./ lme_area_km2;

lFracPD = Plme_mland10 ./ (Plme_mland10 + Dlme_mland10);
cFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10l=log10(slme_mland10+eps);
l10lF=log10(Flme_mland10+eps);
l10lP=log10(Plme_mland10+eps);
l10lD=log10(Dlme_mland10+eps);

l10c=log10(slme_mcatch10+eps);
l10cF=log10(Flme_mcatch10+eps);
l10cP=log10(Plme_mcatch10+eps);
l10cD=log10(Dlme_mcatch10+eps);

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

% Drop Arctic, Antarctic, Hawaii, Australia -------------------------
%% LANDINGS ---------------------------
% Stats
%r
rall=corr(l10l(keep),l10p(keep));
rF=corr(l10lF(keep),l10pF(keep));
rP=corr(l10lP(keep),l10pP(keep));
rD=corr(l10lD(keep),l10pD(keep));
rPD=corr(lFracPD(keep),pFracPD(keep));

%root mean square error
o=l10l(keep);
p=l10p(keep);
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);

o=l10lF(keep);
p=l10pF(keep);
n = length(o);
num=nansum((p-o).^2);
rmseF = sqrt(num/n);

o=l10lP(keep);
p=l10pP(keep);
n = length(o);
num=nansum((p-o).^2);
rmseP = sqrt(num/n);

o=l10lD(keep);
p=l10pD(keep);
n = length(o);
num=nansum((p-o).^2);
rmseD = sqrt(num/n);

o=lFracPD(keep);
p=pFracPD(keep);
n = length(o);
num=nansum((p-o).^2);
rmsePD = sqrt(num/n);

%Fmed
Fall=10^(median(l10l(keep)-l10p(keep)));
FF=10^(median(l10lF(keep)-l10pF(keep)));
FP=10^(median(l10lP(keep)-l10pP(keep)));
FD=10^(median(l10lD(keep)-l10pD(keep)));
FPD=10^(median(lFracPD(keep)-pFracPD(keep)));

% Table
land_stat(1,1) = rall;
land_stat(2,1) = rF;
land_stat(3,1) = rP;
land_stat(4,1) = rD;
land_stat(1,2) = rmse;
land_stat(2,2) = rmseF;
land_stat(3,2) = rmseP;
land_stat(4,2) = rmseD;
land_stat(1,3) = Fall;
land_stat(2,3) = FF;
land_stat(3,3) = FP;
land_stat(4,3) = FD;
land_stat(5,1) = rPD;
land_stat(5,2) = rmsePD;
land_stat(5,3) = FPD;

Lstat = array2table(land_stat,'VariableNames',{'r','RMSE','Fmed'},...
    'RowNames',{'All','F','P','D','FracP'});
writetable(Lstat,[dpath 'LME_Reg_landings_stats_' cfile '.csv'],'Delimiter',',',...
    'WriteRowNames',true)
save([dpath 'LME_Reg_landings_stats_' cfile '.mat'],'land_stat')

% Plots
figure(1)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10l(lme),l10p(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-3.5,1.5,['r = ' num2str(rall)])
text(-3.5,1.0,['RMSE = ' num2str(rmse)])
text(-3.5,0.5,['Fmed = ' num2str(Fall)])
axis([-4 2 -4 2])
xlabel('Watson mean of top 10 years')
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
    plot(l10lF(lme),l10pF(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' num2str(rF)])
text(-5.5,1.0,['RMSE = ' num2str(rmseF)])
text(-5.5,0.5,['Fmed = ' num2str(FF)])
axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean F landings/catch')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10lP(lme),l10pP(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' num2str(rP)])
text(-5.5,1.0,['RMSE = ' num2str(rmseP)])
text(-5.5,0.5,['Fmed = ' num2str(FP)])
axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean P landings/catch')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10lD(lme),l10pD(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' num2str(rD)])
text(-5.5,1.0,['RMSE = ' num2str(rmseD)])
text(-5.5,0.5,['Fmed = ' num2str(FD)])
axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean D landings/catch')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_Watson_landings_comp_types_temp_Stock_LELC.png'])

% For ms
figure(2)
subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10l(lme),l10p(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' sprintf('%2.2f',rall)])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmse)])
axis([-6 2 -6 2])
xlabel('Watson')
ylabel('POEM')
title('D. All fishes')

subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10lF(lme),l10pF(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' sprintf('%2.2f',rF)])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseF)])
axis([-6 2 -6 2])
xlabel('Watson')
ylabel('POEM')
title('A. Forage Fishes')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10lP(lme),l10pP(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' sprintf('%2.2f',rP)])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseP)])
axis([-6 2 -6 2])
xlabel('Watson')
ylabel('POEM')
title('B. Large Pelagics')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10lD(lme),l10pD(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' sprintf('%2.2f',rD)])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseD)])
axis([-6 2 -6 2])
xlabel('Watson')
ylabel('POEM')
title('C. Demersals')
% stamp([harv '_' cfile])
%print('-dpng',[ppath 'Clim_',harv,'_Watson_landings_comp_types_temp_Stock_LELC_ms.png'])


%% CATCHES ---------------------------
% Stats
%r
rall=corr(l10c(keep),l10p(keep));
rF=corr(l10cF(keep),l10pF(keep));
rP=corr(l10cP(keep),l10pP(keep));
rD=corr(l10cD(keep),l10pD(keep));
rPD=corr(cFracPD(keep),pFracPD(keep));

%root mean square error
o=l10c(keep);
p=l10p(keep);
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);

o=l10cF(keep);
p=l10pF(keep);
n = length(o);
num=nansum((p-o).^2);
rmseF = sqrt(num/n);

o=l10cP(keep);
p=l10pP(keep);
n = length(o);
num=nansum((p-o).^2);
rmseP = sqrt(num/n);

o=l10cD(keep);
p=l10pD(keep);
n = length(o);
num=nansum((p-o).^2);
rmseD = sqrt(num/n);

o=cFracPD(keep);
p=pFracPD(keep);
n = length(o);
num=nansum((p-o).^2);
rmsePD = sqrt(num/n);

%Fmed
Fall=10^(median(l10c(keep)-l10p(keep)));
FF=10^(median(l10cF(keep)-l10pF(keep)));
FP=10^(median(l10cP(keep)-l10pP(keep)));
FD=10^(median(l10cD(keep)-l10pD(keep)));
FPD=10^(median(cFracPD(keep)-pFracPD(keep)));

% Table
catch_stat(1,1) = rall;
catch_stat(2,1) = rF;
catch_stat(3,1) = rP;
catch_stat(4,1) = rD;
catch_stat(1,2) = rmse;
catch_stat(2,2) = rmseF;
catch_stat(3,2) = rmseP;
catch_stat(4,2) = rmseD;
catch_stat(1,3) = Fall;
catch_stat(2,3) = FF;
catch_stat(3,3) = FP;
catch_stat(4,3) = FD;
catch_stat(5,1) = rPD;
catch_stat(5,2) = rmsePD;
catch_stat(5,3) = FPD;

Cstat = array2table(catch_stat,'VariableNames',{'r','RMSE','Fmed'},...
    'RowNames',{'All','F','P','D','FracP'});
writetable(Cstat,[dpath 'LME_Reg_catch_stats_' cfile '.csv'],'Delimiter',',',...
    'WriteRowNames',true)
save([dpath 'LME_Reg_catch_stats_' cfile '.mat'],'catch_stat')

% Plots
figure(3)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10c(lme),l10p(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-3.5,1.5,['r = ' num2str(rall)])
text(-3.5,1.0,['RMSE = ' num2str(rmse)])
text(-3.5,0.5,['Fmed = ' num2str(Fall)])
axis([-4 2 -4 2])
xlabel('Watson mean of top 10 years')
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
    plot(l10cF(lme),l10pF(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' num2str(rF)])
text(-5.5,1.0,['RMSE = ' num2str(rmseF)])
text(-5.5,0.5,['Fmed = ' num2str(FF)])
axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean F catch')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10cP(lme),l10pP(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' num2str(rP)])
text(-5.5,1.0,['RMSE = ' num2str(rmseP)])
text(-5.5,0.5,['Fmed = ' num2str(FP)])
axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean P catch')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10cD(lme),l10pD(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' num2str(rD)])
text(-5.5,1.0,['RMSE = ' num2str(rmseD)])
text(-5.5,0.5,['Fmed = ' num2str(FD)])
axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean D catch')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_Watson_catch_comp_types_temp_Stock_LELC.png'])

% For ms
figure(4)
subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10c(lme),l10p(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' sprintf('%2.2f',rall)])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmse)])
axis([-6 2 -6 2])
xlabel('Watson')
ylabel('POEM')
title('D. All fishes')

subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10cF(lme),l10pF(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' sprintf('%2.2f',rF)])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseF)])
axis([-6 2 -6 2])
xlabel('Watson')
ylabel('POEM')
title('A. Forage Fishes')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10cP(lme),l10pP(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' sprintf('%2.2f',rP)])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseP)])
axis([-6 2 -6 2])
xlabel('Watson')
ylabel('POEM')
title('B. Large Pelagics')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10cD(lme),l10pD(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' sprintf('%2.2f',rD)])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseD)])
axis([-6 2 -6 2])
xlabel('Watson')
ylabel('POEM')
title('C. Demersals')
% stamp([harv '_' cfile])
%print('-dpng',[ppath 'Clim_',harv,'_Watson_catch_comp_types_temp_Stock_LELC_ms.png'])


% ===================================================================
%% FRACTION LARGE PELAGIC ----------------------------------
% ===================================================================
% Correlation
figure(5)
plot(x,x,'--k');hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(lFracPD(lme),pFracPD(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(lFracPD(lme),pFracPD(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center'); hold on;
end
text(0.05,0.95,['r = ' sprintf('%2.2f',land_stat(5,1))])
text(0.05,0.9,['RMSE = ' sprintf('%2.2f',land_stat(5,2))])
text(0.05,0.85,['Fmed = ' sprintf('%2.2f',land_stat(5,3))])
axis([0 1 0 1])
xlabel('Watson landings')
ylabel('POEM')
title('Fraction Large Pelagics')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_Watson_landings_comp_fracP.png'])

figure(6)
plot(x,x,'--k');hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(cFracPD(lme),pFracPD(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(cFracPD(lme),pFracPD(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center'); hold on;
end
text(0.05,0.95,['r = ' sprintf('%2.2f',catch_stat(5,1))])
text(0.05,0.9,['RMSE = ' sprintf('%2.2f',catch_stat(5,2))])
text(0.05,0.85,['Fmed = ' sprintf('%2.2f',catch_stat(5,3))])
axis([0 1 0 1])
xlabel('Watson catch')
ylabel('POEM')
title('Fraction Large Pelagics')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_Watson_catch_comp_fracP.png'])

%% Catch map
%on grid
tlme = lme_mask_onedeg;
pFracPD_grid = NaN*ones(180,360);
lFracPD_grid = NaN*ones(180,360);
cFracPD_grid = NaN*ones(180,360);
for L=1:66
    lid = find(tlme==L);
    pFracPD_grid(lid) = pFracPD(L);
    lFracPD_grid(lid) = lFracPD(L);
    cFracPD_grid(lid) = cFracPD(L);
end


diffPl = pFracPD_grid - lFracPD_grid;
diffPc = pFracPD_grid - cFracPD_grid;

%% Subplot with maps and corr
figure(7)
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
surfm(geolat_t,geolon_t,diffPl)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('POEM - Watson difference')

% Reg
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lFracPD_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Watson Fraction LP landings')

%Corr
subplot('Position',[0.55 0.075 0.4 0.4])
plot(x,x,'--k');hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(lFracPD(lme),pFracPD(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.75,0.55,['r = ' sprintf('%2.2f',land_stat(5,1))])
text(0.75,0.5,['RMSE = ' sprintf('%2.2f',land_stat(5,2))])
text(0.75,0.45,['Fmed = ' sprintf('%2.2f',land_stat(5,3))])
axis([0 1.05 0 1.05])
xlabel('Watson')
ylabel('POEM')
title('Fraction Large Pelagics')
stamp(cfile)
print('-dpng',[ppath 'Clim_' harv '_LME_fracPD_Watson_landings_comp_subplot.png'])

figure(8)
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
surfm(geolat_t,geolon_t,diffPc)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('POEM - Watson difference')

% Reg
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,cFracPD_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Watson Fraction LP catch')

%Corr
subplot('Position',[0.55 0.075 0.4 0.4])
plot(x,x,'--k');hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(cFracPD(lme),pFracPD(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.75,0.55,['r = ' sprintf('%2.2f',catch_stat(5,1))])
text(0.75,0.5,['RMSE = ' sprintf('%2.2f',catch_stat(5,2))])
text(0.75,0.45,['Fmed = ' sprintf('%2.2f',catch_stat(5,3))])
axis([0 1.05 0 1.05])
xlabel('Watson')
ylabel('POEM')
title('Fraction Large Pelagics')
stamp(cfile)
print('-dpng',[ppath 'Clim_' harv '_LME_fracPD_Watson_catch_comp_subplot.png'])


