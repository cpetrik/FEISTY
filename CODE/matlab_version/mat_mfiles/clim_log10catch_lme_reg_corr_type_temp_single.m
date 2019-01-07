%POEM catch vs. SAUP catch by LME

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/CSV/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'LME_clim_temp.mat']);

%use weighted landings and all catches (landings, IUU, & discards)
%units in tonnes = MT (1e3 kg)
load([spath 'Reg_LME_Catch_annual.mat']);

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')

AREA_OCN = max(area,1);

lfile = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100';
lpath = ['/Volumes/GFDL/NC/Matlab_new_size/' lfile '/'];
load([lpath 'LME_clim_fished03_' lfile '.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;

% plot info
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

%%
% sum of functional types in each
Flme_land_all = nansum(Flme_wland,3);
Plme_land_all = nansum(Plme_wland,3);
Dlme_land_all = nansum(Dlme_wland,3);
Flme_catch_all = nansum(Flme_wall,3);
Plme_catch_all = nansum(Plme_wall,3);
Dlme_catch_all = nansum(Dlme_wall,3);

% 1970-2005 average
id = find(yr>1955 & yr<=2005);

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

% Top 10 yrs
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

%% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);

%Reg does not have 65 or 66
keep = [1:9,11:34,36:38,46:54,59:60,62];

x=-6:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

frate = 0.3;
tfish = num2str(100+int64(10*frate));

cfile = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE10_CC100_lgRE00100_mdRE00100';

harv = 'All_fish03';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

%load([dpath 'LME_clim_fished',harv,'_' cfile '.mat'],'lme_mcatch');
load([dpath 'LME_clim_',harv,'_loop_' cfile '.mat'],'lme_mcatch');

close all

%% POEM LME biomass in MT
plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
% MT/km2
plme_mcatch = plme_mcatch ./ lme_area_km2;
plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;

%% Plots
l10l=log10(slme_mland10+eps);
l10lF=log10(Flme_mland10+eps);
l10lP=log10(Plme_mland10+eps);
l10lD=log10(Dlme_mland10+eps);

l10c=log10(slme_mcatch10+eps);
l10cF=log10(Flme_mcatch10+eps);
l10cP=log10(Plme_mcatch10+eps);
l10cD=log10(Dlme_mcatch10+eps);

l10p=log10(plme_mcatch);
l10pF=log10(plme_Fmcatch);
l10pP=log10(plme_Pmcatch);
l10pD=log10(plme_Dmcatch);

% Landings
figure(1)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:66
    plot(l10l(i),l10p(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
%axis([-6 2 -6 2])
xlabel('Watson mean of top 10 years')
ylabel('POEM mean of Climatology')
title('Mean landings')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:66
    plot(l10lF(i),l10pF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
%axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean F landings')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:66
    plot(l10lP(i),l10pP(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
%axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean P landings')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:66
    plot(l10lD(i),l10pD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
%axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean D landings')
stamp(cfile)
print('-dpng',[ppath 'Clim_fished',harv,'_Watson_log10land_comp_types_temp.png'])


% All catches
figure(2)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:66
    plot(l10c(i),l10p(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
%axis([-6 2 -6 2])
xlabel('Watson mean of top 10 years')
ylabel('POEM mean of Climatology')
title('Mean catches (landings + IUU + discards)')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:66
    plot(l10cF(i),l10pF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
%axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean F catches')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:66
    plot(l10cP(i),l10pP(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
%axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean P catches')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:66
    plot(l10cD(i),l10pD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
%axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean D catches')
stamp(cfile)
print('-dpng',[ppath 'Clim_fished',harv,'_Watson_log10catch_comp_types_temp.png'])



%% Drop Arctic, Antarctic, Hawaii, Australia -------------------------

%LANDINGS ---------------------------
% Stats
%r
rall=corr(l10l(keep),l10p(keep));
rF=corr(l10lF(keep),l10pF(keep));
rP=corr(l10lP(keep),l10pP(keep));
rD=corr(l10lD(keep),l10pD(keep));

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


figure(3)
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
axis([-4 2 -4 2])
xlabel('Watson mean of top 10 years')
ylabel('POEM total mean of Climatology')
title('Mean landings without Polar and Australia')

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
axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean F landings')

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
axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean P landings')

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
axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean D landings')
stamp(cfile)
print('-dpng',[ppath 'Clim_fished',harv,'_Watson_log10land_comp_types_LELC_temp.png'])


%% CATCHES ---------------------------
% Stats
%r
rall=corr(l10c(keep),l10p(keep));
rF=corr(l10cF(keep),l10pF(keep));
rP=corr(l10cP(keep),l10pP(keep));
rD=corr(l10cD(keep),l10pD(keep));

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


figure(4)
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
axis([-4 2 -4 2])
xlabel('Watson mean of top 10 years')
ylabel('POEM total mean of Climatology')
title('Mean catches (all) without Polar and Australia')

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
axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean F catches')

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
axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean P catches')

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
axis([-6 2 -6 2])
xlabel('Watson (log10 MT km^-^2)')
ylabel('POEM (log10 MT km^-^2)')
title('Mean D catches')
stamp(cfile)
print('-dpng',[ppath 'Clim_fished',harv,'_Watson_log10catch_comp_types_LELC_temp.png'])

%% For ms
% figure(5)
% subplot(2,2,4)
% plot(x,x,'--k'); hold on;
% plot(x,x2h,'--b'); hold on;
% plot(x,x2l,'--b'); hold on;
% plot(x,x5h,'--r'); hold on;
% plot(x,x5l,'--r'); hold on;
% for i=1:length(keep)
%     lme=keep(i);
%     plot(l10s(lme),l10p(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
% end
% text(-5.5,1.5,['r = ' sprintf('%2.2f',rall)])
% text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmse)])
% axis([-6 2 -6 2])
% xlabel('Watson')
% ylabel('POEM')
% title('D. All fishes')
% 
% subplot(2,2,1)
% plot(x,x,'--k'); hold on;
% plot(x,x2h,'--b'); hold on;
% plot(x,x2l,'--b'); hold on;
% plot(x,x5h,'--r'); hold on;
% plot(x,x5l,'--r'); hold on;
% for i=1:length(keep)
%     lme=keep(i);
%     plot(l10sF(lme),l10pF(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
% end
% text(-5.5,1.5,['r = ' sprintf('%2.2f',rF)])
% text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseF)])
% axis([-6 2 -6 2])
% xlabel('Watson')
% ylabel('POEM')
% title('A. Forage Fishes')
% 
% subplot(2,2,2)
% plot(x,x,'--k'); hold on;
% plot(x,x2h,'--b'); hold on;
% plot(x,x2l,'--b'); hold on;
% plot(x,x5h,'--r'); hold on;
% plot(x,x5l,'--r'); hold on;
% for i=1:length(keep)
%     lme=keep(i);
%     plot(l10sP(lme),l10pP(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
% end
% text(-5.5,1.5,['r = ' sprintf('%2.2f',rP)])
% text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseP)])
% axis([-6 2 -6 2])
% xlabel('Watson')
% ylabel('POEM')
% title('B. Large Pelagics')
% 
% subplot(2,2,3)
% plot(x,x,'--k'); hold on;
% plot(x,x2h,'--b'); hold on;
% plot(x,x2l,'--b'); hold on;
% plot(x,x5h,'--r'); hold on;
% plot(x,x5l,'--r'); hold on;
% for i=1:length(keep)
%     lme=keep(i);
%     plot(l10lD(lme),l10pD(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
% end
% text(-5.5,1.5,['r = ' sprintf('%2.2f',rD)])
% text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseD)])
% axis([-6 2 -6 2])
% xlabel('Watson')
% ylabel('POEM')
% title('C. Demersals')
% % stamp(cfile)
% print('-dpng',[ppath 'Clim_fished',harv,'_Watson_log10catch_comp_types_LELC_temp_ms.png'])
% 
% 
% 
%% ID LMEs
figure(7)
subplot(2,1,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:66
    text(l10cF(i),l10pF(i),num2str(i)); hold on;
end
axis([-6 1 -2 0.5])
% axis([-6 1 -6 1])
% axis('square')
xlabel('Watson F catch (log10 MT km^-^2)')
ylabel('POEM F catch (log10 MT km^-^2)')
title('Mean F catch all')

subplot(2,1,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    text(l10cF(lme),l10pF(lme),num2str(lme)); hold on;
end
axis([-6 1 -2 0.5])
% axis([-6 1 -6 1])
% axis('square')
xlabel('Watson F catch (log10 MT km^-^2)')
ylabel('POEM F catch (log10 MT km^-^2)')
title('Mean F catch no LELC')
print('-dpng',[ppath 'Clim_',harv,'_Watson_log10catch_comp_F_LMEs.png'])


figure(8)
subplot(2,1,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:66
    text(l10lF(i),l10pF(i),num2str(i)); hold on;
end
axis([-6 1 -2 0.5])
% axis([-6 1 -6 1])
% axis('square')
xlabel('Watson F landings (log10 MT km^-^2)')
ylabel('POEM F catch (log10 MT km^-^2)')
title('Mean F catch all')

subplot(2,1,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,'--b'); hold on;
plot(x,x2l,'--b'); hold on;
plot(x,x5h,'--r'); hold on;
plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    text(l10lF(lme),l10pF(lme),num2str(lme)); hold on;
end
axis([-6 1 -2 0.5])
% axis([-6 1 -6 1])
% axis('square')
xlabel('Watson F landings (log10 MT km^-^2)')
ylabel('POEM F catch (log10 MT km^-^2)')
title('Mean F catch no LELC')
print('-dpng',[ppath 'Clim_',harv,'_Watson_log10land_comp_F_LMEs.png'])


 
