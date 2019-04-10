%POEM catch vs. SAUP catch by LME
%Use same methods as Stock et al. 2017 to reduce SAUP dataset

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/Global_Fishing_Watch/fishing_effort/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'LME_clim_temp_zoop_det.mat']);

%GFW
load([spath 'LME_GFW_pseine_hrs_all_yrs_1degree.mat']);

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

cmYOR=cbrewer('seq','YlOrRd',28);
cmRP=cbrewer('seq','RdPu',28);
cmPR=cbrewer('seq','PuRd',28);

% Assign a color to each LME based on temp
tmap=colormap(jet(66));
%tmap=cmocean('thermal',66);
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);
close all

load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

%% POEM LME biomass in MT
plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
% MT/km2
plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;

l10pF=log10(plme_Fmcatch);
%
l10v  = log10(lme_sum(:,1)+eps);
l10vk = log10(lme_sum(:,3)+eps);
l10f  = log10(lme_sum(:,2)+eps);
l10fk = log10(lme_sum(:,4)+eps);
ml10v  = log10(lme_mean(:,1)+eps);
ml10vk = log10(lme_mean(:,3)+eps);
ml10f  = log10(lme_mean(:,2)+eps);
ml10fk = log10(lme_mean(:,4)+eps);

%% ----------------------- Total of all years --------------------------
% All measures
figure(1)
subplot(2,2,1)
for i=1:length(keep)
    lme=keep(i);
    plot(l10v(lme),l10pF(lme),'.k','MarkerSize',20,'color',tmap(tid(lme,2),:)); hold on;
end
axis([1 7 -2.5 0.5])
title('Vessel hrs')
ylabel('POEM')
xlabel('GFW')

subplot(2,2,2)
for i=1:length(keep)
    lme=keep(i);
    plot(l10vk(lme),l10pF(lme),'.k','MarkerSize',20,'color',tmap(tid(lme,2),:)); hold on;
end
axis([-5 2 -2.5 0.5])
title('Vessel hrs km^-^2')
ylabel('POEM')
xlabel('GFW')

subplot(2,2,3)
for i=1:length(keep)
    lme=keep(i);
    plot(l10f(lme),l10pF(lme),'.k','MarkerSize',20,'color',tmap(tid(lme,2),:)); hold on;
end
axis([0 7 -2.5 0.5])
title('Fishing hrs')
ylabel('POEM')
xlabel('GFW')

subplot(2,2,4)
for i=1:length(keep)
    lme=keep(i);
    plot(l10fk(lme),l10pF(lme),'.k','MarkerSize',20,'color',tmap(tid(lme,2),:)); hold on;
end
axis([-6 1 -2.5 0.5])
title('Fishing hrs km^-^2')
ylabel('POEM')
xlabel('GFW')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_GFW_tot_Fcomp_temp_Stock_LELC.png'])

%% ID LMEs
figure(2)
for i=1:length(keep)
    lme=keep(i);
    plot(l10fk(lme),l10pF(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10fk(lme),l10pF(lme),num2str(lme),...
         'Color','k','HorizontalAlignment','center'); hold on;
end
axis([-6 1 -2.5 0.5])
title('Total Fishing hrs km^-^2')
ylabel('POEM')
xlabel('GFW')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_GFW_tot_fish_hrskm2_Fcomp_temp_Stock_LELC_ID.png'])

%% Linear regression
keep2 = [keep(1:15);keep(17:41);keep(43:45)];
x=l10fk(keep2);
y=l10pF(keep2);
[pfit,S] = polyfit(x,y,1);
X = -6.5:0.25:1;
[Y,delta] = polyconf(pfit,X,S);

%
figure(3)
for i=1:length(keep2)
    lme=keep2(i);
    plot(l10fk(lme),l10pF(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10fk(lme),l10pF(lme),num2str(lme),...
         'Color','k','HorizontalAlignment','center'); hold on;
end
plot(X,Y,'k'); hold on
plot(X,Y+delta,'--k'); hold on
plot(X,Y-delta,'--k'); hold on
%axis([-6 1 -2.5 0.5])

%% ----------------------- Mean of all years --------------------------
% All measures
figure(4)
subplot(2,2,1)
for i=1:length(keep)
    lme=keep(i);
    plot(ml10v(lme),l10pF(lme),'.k','MarkerSize',20,'color',tmap(tid(lme,2),:)); hold on;
end
axis([1 7 -2.5 0.5])
title('Vessel hrs')
ylabel('POEM')
xlabel('GFW')

subplot(2,2,2)
for i=1:length(keep)
    lme=keep(i);
    plot(ml10vk(lme),l10pF(lme),'.k','MarkerSize',20,'color',tmap(tid(lme,2),:)); hold on;
end
axis([-5 2 -2.5 0.5])
title('Vessel hrs km^-^2')
ylabel('POEM')
xlabel('GFW')

subplot(2,2,3)
for i=1:length(keep)
    lme=keep(i);
    plot(ml10f(lme),l10pF(lme),'.k','MarkerSize',20,'color',tmap(tid(lme,2),:)); hold on;
end
axis([0 7 -2.5 0.5])
title('Fishing hrs')
ylabel('POEM')
xlabel('GFW')

subplot(2,2,4)
for i=1:length(keep)
    lme=keep(i);
    plot(ml10fk(lme),l10pF(lme),'.k','MarkerSize',20,'color',tmap(tid(lme,2),:)); hold on;
end
axis([-6 1 -2.5 0.5])
title('Fishing hrs km^-^2')
ylabel('POEM')
xlabel('GFW')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_GFW_mean_Fcomp_temp_Stock_LELC.png'])

%% ID LMEs
figure(5)
for i=1:length(keep)
    lme=keep(i);
    plot(ml10fk(lme),l10pF(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(ml10fk(lme),l10pF(lme),num2str(lme),...
         'Color','k','HorizontalAlignment','center'); hold on;
end
axis([-6 1 -2.5 0.5])
title('Mean Fishing hrs km^-^2')
ylabel('POEM')
xlabel('GFW')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_GFW_mean_fish_hrskm2_Fcomp_temp_Stock_LELC_ID.png'])

%% Linear regression
keep2 = [keep(1:15);keep(17:41);keep(43:45)];
x=ml10fk(keep2);
y=l10pF(keep2);
[pfit,S] = polyfit(x,y,1);
X = -6.5:0.25:1;
[Y,delta] = polyconf(pfit,X,S);

%
figure(6)
for i=1:length(keep2)
    lme=keep2(i);
    plot(ml10fk(lme),l10pF(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(ml10fk(lme),l10pF(lme),num2str(lme),...
         'Color','k','HorizontalAlignment','center'); hold on;
end
plot(X,Y,'k'); hold on
plot(X,Y+delta,'--k'); hold on
plot(X,Y-delta,'--k'); hold on
%axis([-6 1 -2.5 0.5])

%% Scale both from 0 to 1
eff = (x - min(x)) ./ (max(x) - min(x));
cat = (y - min(y)) ./ (max(y) - min(y));

z=-8:0.5:8;
z2h = z+log10(2);
z2l = z-log10(2);
z5h = z+log10(5);
z5l = z-log10(5);

figure(7)
for i=1:length(keep2)
    lme=keep2(i);
    plot(eff(i),cat(i),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(eff(i),cat(i),num2str(lme),...
         'Color','k','HorizontalAlignment','center'); hold on;
end
plot(z,z,'--k'); hold on;
% plot(z,z2h,'--b'); hold on;
% plot(z,z2l,'--b'); hold on;
% plot(z,z5h,'--r'); hold on;
% plot(z,z5l,'--r'); hold on;
axis([0 1.1 0 1.1])

